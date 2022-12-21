#pragma once

#include "blur.h"
#include "helpers.h"
#include "randgen.h"
#include "SimplexNoise.h"

#include <stdlib.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

struct World {

    enum Tile { Soil, Rock, Water, Undefined };
    enum Biome { Plain, Shrubland, Forest, Swamp, Rocks, Underwater, Count };

    struct Area {
        Biome biome;
        Vec2 center;
    };

    static constexpr float k_altitude_max = 1.0f;
    static constexpr float k_altitude_min = 0.0f;
    static constexpr float k_lake_altitude = 0.2f;
    static constexpr float k_no_water = 0.0f;
    static constexpr float k_no_slope = 0.0f;

    World(uint32_t w_width, uint32_t w_height) {
        tiles_buffer.resize(w_width * w_height, Undefined);
        altitude_buffer.resize(w_width * w_height, 0.0f);
        altitude_water_buffer.resize(w_width * w_height, 0.0f);
        water_buffer.resize(w_width * w_height, 0.0f);
        underground_water_buffer.resize(w_width * w_height, 0.0f);
        normals.resize(w_width * w_height);
        width = w_width;
        height = w_height;
    }

    void init() {
        // Randomly initialize the world.
        for(uint32_t y = 0; y < height; ++y) {
            for(uint32_t x = 0; x < width; ++x) {

                // Clears the map with Land everywhere.
                tiles_buffer[x + y * width] = Soil;

                // Computes the tile altitude (using simplex noise).
                float freq = 1.7f;
                float o0 = SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);
                freq *= 2.0f;
                float o1 = SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);
                freq *= 2.0f;
                float o2 = SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);
                freq *= 2.0f;
                float o3 = SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);

                // Transforms the altitude to the range [0, 1].
                float altitude = (o0 + o1 * 0.5f + o2 * 0.25f + o3 * 0.125f) / 1.625f;
                altitude = (1.0f + altitude) * 0.5f;
                altitude_buffer[x + y * width] = altitude;

                // Clears the water buffer.
                water_buffer[x + y * width] = 0.0f;
            }
        }

        // Computes terrain's normals.
        for(uint32_t y = 0; y < height; ++y) {
            for(uint32_t x = 0; x < width; ++x) {

                float altitude = altitude_at(x, y);
                float dhx = altitude_at(x + 1, y) - altitude;
                float dhy = altitude_at(x, y - 1) - altitude;

                const float dx = 10.0f / width;
                const float dy = 10.0f / height;

                // normal = normalize(
                //  dx             0
                //  0       x      dy
                //  dhx            dhy
                // )
                Vec3 normal;
                normal.x = -dy * dhx;
                normal.y = -dhy * dx;
                normal.z = dx * dy;
                float norm = std::sqrtf(normal. x * normal.x + normal.y * normal.y + normal.z * normal.z);
                normal.x /= norm;
                normal.y /= norm;
                normal.z /= norm;

                normals[x + y * width] = normal;

                // if(altitude > 0.25f && normal.z < 0.707f) {
                //     tiles_buffer[x + y * width] = Rock;
                // }
            }
        }
    }

    std::string to_string(Biome biome) {
        switch(biome) {
            case Forest: return "Forest";
            case Rocks: return "Rocks";
            case Plain: return "Plain";
            case Shrubland: return "Shrubland";
            case Swamp: return "Swamp";
            case Underwater: return "Underwater";
            default: return "Undefined";
        }
    }

    bool is_inside(int32_t x, int32_t y) const {
        if(x < 0 || y < 0 || x >= width || y >= height) {
            return false;
        }

        return true;
    }

    Tile at(uint32_t x, uint32_t y) const {
        if(!is_inside(x, y)) {
            return Undefined;
        }

        return tiles_buffer[x + y * width];
    }

    Tile& at(uint32_t x, uint32_t y) {
        static Tile g_tile_undefined = Undefined;
        if(!is_inside(x, y)) {
            return g_tile_undefined;
        }

        return tiles_buffer[x + y * width];
    }

    float altitude_at(uint32_t x, uint32_t y) const {
        if(!is_inside(x, y)) {
            return k_altitude_max;
        }

        return altitude_buffer[x + y * width];
    }

    float water_altitude_at(uint32_t x, uint32_t y) const {
        if(!is_inside(x, y)) {
            return k_altitude_min;
        }

        return altitude_water_buffer[x + y * width];
    }

    float water_at(uint32_t x, uint32_t y) const {
        if(!is_inside(x, y)) {
            return k_no_water;
        }

        return water_buffer[x + y * width];
    }

    float underground_water_at(uint32_t x, uint32_t y) const {
        if(!is_inside(x, y)) {
            return k_no_water;
        }

        return underground_water_buffer[x + y * width];
    }

    float slope_at(uint32_t x, uint32_t y) const {

        if(!is_inside(x, y)) {
            return k_no_slope;
        }

        float slope = 0.0f;
        float altitude = altitude_at(x, y);
        slope += altitude_at(x - 1, y - 1) - altitude;
        slope += altitude_at(x, y - 1) - altitude;
        slope += altitude_at(x + 1, y - 1) - altitude;
        slope += altitude_at(x - 1, y) - altitude;
        slope += altitude_at(x + 1, y) - altitude;
        slope += altitude_at(x - 1, y + 1) - altitude;
        slope += altitude_at(x, y + 1) - altitude;
        slope += altitude_at(x + 1, y + 1) - altitude;
        return slope * 0.125f * 1000.0f;
    }

    const Vec3& normal_at(uint32_t x, uint32_t y) const {
        if(!is_inside(x, y)) {
            static constexpr Vec3 k_no_normal = { 0.0f, 0.0f, 0.0f };
            return k_no_normal;
        }

        return normals[x + y * width];
    }

    template<Tile T>
    uint32_t count_neighbours_of_type(uint32_t x, uint32_t y) const {

        uint32_t res = 0;

        for(int32_t i = -1; i <= 1; ++i) {
            int32_t nx = x + i;
            if(nx < 0 || nx >= width) {
                continue;
            }
            for(int32_t j = -1; j <= 1; ++j) {

                int32_t ny = y + j;
                if(ny < 0 || ny >= height) {
                    continue;
                }

                if(nx == x && ny == y) {
                    continue;
                }

                if(tiles_buffer[nx + ny * width] == T) {
                    ++res;
                }
            }
        }

        return res;
    }

    void compute_rivers_and_lakes() {

        printf("Computing rivers and lakes...\n");

        const size_t drops_count = 80;
        std::array<float, 8> gradient;
        for(size_t d = 0; d < drops_count; ++d) {
            printf("\r%lu%%", d * 100 / drops_count);
            fflush(NULL);

            int32_t x = rand() % width;
            int32_t y = rand() % height;

            while(is_inside(x, y)) {

                auto altitude = altitude_at(x, y);

                // Computes the flow direction's pdf.
                float pdf = 0.0f;
                gradient[0] = clamp(altitude - altitude_at(x - 1, y - 1), 0.0f, 1.0f);
                pdf += gradient[0];
                gradient[1] = clamp(altitude - altitude_at(x, y - 1), 0.0f, 1.0f);
                pdf += gradient[1];
                gradient[2] = clamp(altitude - altitude_at(x + 1, y - 1), 0.0f, 1.0f);
                pdf += gradient[2];
                gradient[3] = clamp(altitude - altitude_at(x - 1, y), 0.0f, 1.0f);
                pdf += gradient[3];
                gradient[4] = clamp(altitude - altitude_at(x + 1, y), 0.0f, 1.0f);
                pdf += gradient[4];
                gradient[5] = clamp(altitude - altitude_at(x - 1, y + 1), 0.0f, 1.0f);
                pdf += gradient[5];
                gradient[6] = clamp(altitude - altitude_at(x, y + 1), 0.0f, 1.0f);
                pdf += gradient[6];
                gradient[7] = clamp(altitude - altitude_at(x + 1, y + 1), 0.0f, 1.0f);
                pdf += gradient[7];

                if(pdf == 0.0f) 
                    break;

                for(size_t g = 0; g < 8; ++g) {
                    gradient[g] /= pdf;
                }

                // Selects a neighbour based on the pdf.
                float r = rand() / (float)RAND_MAX;
                float prob = 0.0f;
                size_t g = 0;
                for(; g < 8; ++g) {
                    prob += gradient[g];
                    if(r < prob) {
                        break;
                    }
                }

                if(g == 8)
                    break;

                water_buffer[x + y * width] = 1.0f;

                // Moves to the selected neighbour.
                switch(g) {
                    case 0: x -= 1; y -= 1; break;
                    case 1: y -= 1; break;
                    case 2: x += 1; y -= 1; break;
                    case 3: x -= 1; break;
                    case 4: x += 1; break;
                    case 5: x -= 1; y += 1; break;
                    case 6: y -= 1; break;
                    case 7: x += 1; y += 1; break;
                }
            }
        }

        // Sets the Water tiles.
        for(size_t i = 0; i < water_buffer.size(); ++i) {

            // Lakes (and ponds).
            if(altitude_buffer[i] < k_lake_altitude) {
                water_buffer[i] = 1.0f;
            }

            if(water_buffer[i] != k_no_water) {
                tiles_buffer[i] = Water;
            }
        }

        // Dilates the water where needed (to avoid isolated water tiles).
        auto new_tiles = tiles_buffer;
        for(uint32_t y = 0; y < height; ++y) {
            for(uint32_t x = 0; x < width; ++x) {
                new_tiles[x + y * width] = at(x, y);
                if(new_tiles[x + y * width] != Water) {
                    if(count_neighbours_of_type<Water>(x, y) >= 2) {
                        new_tiles[x + y * width] = Water;
                        water_buffer[x + y * width] = 1.0f;
                    }
                }
            }
        }
        tiles_buffer = new_tiles;

        // Process the water tiles' heights to ensure that they are never higher
        // in altitude than their non-water tiles.
        for(uint32_t y = 0; y < height; ++y) {
            for(uint32_t x = 0; x < width; ++x) {
                if(at(x, y) == Water) {
                    float ground_offset = 0.05f;
                    float min_terrain_height = altitude_at(x, y);
                    for(int32_t ny = y - 1; ny <= y + 1; ++ny) {
                        for(int32_t nx = x - 1; nx <= x + 1; ++nx) {
                            if(at(nx, ny) == Water) {
                                continue;
                            }
                            min_terrain_height = std::min(min_terrain_height, altitude_at(nx, ny));
                        }
                    }

                    altitude_buffer[x + y * width] = min_terrain_height;

                    // Altitude of the riverbed should be displaced down to be lower the river surface.
                    if (altitude_buffer[x + y * width] <= k_lake_altitude) {
                        ground_offset = 0.15f;
                    }

                    altitude_buffer[x + y * width] = std::max(min_terrain_height - ground_offset, 0.0f);

                    // Altitude of the water is a bit below the neighbour ground tiles.
                    altitude_water_buffer[x + y * width] = std::max(min_terrain_height - ground_offset * 0.1f, 0.0f);

                } else { 

                    altitude_water_buffer[x + y * width] = std::max(altitude_at(x, y) - 0.05f, 0.0f);
                }
            }
        }

        // Propagates water through the ground.
        auto in = water_buffer.data();
        auto out = underground_water_buffer.data();
        Blur::fast_gaussian_blur(in, out, width, height, 13.0f);

        printf("\r100%%\ndone!\n");
    }

    void compute_areas(size_t count) {

        size_t length = std::max<size_t>(1, sqrt(count));
        count = length * length;

        std::array<uint32_t, Biome::Count> stats;
        stats.fill(0);

        for(size_t s = 0; s < count; ++s) {

            Area area;
            area.center = RandGen::sample_2d(s, length, length, 0);
            area.center.x *= width;
            area.center.y *= height;

            float altitude = altitude_at(area.center.x, area.center.y);
            float slope = std::fabs(slope_at(area.center.x, area.center.y));
            float water = underground_water_at(area.center.x, area.center.y);

            float prob = rand() / (float)RAND_MAX;

            if(slope < 0.5f) {
                if(water > 0.1f && prob > water) {
                    area.biome = Swamp;
                } else if(prob < water * 50.0f) {
                    area.biome = Forest;
                } else if(prob < water * 200.0f) {
                    area.biome = Shrubland;
                } else {
                    area.biome = Plain;
                }
            } else {
                if(prob < altitude) {
                    area.biome = Rocks;
                } else {
                    area.biome = Plain;
                }
            }

            areas.push_back(area);

            ++stats[area.biome];
        }

        for(size_t b = 0; b < Biome::Count; ++b) {
            printf("Biome %s: %d\n", to_string((Biome)b).c_str(), stats[b]);
        }

        printf("%lu areas have been generated.\n", areas.size());
    }

    const Area& get_area(uint32_t x, uint32_t y) const {
        static Area g_null_area = { Plain, { 0.f, 0.f} };
        if(areas.empty()) {
            return g_null_area;
        }

        size_t nearest = 0;
        float best_dist = std::numeric_limits<float>::max();
        for(size_t a = 0; a < areas.size(); ++a) {
            float dist = sqrtf((x - areas[a].center.x) * (x - areas[a].center.x) + (y - areas[a].center.y) * (y - areas[a].center.y));
            if(dist < best_dist) {
                best_dist = dist;
                nearest = a;
            }
        }

        return areas[nearest];
    }

    const std::vector<float>& get_altitude_buffer() const {
        return altitude_buffer;
    }

    std::vector<Tile> tiles_buffer;
    std::vector<float> altitude_buffer;
    std::vector<float> altitude_water_buffer;
    std::vector<float> water_buffer;
    std::vector<float> underground_water_buffer;
    std::vector<Area> areas;
    std::vector<Vec3> normals;
    uint32_t width, height;
};
