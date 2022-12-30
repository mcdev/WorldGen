#pragma once

#include "blur.h"
#include "helpers.h"
#include "randgen.h"
#include "SimplexNoise.h"

#include "Eigen/Dense"

#include <stdlib.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

const float k_altitude_max = 1.0f;
const float k_altitude_min = 0.0f;
const float k_lake_altitude = 0.2f;
const float k_no_water = 0.0f;
const float k_no_slope = 0.0f;
const Vec3 k_default_current = { 0.0f, 0.0f, 0.0f };
const Vec3 k_no_gradient = { 0.0f, 0.0f, 0.0f };
const Vec3 k_default_normal = { 0.0f, 0.0f, 1.0f };

struct World {

    enum Tile { Soil, Rock, Water, Undefined };
    enum Biome { Plain, Shrubland, Forest, Swamp, Rocks, Underwater, Count };

    struct Area {
        Biome biome;
        Vec2 center;
    };

    World(int32_t w_width, int32_t w_height) {
        tiles_buffer.resize(w_width * w_height, Undefined);
        altitude_buffer.resize(w_width * w_height, 0.0f);
        altitude_water_buffer.resize(w_width * w_height, 0.0f);
        water_buffer.resize(w_width * w_height, 0.0f);
        water_current_buffer.resize(w_width * w_height * 3, 0.0f);
        underground_water_buffer.resize(w_width * w_height, 0.0f);
        normals.resize(w_width * w_height);
        width = clamp(w_width, 0, 4096);
        height = clamp(w_height, 0, 4096);
    }

    void init() {
        // Randomly initialize the world.

        for(int32_t y = 0; y < height; ++y) {
            for(int32_t x = 0; x < width; ++x) {

                // Clear the map with Land everywhere.
                tiles_buffer[x + y * width] = Soil;

                // Compute the tile altitude (using simplex noise).
                float freq = 1.7f;
                float o0 = SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);
                freq *= 2.0f;
                float o1 = SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);
                freq *= 2.0f;
                float o2 = SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);
                freq *= 2.0f;
                float o3 = SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);

                // Transform the altitude to the range [0, 1].
                float altitude = (o0 + o1 * 0.5f + o2 * 0.25f + o3 * 0.125f) / 1.625f;
                altitude = (1.0f + altitude) * 0.5f;
                altitude_buffer[x + y * width] = altitude;

                // Clear the water buffer.
                water_buffer[x + y * width] = 0.0f;

                // Set default (small) current everywhere.
                water_current_buffer[3 * (x + y * width) + 0] -= k_default_current.x;
                water_current_buffer[3 * (x + y * width) + 1] -= k_default_current.y;
                water_current_buffer[3 * (x + y * width) + 2] += k_default_current.z;
            }
        }

        // Compute terrain's normals.
        for(int32_t y = 0; y < height; ++y) {
            for(int32_t x = 0; x < width; ++x) {

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

    void clip(int32_t& x, int32_t& y) const {
        x = clamp(x, 0, width - 1);
        y = clamp(y, 0, height- 1);
    }

    Tile at(int32_t x, int32_t y) const {
        clip(x, y);
        return tiles_buffer[x + y * width];
    }

    Tile& at(int32_t x, int32_t y) {
        clip(x, y);
        return tiles_buffer[x + y * width];
    }

    float altitude_at(int32_t x, int32_t y) const {
        clip(x, y);
        return altitude_buffer[x + y * width];
    }

    float average_altitude_at(int32_t x, int32_t y, int32_t radius) const {
        float res = 0.0f;

        for (int32_t dy = -radius; dy <= radius; ++dy) {
            
            int32_t cy = y + dy;

            for (int32_t dx = -radius; dx <= radius; ++dx) {
                
                int32_t cx = x + dx;

                clip(cx, cy);
                res += altitude_at(cx, cy);
            }
        }
        return res / ((2 * radius + 1) * (2 * radius + 1));
    }

    float min_altitude_at(int32_t x, int32_t y, int32_t radius = 0) const {
        float res = k_altitude_max;

        for (int32_t dy = -radius; dy <= radius; ++dy) {
            
            int32_t cy = y + dy;
            
            for (int32_t dx = -radius; dx <= radius; ++dx) {
            
                int32_t cx = x + dx;
                clip(cx, cy);
                res = std::min(res, altitude_at(cx, cy));
            }
        }
        return res;
    }

    float max_altitude_at(int32_t x, int32_t y, int32_t radius = 0) const {
        float res = k_altitude_min;

        for (int32_t dy = -radius; dy <= radius; ++dy) {
            
            int32_t cy = y + dy;

            for (int32_t dx = -radius; dx <= radius; ++dx) {
                
                int32_t cx = x + dx;
                clip(cx, cy);
                res = std::max(res, altitude_at(cx, cy));
            }
        }
        return res;
    }

    float water_altitude_at(int32_t x, int32_t y) const {
        clip(x, y);
        return altitude_water_buffer[x + y * width];
    }

    float average_water_altitude_at(int32_t x, int32_t y, int32_t radius) const {
        float res = 0.0f;

        for (int32_t dy = -radius; dy <= radius; ++dy) {
            
            int32_t cy = y + dy;

            for (int32_t dx = -radius; dx <= radius; ++dx) {
                
                int32_t cx = x + dx;
                clip(cx, cy);
                res += water_altitude_at(cx, cy);
            }
        }
        return res / ((2 * radius + 1) * (2 * radius + 1));
    }

    float min_water_altitude_at(int32_t x, int32_t y, int32_t radius = 0) const {
        float res = k_altitude_max;

        for (int32_t dy = -radius; dy <= radius; ++dy) {
            
            int32_t cy = y + dy;

            for (int32_t dx = -radius; dx <= radius; ++dx) {
                int32_t cx = x + dx;
                clip(cx, cy);
                res = std::min(res, water_altitude_at(cx, cy));
            }
        }
        return res;
    }

    float max_water_altitude_at(int32_t x, int32_t y, int32_t radius = 0) const {
        float res = k_altitude_min;

        for (int32_t dy = -radius; dy <= radius; ++dy) {
            
            int32_t cy = y + dy;

            for (int32_t dx = -radius; dx <= radius; ++dx) {
                int32_t cx = x + dx;
                clip(cx, cy);
                res = std::max(res, water_altitude_at(cx, cy));
            }
        }
        return res;
    }

    float water_at(int32_t x, int32_t y) const {
        clip(x, y);
        return water_buffer[x + y * width];
    }

    Vec3 water_current_at(int32_t x, int32_t y) const {
        clip(x, y);
        return { water_current_buffer[3 * (x + y * width) + 0], water_current_buffer[3 * (x + y * width) + 1], water_current_buffer[3 * (x + y * width) + 2] };
    }

    float average_waterbed_altitude_at(int32_t x, int32_t y, int32_t radius) const {
        
        float res = 0.0f;
        int32_t count = 0;

        for (int32_t dy = -radius; dy <= radius; ++dy) {

            int32_t cy = y + dy;

            for (int32_t dx = -radius; dx <= radius; ++dx) {

                int32_t cx = x + dx;
                clip(cx, cy);

                if (at(cx, cy) != Water) {
                    continue;
                }

                Vec3 normal = normal_at(cx, cy);
                float slope = dot(normal, z_axis);

                //if (slope > 0.95f) {
                //    continue;
                //}

                res += average_altitude_at(cx, cy, 3);
                ++count;
            }
        }

        return count > 0 ? res / count : k_lake_altitude;
    }

    float underground_water_at(int32_t x, int32_t y) const {
        clip(x, y);
        return underground_water_buffer[x + y * width];
    }

    float slope_at(int32_t x, int32_t y) const {
        clip(x, y);

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

    std::pair<Vec3, Vec3> best_plane_from_points(const std::vector<Vec3>& c)
    {
        // Copy coordinates to  matrix in Eigen format.
        size_t num_points = c.size();

        if (num_points == 0) {
            return std::make_pair(zero3, k_default_normal);
        }

        Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic > coord(3, num_points);
        for (size_t i = 0; i < num_points; ++i) {
            coord.col(i)(0) = c[i].x;
            coord.col(i)(1) = c[i].y;
            coord.col(i)(2) = c[i].z;
        }

        // Calculate centroid.
        Vec3 centroid(coord.row(0).mean(), coord.row(1).mean(), coord.row(2).mean());

        // Subtract centroid.
        coord.row(0).array() -= centroid.x; coord.row(1).array() -= centroid.y; coord.row(2).array() -= centroid.z;

        // We only need the left-singular matrix here.
        // http://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
        auto svd = coord.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        Vec3 plane_normal;
        plane_normal.x = svd.matrixU().rightCols<1>()(0);
        plane_normal.y = svd.matrixU().rightCols<1>()(1);
        plane_normal.z = svd.matrixU().rightCols<1>()(2);
        return std::make_pair(centroid, plane_normal);
    }

    void best_plane_at(int32_t x, int32_t y, int32_t radius, Vec3& normal, Vec3& point) {

        // Collect points.
        std::vector<Vec3> points;
        for (int32_t dy = -radius; dy <= radius; dy += radius) {

            int32_t cy = y + dy;

            for (int32_t dx = -radius; dx <= radius; dx += radius) {
                int32_t cx = x + dx;
        
                //if (water_at(cx, cy) == 1.0f) {
                //    continue;
                //}
                clip(cx, cy);

                Vec3 pt;
                pt.x = dx;
                pt.y = dy;
                pt.z = altitude_at(cx, cy) * 36.0f;
                points.push_back(pt);
            }
        }
        
        std::pair<Vec3, Vec3> plane = best_plane_from_points(points);
        point = plane.first;
        normal = plane.second;
    }

    const Vec3& normal_at(int32_t x, int32_t y) const {
        clip(x, y);
        return normals[x + y * width];
    }

    template<Tile T>
    uint32_t count_neighbours_of_type(int32_t x, int32_t y) const {

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

    Vec2 water_gradient_at(int32_t x, int32_t y, int32_t radius = 1) {

        static const Vec2 dir[] = {
            { -0.707f, -0.707f },
            { +0.000f, -1.000f },
            { +0.707f, -0.707f },
            { -1.000f, +0.000f },
            { +1.000f, +0.000f },
            { -0.707f, +0.707f },
            { +0.000f, +1.000f },
            { +0.707f, +0.707f },
        };

        static const Vec2 offset[] = {
            { -1, -1},
            { +0, -1 },
            { +1, -1 },
            { -1, +0 },
            { +1, +0 },
            { -1, +1 },
            { +0, +1 },
            { +1, +1 },
        };

        // Compute the gradient vector.
        Vec2 gradient = zero2;
        auto altitude = water_altitude_at(x, y);

        for (auto d = 0; d < 8; ++d) {
            int32_t cx = x + offset[d].x * radius;
            int32_t cy = y + offset[d].y * radius;
            if (water_at(cx, cy) != 1.0f) {
                continue;
            }
            float z_diff = water_altitude_at(cx, cy) - altitude;
            gradient += dir[d] * z_diff;
        }

        return gradient;
    }

    void compute_rivers_and_lakes() {

        printf("Computing rivers and lakes...\n");

        const uint32_t drops_count = 80;
        std::array<float, 8> gradient;
        
        Vec2 dir[] = {
            { -0.707f, -0.707f },
            { +0.000f, -1.000f },
            { +0.707f, -0.707f },
            { -1.000f, +0.000f },
            { +1.000f, +0.000f },
            { -0.707f, +0.707f },
            { +0.000f, +1.000f },
            { +0.707f, +0.707f },
        };

        printf("Computing water flow...\n");
        for(uint32_t d = 0; d < drops_count; ++d) {
            printf("\r%d%%", d * 100 / drops_count);
            fflush(NULL);

            int32_t x = rand() % width;
            int32_t y = rand() % height;

            Vec3 current = k_default_current;
            while(is_inside(x, y)) {

                auto altitude = altitude_at(x, y);

                // Compute the flow direction's pdf.
                const float k_min_alt_diff = 0.001f;
                float pdf = 0.0f;
                gradient[0] = remap(altitude - altitude_at(x - 1, y - 1), k_min_alt_diff, 1.0f, 0.0f, 1.0f);
                pdf += gradient[0];
                gradient[1] = remap(altitude - altitude_at(x    , y - 1), k_min_alt_diff, 1.0f, 0.0f, 1.0f);
                pdf += gradient[1];
                gradient[2] = remap(altitude - altitude_at(x + 1, y - 1), k_min_alt_diff, 1.0f, 0.0f, 1.0f);
                pdf += gradient[2];
                gradient[3] = remap(altitude - altitude_at(x - 1, y    ), k_min_alt_diff, 1.0f, 0.0f, 1.0f);
                pdf += gradient[3];
                gradient[4] = remap(altitude - altitude_at(x + 1, y    ), k_min_alt_diff, 1.0f, 0.0f, 1.0f);
                pdf += gradient[4];
                gradient[5] = remap(altitude - altitude_at(x - 1, y + 1), k_min_alt_diff, 1.0f, 0.0f, 1.0f);
                pdf += gradient[5];
                gradient[6] = remap(altitude - altitude_at(x    , y + 1), k_min_alt_diff, 1.0f, 0.0f, 1.0f);
                pdf += gradient[6];
                gradient[7] = remap(altitude - altitude_at(x + 1, y + 1), k_min_alt_diff, 1.0f, 0.0f, 1.0f);
                pdf += gradient[7];

                for(size_t g = 0; g < 8; ++g) {
                    gradient[g] /= pdf;
                }

                // Select a neighbour based on the pdf.
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

                // Compute the selected neighbour coordinates.
                switch(g) {
                    case 0: 
                        x -= 1; 
                        y -= 1; 
                        break;
                    case 1: 
                        y -= 1; 
                        break;
                    case 2: 
                        x += 1; 
                        y -= 1; 
                        break;
                    case 3: 
                        x -= 1; 
                        break;
                    case 4: 
                        x += 1; 
                        break;
                    case 5: 
                        x -= 1; 
                        y += 1; 
                        break;
                    case 6: 
                        y += 1; 
                        break;
                    case 7: 
                        x += 1; 
                        y += 1; 
                        break;
                }
            }
        }
        printf("\r100%%\ndone!\n");

        // Set the Water tiles.
        for(size_t i = 0; i < water_buffer.size(); ++i) {

            // Lakes (and ponds).
            if(altitude_buffer[i] < k_lake_altitude) {
                water_buffer[i] = 1.0f;
            }

            if(water_buffer[i] != k_no_water) {
                tiles_buffer[i] = Water;
            }
        }

        // Dilate the water where needed (to avoid isolated water tiles).
        printf("Removing small islands...\n");
        auto new_tiles = tiles_buffer;
        for(int32_t y = 0; y < height; ++y) {
            for(int32_t x = 0; x < width; ++x) {
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

        // Erode the terrain where water is present.
        printf("Eroding the terrain...\n");
        for(int32_t y = 0; y < height; ++y) {
            for(int32_t x = 0; x < width; ++x) {
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

                    // Altitude of the riverbed should be displaced down to be lower than the river surface.
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

        // Compute the water surface.
        printf("Computing water surface...\n");
        for (int32_t y = 0; y < height; ++y) {
            for (int32_t x = 0; x < width; ++x) {
                if (at(x, y) != Water) {
                    // Skip non water tiles.
                    continue;
                }

                //altitude_water_buffer[x + y * width] = average_waterbed_altitude_at(x, y, 2);
                altitude_water_buffer[x + y * width] = std::max(k_lake_altitude, std::max(altitude_at(x, y), average_waterbed_altitude_at(x, y, 5)) + normal_at(x, y).z * 0.0f);
            }
        }

        // Compute the water current.
        printf("Computing water current...\n");
        for (int32_t y = 0; y < height; ++y) {
            for (int32_t x = 0; x < width; ++x) {
                if (water_at(x, y) != 1.0f) {
                    // Skip non water tiles.
                    continue;
                }

                // Compute the gradient at the current location.
                Vec3 gradient = water_gradient_at(x, y, 2);

                // If no gradient, try with a bigger radius.
                if (gradient.length() == 0.0f) {
                    gradient = water_gradient_at(x, y, 8);
                }

                // If still no gradient (large flat water surface case), fallback to a noise-based flow map.
                if (gradient.length() == 0.0f) {
                    float freq = 16.0f;
                    float angle = 1.57f * SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);
                    gradient = { cos(angle) * 0.02f, sin(angle) * 0.02f, 0.0f };
                }

                gradient *= 100.0f;

                water_current_buffer[3 * (x + y * width) + 0] = -gradient.x;
                water_current_buffer[3 * (x + y * width) + 1] = -gradient.y;
                water_current_buffer[3 * (x + y * width) + 2] = 1.0f;
            }
        }

        bool blur_current = true;
        if (blur_current) {
            printf("Blurring water current...\n");

            // Blur the current vectors to make them smoother.
            std::vector<float> c_in, c_out;
            c_in.resize(water_current_buffer.size() / 3);
            c_out.resize(water_current_buffer.size() / 3);

            for (auto i = 0; i < c_in.size(); ++i) {
                c_in[i] = water_current_buffer[3 * i + 0];
            }
            auto src = c_in.data();
            auto dst = c_out.data();
            Blur::fast_gaussian_blur(src, dst, width, height, 1.0f);
            for (auto i = 0; i < c_in.size(); ++i) {
                water_current_buffer[3 * i + 0] = c_out[i];
            }

            for (auto i = 0; i < c_in.size(); ++i) {
                c_in[i] = water_current_buffer[3 * i + 1];
            }
            src = c_in.data();
            dst = c_out.data();
            Blur::fast_gaussian_blur(src, dst, width, height, 1.0f);
            for (auto i = 0; i < c_in.size(); ++i) {
                water_current_buffer[3 * i + 1] = c_out[i];
            }
        }

        // Propagate water through the ground.
        auto in = water_buffer.data();
        auto out = underground_water_buffer.data();
        Blur::fast_gaussian_blur(in, out, width, height, 13.0f);
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

    const Area& get_area(int32_t x, int32_t y) const {
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
    std::vector<float> water_current_buffer;
    std::vector<float> underground_water_buffer;
    std::vector<Area> areas;
    std::vector<Vec3> normals;
    int32_t width, height;
};
