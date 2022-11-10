#pragma once

#include "SimplexNoise.h"
#include "helpers.h"

#include <stdlib.h>
#include <array>

enum Tile { Undefined, Land, Water };
enum Biome { Plain, Shrubland, Forest, Swamp };

struct World {

    static constexpr float k_altitude_max = 1.0f;
    static constexpr float k_altitude_min = -0.2f;
    static constexpr float k_lake_altitude = -0.01f;
    static constexpr float k_no_water = 0.0f;

    World(uint32_t w_width, uint32_t w_height) {
        data.resize(w_width * w_height, Undefined);
        altitude_buffer.resize(w_width * w_height, 0.0f);
        water_buffer.resize(w_width * w_height, 0.0f);
        width = w_width;
        height = w_height;
    }

    void init() {
        // Randomly initialize the world.
        for(uint32_t y = 0; y < height; ++y) {
            for(uint32_t x = 0; x < width; ++x) {

                // Clears the map with Land everywhere.
                data[x + y * width] = Land;

                // Computes the tile altitude (using simplex noise).
                float freq = 1.7f;
                float o0 = SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);
                freq *= 2.0f;
                float o1 = SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);
                freq *= 2.0f;
                float o2 = SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);
                freq *= 2.0f;
                float o3 = SimplexNoise::noise(freq * x / (float)width, freq * y / (float)height);

                // Transforms the altitude to the range [-1, 1].
                float altitude = (o0 + o1 * 0.5f + o2 * 0.25f + o3 * 0.125f) / 1.625f;
                altitude = (1.0f + altitude) * 0.5f;
                // Transforms the altitude to the range [k_altitude_min, k_altitude_max].
                altitude = k_altitude_min * (1.0f - altitude) + k_altitude_max * altitude;
                altitude_buffer[x + y * width] = altitude;

                // Clears the water buffer.
                water_buffer[x + y * width] = 0.0f;
            }
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

        return data[x + y * width];
    }

    Tile& at(uint32_t x, uint32_t y) {
        static Tile g_tile_undefined = Undefined;
        if(!is_inside(x, y)) {
            return g_tile_undefined;
        }

        return data[x + y * width];
    }

    float altitude_at(uint32_t x, uint32_t y) const {
        if(!is_inside(x, y)) {
            return k_altitude_max;
        }

        return altitude_buffer[x + y * width];
    }

    float water_at(uint32_t x, uint32_t y) const {
        if(!is_inside(x, y)) {
            return k_no_water;
        }

        return water_buffer[x + y * width];
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

                if(data[nx + ny * width] == T) {
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

                water_buffer[x + y * width] += 1.0f;

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

        // Lakes (and ponds).
        for(size_t i = 0; i < water_buffer.size(); ++i) {
            if(altitude_buffer[i] < k_lake_altitude) {
                water_buffer[i] += 1.0f;
            }
        }

        printf("\r100%%\ndone!\n");
    }

    std::vector<Tile> data;
    std::vector<float> altitude_buffer;
    std::vector<float> water_buffer;
    uint32_t width, height;
};
