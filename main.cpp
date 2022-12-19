
#include "lodepng.h"
#include "world.h"
#include "helpers.h"

#include <cmath>
#include <stdio.h>

float lerp(float min, float max, float t) {
    if(t <= 0.0f) return min;
    else if(t >= 1.0f) return max;
    else return min * (1.0f - t) + max * t;
}

uint8_t raw_encode_altitude(float altitude, uint32_t min_height, uint32_t max_height) {
    float raw_altitude = lerp(min_height, max_height, altitude);
    raw_altitude = (raw_altitude - min_height) / (max_height - min_height);
    return 255.0f * clamp(raw_altitude, 0.0f, 1.0f);
}

void raw_write(const char *filename, const std::vector<uint8_t>& data) {
	FILE *fp = fopen(filename, "wb");
	if (fp == NULL) {
        return;
    }

	fwrite(data.data(), data.size(), 1, fp);
	fclose(fp);
}

void tga_write(const char *filename, uint32_t width, uint32_t height, uint8_t *dataBGRA, uint8_t dataChannels=4, uint8_t fileChannels=3)
{
	FILE *fp = NULL;
	// MSVC prefers fopen_s, but it's not portable
	//fp = fopen(filename, "wb");
	fp = fopen(filename, "wb");
	if (fp == NULL) return;

	// You can find details about TGA headers here: http://www.paulbourke.net/dataformats/tga/
	uint8_t header[18] = { 0,0,2,0,0,0,0,0,0,0,0,0, (uint8_t)(width%256), (uint8_t)(width/256), (uint8_t)(height%256), (uint8_t)(height/256), (uint8_t)(fileChannels*8), 0x20 };
	fwrite(&header, 18, 1, fp);

	for (uint32_t i = 0; i < width*height; i++)
	{
		for (uint32_t b = 0; b < fileChannels; b++)
		{
			fputc(dataBGRA[(i*dataChannels) + (b%dataChannels)], fp);
		}
	}
	fclose(fp);
}

int main(int argc, char** argv) {

    const uint32_t w_width = 512;
    const uint32_t w_height = 512;

    if(argc != 3) {
        printf("Usage is wgen <min_altitude> <max_altitude> in meters.\n");
        exit(0);
    }

    float min_height = atof(argv[1]);
    float max_height = atof(argv[2]);

    if(max_height < min_height) {
        printf("<min_height> must be lower than <max_height>.\n");
        exit(0);
    }

    World world(w_width, w_height);
    world.init();
    world.compute_rivers_and_lakes();
    world.compute_areas(256);

    // Translates the world's map into an image.
    Color c_land = { 0, 255, 0, 255 };
    Color c_river = { 0, 0, 255, 255 };
    Color c_rock = { 75, 75, 75, 255 };
    Color c_shrubland = { 113, 85, 60, 255 };
    Color c_forest = { 58, 77, 25, 255 };
    Color c_plain = { 89, 114, 57, 255 };
    Color c_swamp = { 20, 249, 182, 255 };
    Color c_water;

    std::vector<uint8_t> image;
    for(uint32_t y = 0; y < w_height; ++y) {
        for(uint32_t x = 0; x < w_width; ++x) {
            float altitude = world.altitude_at(x, y);

            c_water = c_river;
            c_land.r = clamp(altitude * 255.0f, 0.0f, 255.0f);
            c_land.g = c_land.b = c_land.r;

            World::Tile t = world.at(x, y);

            if(altitude < World::k_lake_altitude) {
                uint8_t green = 255.0f * fabsf((altitude - World::k_altitude_min) * 1.5f);
                c_water.b = c_water.b > green ? c_water.b - green : 0;
                c_water.g = green;
            }

            const World::Area& area = world.get_area(x, y);
            if(t == World::Soil) {
                c_land.r = clamp(area.biome / (float)World::Biome::Count * 255.0f, 0.0f, 255.0f);
                c_land.g = c_land.b = c_land.r;

                switch(area.biome) {
                    case World::Swamp: {
                        c_land = c_swamp;
                        break;
                    }
                    case World::Count: {
                        c_land.r = 1.0f;
                        c_land.g = c_land.b = 0.0f;
                        break;
                    }
                    case World::Plain: {
                        c_land = c_plain;
                        break;
                    }
                    case World::Shrubland: {
                        c_land = c_shrubland;
                        break;
                    }
                    case World::Rocks: {
                        c_land = c_rock;
                        break;
                    }
                    case World::Forest: {
                        c_land = c_forest;
                        break;
                    }
                }
            }

            float h = 0.5f + 0.5f * altitude;
            //c_land.r = clamp(altitude * 255.0f, 0.0f, 255.0f);
            const Vec3& normal = world.normal_at(x, y);
            c_land.r *= normal.z * h;
            c_land.g *= normal.z * h;
            c_land.b *= normal.z * h;

            switch(t) {
                case World::Soil: image.push_back(c_land.r); image.push_back(c_land.g); image.push_back(c_land.b); image.push_back(c_land.a); break;
                case World::Rock: image.push_back(c_rock.r); image.push_back(c_rock.g); image.push_back(c_rock.b); image.push_back(c_rock.a); break;
                case World::Water: image.push_back(c_water.r); image.push_back(c_water.g); image.push_back(c_water.b); image.push_back(c_water.a); break;
                default: break;
            }
        }
    }

    std::vector<uint8_t> png;
    unsigned error = lodepng::encode(png, image, w_width, w_height);
    if(!error) lodepng::save_file(png, "world.png");

    // Exports the terrain's height maps in raw format.
    const std::vector<float> altitude_buffer = world.get_altitude_buffer();
    std::vector<uint8_t> raw_heightmap;
    raw_heightmap.resize(altitude_buffer.size());
    std::vector<uint8_t> raw_heightmap_water;
    raw_heightmap_water.resize(altitude_buffer.size());
    uint32_t i = 0;
    for(uint32_t y = 0; y < w_height; ++y) {
        for(uint32_t x = 0; x < w_width; ++x) {

            float altitude = lerp(min_height, max_height, world.altitude_at(y, x));
            altitude = (altitude - min_height) / (max_height - min_height);
            255.0f * clamp(altitude, 0.0f, 1.0f);
            raw_heightmap[i] = raw_encode_altitude(altitude, min_height, max_height);

            altitude = lerp(min_height, max_height, world.water_altitude_at(y, x));
            altitude = (altitude - min_height) / (max_height - min_height);
            raw_heightmap_water[i] = 255.0f * clamp(altitude, 0.0f, 1.0f);

            ++i;
        }
    }
    uint32_t amplitude = max_height - min_height;
    std::string heightmap_name = "heightmap_" + std::to_string(amplitude) + ".raw";
    raw_write(heightmap_name.c_str(), raw_heightmap);

    heightmap_name = "heightmap_water_" + std::to_string(amplitude) + ".raw";
    raw_write(heightmap_name.c_str(), raw_heightmap_water);

    // Exports the tiles/height/biomes map.
    std::vector<uint8_t> biomes_map;
    Color biome;
    for(uint32_t y = 0; y < w_height; ++y) {
        for(uint32_t x = 0; x < w_width; ++x) {
            auto area = world.get_area(x, y);
            auto tile = world.at(x, y);
            biome.r = tile;
            biome.g = raw_heightmap[y + x * w_width];
            biome.b = area.biome;
            biome.a = 255;

            biomes_map.push_back(biome.r);
            biomes_map.push_back(biome.g);
            biomes_map.push_back(biome.b);
            biomes_map.push_back(biome.a);
        }
    }
    raw_write("biomes.bytes", biomes_map);

    return 0;
}