
#include "lodepng.h"
#include "world.h"
#include "helpers.h"

#include <cmath>
#include <stdio.h>

static uint8_t raw_encode_altitude(float altitude, float min_height, float max_height) {
    float raw_altitude = lerp(min_height, max_height, altitude);
    raw_altitude = (raw_altitude - min_height) / (max_height - min_height);
    return 255.0f * clamp(raw_altitude, 0.0f, 1.0f);
}

static void raw_write(const char *filename, const std::vector<uint8_t>& data) {
    FILE* fp;
    fopen_s(&fp, filename, "wb");
	if (fp == NULL) {
        return;
    }

	fwrite(data.data(), data.size(), 1, fp);
	fclose(fp);
}

static void tga_write(const char *filename, uint32_t width, uint32_t height, uint8_t *dataBGRA, uint8_t dataChannels=4, uint8_t fileChannels=3)
{
	FILE *fp = NULL;
	fopen_s (&fp, filename, "wb");
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

static float average_heights(const std::vector<uint8_t>& data, uint32_t width, uint32_t height, int x0, int y0, int x1, int y1) {

    float count = 0.0f;
    float h = 0.0f;
    for(int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            if ((x < 0) || (y < 0) || (x >= width) || (y >= height)) {
                continue;
            }
            ++count;
            h += data[y + x * width];
        }
    }

    return h / count;
}

static void obj_write(const char* filename, const std::vector<uint8_t>& data, uint32_t width, uint32_t height) {
    FILE* fp;
    fopen_s(&fp, filename, "wt");
    if (fp == NULL) {
        return;
    }

    // Create the geometry from the height map.
    std::vector<float> vertices;
    vertices.resize(3 * (width + 1) * (height + 1));

    std::vector<float> normals;
    normals.resize(3 * (width + 1) * (height + 1), 0.0f);

    std::vector<uint32_t> indices;
    indices.resize(2 * 3 * width * height);

    std::vector<float> uv0;
    uv0.resize(2 * (width + 1) * (height + 1));

    uint32_t quad = 0;
    for (uint32_t y = 0; y < height; ++y) {
        for (uint32_t x = 0; x < width; ++x) {

            // Output vertices.
            int xl = (int)x - 1;
            int xr = (int)x + 1;
            int ya = (int)y - 1;
            int yb = (int)y + 1;

            uint32_t i0 = x + y * (width + 1);
            uint32_t i1 = i0 + 1;
            uint32_t i3 = x + (y + 1) * (width + 1);
            uint32_t i2 = i3 + 1;

            float z0 = average_heights(data, width, height, xl, ya, x,  y);
            vertices[3 * i0 + 0] = x / (float)width;
            vertices[3 * i0 + 1] = z0 / 255.0f;
            vertices[3 * i0 + 2] = y / (float)height;

            if (x == width - 1) {
                float z1 = average_heights(data, width, height, x, ya, xr, y);
                vertices[3 * i1 + 0] = 1.0f;
                vertices[3 * i1 + 1] = z1 / 255.0f;
                vertices[3 * i1 + 2] = y / (float)height;
            }

            if (y == height - 1) {
                float z3 = average_heights(data, width, height, xl, y, x, yb);
                vertices[3 * i3 + 0] = x / (float)width;
                vertices[3 * i3 + 1] = z3 / 255.0f;
                vertices[3 * i3 + 2] = 1.0f;

                if (x == width - 1) {
                    float z2 = average_heights(data, width, height, x, y, xr, yb);
                    vertices[3 * i2 + 0] = 1.0f;
                    vertices[3 * i2 + 1] = z2 / 255.0f;
                    vertices[3 * i2 + 2] = 1.0f;
                }
            }

            // Output triangle indices.
            indices[3 * 2 * quad + 0] = i0;
            indices[3 * 2 * quad + 1] = i3;
            indices[3 * 2 * quad + 2] = i2;
            indices[3 * 2 * quad + 3] = i0;
            indices[3 * 2 * quad + 4] = i2;
            indices[3 * 2 * quad + 5] = i1;

            // Output texture coordinates.
            uv0[2 * i0 + 0] = x / (float)width;
            uv0[2 * i0 + 1] = y / (float)height;

            uv0[2 * i1 + 0] = (x + 1) / (float)width;
            uv0[2 * i1 + 1] = y / (float)height;

            uv0[2 * i2 + 0] = (x + 1) / (float)width;
            uv0[2 * i2 + 1] = (y + 1) / (float)height;

            uv0[2 * i3 + 0] = x / (float)width;
            uv0[2 * i3 + 1] = (y + 1) / (float)height;

            ++quad;
        }
    }

    // Compute vertex normals.
    for (uint32_t f = 0; f < indices.size(); f += 3) {
        auto i0 = indices[f + 0];
        auto i1 = indices[f + 1];
        auto i2 = indices[f + 2];

        Vec3 v0;
        v0.x = vertices[3 * i0 + 0];
        v0.y = vertices[3 * i0 + 1];
        v0.z = vertices[3 * i0 + 2];

        Vec3 v1;
        v1.x = vertices[3 * i1 + 0];
        v1.y = vertices[3 * i1 + 1];
        v1.z = vertices[3 * i1 + 2];

        Vec3 v2;
        v2.x = vertices[3 * i2 + 0];
        v2.y = vertices[3 * i2 + 1];
        v2.z = vertices[3 * i2 + 2];

        Vec3 e01 = v1 - v0;
        Vec3 e02 = v2 - v0;

        Vec3 n = normalize(cross(e01, e02));

        normals[3 * i0 + 0] += n.x;
        normals[3 * i0 + 1] += n.y;
        normals[3 * i0 + 2] += n.z;

        normals[3 * i1 + 0] += n.x;
        normals[3 * i1 + 1] += n.y;
        normals[3 * i1 + 2] += n.z;

        normals[3 * i2 + 0] += n.x;
        normals[3 * i2 + 1] += n.y;
        normals[3 * i2 + 2] += n.z;
    }

    for (uint32_t n = 0; n < normals.size(); n += 3) {
        Vec3 nrm { normals[n + 0], normals[n + 1], normals[n + 2] };
        nrm = normalize(nrm);
        normals[n + 0] = nrm.x;
        normals[n + 1] = nrm.y;
        normals[n + 2] = nrm.z;
    }

    fprintf(fp, "# Obj file generateds by wgen.exe.\n");

    // Output the geometry to the file.
    for (auto v = 0; v < vertices.size(); v += 3) {
        fprintf(fp, "v %f %f %f\n", vertices[v + 0], vertices[v + 1], vertices[v + 2]);
    }
    
    for (auto vn = 0; vn < normals.size(); vn += 3) {
        fprintf(fp, "vn %f %f %f\n", normals[vn + 0], normals[vn + 1], normals[vn + 2]);
    }

    for (auto uv = 0; uv < uv0.size(); uv += 2) {
        fprintf(fp, "vt %f %f\n", uv0[uv + 0], uv0[uv + 1]);
    }
    
    // OBJ indices start at 1.
    for (auto i = 0; i < indices.size(); i += 3) {
        fprintf(fp, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", 1 + indices[i + 0], 1 + indices[i + 0], 1 + indices[i + 0], 1 + indices[i + 1], 1 + indices[i + 1], 1 + indices[i + 1], 1 + indices[i + 2], 1 + indices[i + 2], 1 + indices[i + 2]);
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

    printf("Generating terrain in [%f, %f]...\n", min_height, max_height);

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
            raw_heightmap[i] = raw_encode_altitude(world.altitude_at(x, y), min_height, max_height);
            raw_heightmap_water[i] = raw_encode_altitude(world.water_altitude_at(x, y), min_height, max_height);
            ++i;
        }
    }
    uint32_t amplitude = max_height - min_height;
    std::string heightmap_name = "heightmap_" + std::to_string(amplitude) + ".raw";
    raw_write(heightmap_name.c_str(), raw_heightmap);

    heightmap_name = "heightmap_water_" + std::to_string(amplitude) + ".raw";
    raw_write(heightmap_name.c_str(), raw_heightmap_water);

    std::string obj_name = "heightmap_water_" + std::to_string(amplitude) + ".obj";
    obj_write(obj_name.c_str(), raw_heightmap_water, w_width, w_height);

    // Exports the tiles/height/biomes map.
    std::vector<uint8_t> biomes_map;
    Color biome;
    for(uint32_t y = 0; y < w_height; ++y) {
        for(uint32_t x = 0; x < w_width; ++x) {
            auto area = world.get_area(x, y);
            auto tile = world.at(x, y);
            biome.r = tile;
            biome.g = raw_heightmap[y + x * w_width];
            biome.b = tile == World::Water ? World::Underwater : area.biome;
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