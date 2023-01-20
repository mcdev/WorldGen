
#include "lodepng.h"
#include "world.h"
#include "helpers.h"

#include <cmath>
#include <map>
#include <stdio.h>

static uint8_t raw_encode_altitude(float altitude, float min_height, float max_height) {
    float raw_altitude = lerp(min_height, max_height, altitude);
    raw_altitude = (raw_altitude - min_height) / (max_height - min_height);
    return 255.0f * clamp(raw_altitude, 0.0f, 1.0f);
}

static void raw_write(const std::string& filename, const std::vector<uint8_t>& data) {
    FILE* fp;
    fopen_s(&fp, filename.c_str(), "wb");
	if (fp == NULL) {
        return;
    }

	fwrite(data.data(), data.size(), 1, fp);
	fclose(fp);
}

static void tga_write(const std::string& filename, uint32_t width, uint32_t height, uint8_t *dataBGRA, uint8_t dataChannels=4, uint8_t fileChannels=3)
{
	FILE *fp = NULL;
	fopen_s (&fp, filename.c_str(), "wb");
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

static Vec2 average_current(const World& world, int x0, int y0, int x1, int y1) {

    Vec3 h = { 0.0f, 0.0f, 0.0f };
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            h = h + world.water_current_at(y, x);
        }
    }

    Vec2 res = h.xy();
    if (h.z != 0.0f)
        res /= h.z;

    return res;
}

static uint32_t translate_indice(std::map<uint32_t, uint32_t>& indices_map, uint32_t indice,
    const std::vector<float>& vertices, std::vector<float>& new_vertices,
    const std::vector<float>& uv0, std::vector<float>& new_uv0,
    const std::vector<float>& normals, std::vector<float>& new_normals) {

    auto ite = indices_map.find(indice);
    if (ite != indices_map.end()) {
        return ite->second;
    }

    uint32_t new_indice = new_vertices.size() / 3;
    indices_map.insert(std::make_pair(indice, new_indice));

    new_vertices.push_back(vertices[3 * indice + 0]);
    new_vertices.push_back(vertices[3 * indice + 1]);
    new_vertices.push_back(vertices[3 * indice + 2]);

    new_uv0.push_back(uv0[2 * indice + 0]);
    new_uv0.push_back(uv0[2 * indice + 1]);

    new_normals.push_back(normals[3 * indice + 0]);
    new_normals.push_back(normals[3 * indice + 1]);
    new_normals.push_back(normals[3 * indice + 2]);

    return new_indice;
}

static void obj_write(const char* filename, const std::vector<uint8_t>& data, const World& world) {
    FILE* fp;
    fopen_s(&fp, filename, "wt");
    if (fp == NULL) {
        return;
    }

    // Create the geometry from the height map.
    std::vector<float> vertices;
    vertices.resize(3 * (world.width + 1) * (world.height + 1));

    std::vector<float> normals;
    normals.resize(3 * (world.width + 1) * (world.height + 1), 0.0f);

    std::vector<uint32_t> indices;
    indices.resize(2 * 3 * world.width * world.height);

    std::vector<float> uv0;
    uv0.resize(2 * (world.width + 1) * (world.height + 1));

    uint32_t quad = 0;
    for (uint32_t y = 0; y < world.height; ++y) {
        for (uint32_t x = 0; x < world.width; ++x) {

            // Output vertices.
            int xl = (int)x - 1;
            int xr = (int)x + 1;
            int ya = (int)y - 1;
            int yb = (int)y + 1;

            uint32_t i0 = x + y * (world.width + 1);
            uint32_t i1 = i0 + 1;
            uint32_t i3 = x + (y + 1) * (world.width + 1);
            uint32_t i2 = i3 + 1;

            float z0 = average_heights(data, world.width, world.height, xl, ya, x,  y);
            vertices[3 * i0 + 0] = x / (float)world.width;
            vertices[3 * i0 + 1] = z0 / 255.0f;
            vertices[3 * i0 + 2] = y / (float)world.height;

            Vec2 c0 = average_current(world, xl, ya, x, y);
            uv0[2 * i0 + 0] = c0.x;
            uv0[2 * i0 + 1] = c0.y;

            if (x == world.width - 1) {
                float z1 = average_heights(data, world.width, world.height, x, ya, xr, y);
                vertices[3 * i1 + 0] = 1.0f;
                vertices[3 * i1 + 1] = z1 / 255.0f;
                vertices[3 * i1 + 2] = y / (float)world.height;

                Vec2 c1 = average_current(world, x, ya, xr, y);
                uv0[2 * i1 + 0] = c1.x;
                uv0[2 * i1 + 1] = c1.y;
            }

            if (y == world.height - 1) {
                float z3 = average_heights(data, world.width, world.height, xl, y, x, yb);
                vertices[3 * i3 + 0] = x / (float)world.width;
                vertices[3 * i3 + 1] = z3 / 255.0f;
                vertices[3 * i3 + 2] = 1.0f;

                Vec2 c3 = average_current(world, xl, y, x, yb);
                uv0[2 * i3 + 0] = c3.x;
                uv0[2 * i3 + 1] = c3.y;

                if (x == world.width - 1) {
                    float z2 = average_heights(data, world.width, world.height, x, y, xr, yb);
                    vertices[3 * i2 + 0] = 1.0f;
                    vertices[3 * i2 + 1] = z2 / 255.0f;
                    vertices[3 * i2 + 2] = 1.0f;

                    Vec2 c2 = average_current(world, x, y, xr, yb);
                    uv0[2 * i2 + 0] = c2.x;
                    uv0[2 * i2 + 1] = c2.y;
                }
            }

            // Output triangle indices.
            indices[3 * 2 * quad + 0] = i0;
            indices[3 * 2 * quad + 1] = i3;
            indices[3 * 2 * quad + 2] = i2;
            indices[3 * 2 * quad + 3] = i0;
            indices[3 * 2 * quad + 4] = i2;
            indices[3 * 2 * quad + 5] = i1;

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

        Vec3 n = cross(e01, e02);
        n.normalize();

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

    for (auto n = 0; n < normals.size(); n += 3) {
        Vec3 nrm { normals[n + 0], normals[n + 1], normals[n + 2] };
        nrm.normalize();
        normals[n + 0] = nrm.x;
        normals[n + 1] = nrm.y;
        normals[n + 2] = nrm.z;
    }

    // Remove unnecessary quads from the mesh (those under the ground).
    std::vector<float> new_vertices;
    std::vector<float> new_uv0;
    std::vector<float> new_normals;
    std::vector<uint32_t> new_indices;
    std::map<uint32_t, uint32_t> indices_map;

    for (uint32_t y = 0; y < world.height; ++y) {
        for (uint32_t x = 0; x < world.width; ++x) {

            if (world.max_water_altitude_at(y, x, 1) < world.min_altitude_at(y, x, 1)) {
                // Skip underground quads.
                continue;
            }

            uint32_t i0 = x + y * (world.width + 1);
            uint32_t i1 = i0 + 1;
            uint32_t i3 = x + (y + 1) * (world.width + 1);
            uint32_t i2 = i3 + 1;

            i0 = translate_indice(indices_map, i0, vertices, new_vertices, uv0, new_uv0, normals, new_normals);
            i1 = translate_indice(indices_map, i1, vertices, new_vertices, uv0, new_uv0, normals, new_normals);
            i2 = translate_indice(indices_map, i2, vertices, new_vertices, uv0, new_uv0, normals, new_normals);
            i3 = translate_indice(indices_map, i3, vertices, new_vertices, uv0, new_uv0, normals, new_normals);

            // Output triangle indices.
            new_indices.push_back(i0);
            new_indices.push_back(i3);
            new_indices.push_back(i2);
            new_indices.push_back(i0);
            new_indices.push_back(i2);
            new_indices.push_back(i1);
        }
    }

    printf("%.02f%% of the water mesh has been filtered.\n", (float)(indices.size() - new_indices.size()) * 100.0f / indices.size());

    fprintf(fp, "# Obj file generated by wgen.exe.\n");

    // Output the geometry to the file.
    for (auto v = 0, c = 0; v < new_vertices.size(); v += 3, c += 2) {
        fprintf(fp, "v %f %f %f\n", new_vertices[v + 0], new_vertices[v + 1], new_vertices[v + 2]);
    }
    
    for (auto vn = 0; vn < new_normals.size(); vn += 3) {
        fprintf(fp, "vn %f %f %f\n", new_normals[vn + 0], new_normals[vn + 1], new_normals[vn + 2]);
    }

    for (auto uv = 0; uv < new_uv0.size(); uv += 2) {
        fprintf(fp, "vt %f %f\n", new_uv0[uv + 0], new_uv0[uv + 1]);
    }
    
    // OBJ indices start at 1.
    for (auto i = 0; i < new_indices.size(); i += 3) {
        fprintf(fp, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", 
            1 + new_indices[i + 0], 1 + new_indices[i + 0], 1 + new_indices[i + 0], 
            1 + new_indices[i + 1], 1 + new_indices[i + 1], 1 + new_indices[i + 1], 
            1 + new_indices[i + 2], 1 + new_indices[i + 2], 1 + new_indices[i + 2]);
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

            if(altitude < k_lake_altitude) {
                uint8_t green = 255.0f * fabsf((altitude - k_altitude_min) * 1.5f);
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
                        c_land.r = 255;
                        c_land.g = c_land.b = 0;
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

    const std::string k_path = "D:\\101games\\Unity\\Wild01\\Assets\\Terrains\\";

    std::vector<uint8_t> png;
    unsigned error = lodepng::encode(png, image, w_width, w_height);
    if(!error) lodepng::save_file(png, k_path + "world.png");

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
    std::string heightmap_name = k_path + "heightmap_" + std::to_string(amplitude) + ".raw";
    raw_write(heightmap_name.c_str(), raw_heightmap);

    heightmap_name = k_path + "heightmap_water_" + std::to_string(amplitude) + ".raw";
    raw_write(heightmap_name.c_str(), raw_heightmap_water);

    std::string obj_name = k_path + "heightmap_water_" + std::to_string(amplitude) + ".obj";
    obj_write(obj_name.c_str(), raw_heightmap_water, world);

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
    raw_write(k_path + "biomes.bytes", biomes_map);

    printf("Hit ENTER to quit.\n");
    //getc(stdin);

    return 0;
}