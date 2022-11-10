
#include "lodepng.h"
#include "world.h"
#include "helpers.h"

#include <math.h>


int main(int argc, char** argv) {

    const uint32_t w_width = 512;
    const uint32_t w_height = 512;

    World world(w_width, w_height);
    world.init();
    world.compute_rivers_and_lakes();

    // Translates the world's map into an image.
    Color c_land = { 0, 255, 0, 255 };
    Color c_river = { 0, 0, 255, 255 };
    Color c_water;

    std::vector<uint8_t> image;
    for(uint32_t y = 0; y < w_height; ++y) {
        for(uint32_t x = 0; x < w_width; ++x) {
            float altitude = world.altitude_at(x, y);

            c_water = c_river;
            c_land.r = clamp(altitude * 255.0f, 0.0f, 255.0f);
            c_land.g = c_land.b = c_land.r;

            Tile t = world.at(x, y);
            if(world.water_at(x, y) != World::k_no_water) {
                t = Water;
            }

            if(altitude < World::k_lake_altitude) {
                uint8_t green = 255.0f * fabsf((altitude - World::k_altitude_min) * 3.0f);
                c_water.b = c_water.b > green ? c_water.b - green : 0;
                c_water.g = green;
            }

            switch(t) {
                case Land: image.push_back(c_land.r); image.push_back(c_land.g); image.push_back(c_land.b); image.push_back(c_land.a); break;
                case Water: image.push_back(c_water.r); image.push_back(c_water.g); image.push_back(c_water.b); image.push_back(c_water.a); break;
                default: break;
            }
        }
    }

    std::vector<uint8_t> png;
    unsigned error = lodepng::encode(png, image, w_width, w_height);
    if(!error) lodepng::save_file(png, "world.png");

    return 0;
}