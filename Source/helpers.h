#pragma once

#include <stdlib.h>

template <class T>
T clamp(T v, T min, T max) {
    if(v > max) return max;
    if(v < min) return min;
    return v;
}

template <class T>
T lerp(T min, T max, float t) {
    if (t <= 0.0f) return min;
    else if (t >= 1.0f) return max;
    else return min * (1.0f - t) + max * t;
}

struct Color { uint8_t r, g, b, a; };

struct Vec2 { float x, y; };
struct Vec3 { float x, y, z; };
