#pragma once

#include <stdlib.h>

template <class T>
T clamp(T v, T min, T max) {
    if(v > max) return max;
    if(v < min) return min;
    return v;
}

struct Color { uint8_t r, g, b, a; };
