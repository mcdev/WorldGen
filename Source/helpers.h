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

Vec3 operator - (const Vec3& i, const Vec3& j) {
    return { i.x - j.x, i.y - j.y, i.z - j.z };
}

Vec3 operator / (const Vec3& v, float f) {
    return { v.x / f, v.y / f, v.z / f };
}

Vec3 cross(const Vec3& i, const Vec3& j) {
    return { i.y * j.z - j.y * i.z, i.z * j.x - j.z * i.x, i.x * j.y - j.x * i.y };
}

float length(const Vec3& v) {
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vec3 normalize(const Vec3& v) {
    return v / length(v);
}