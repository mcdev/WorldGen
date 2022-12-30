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

template <class T>
T remap(T v, T min, T max, T new_min, T new_max) {
    return new_min + (clamp(v, min, max) - min) / (max - min) * new_max;
}

float saturate(float v) {
    return clamp(v, 0.0f, 1.0f);
}

struct Color { uint8_t r, g, b, a; };

struct Vec2 { 
    void operator += (const Vec2& other) {
        x += other.x;
        y += other.y;
    }

    void operator /= (const float f) {
        x /= f;
        y /= f;
    }

    float length() {
        return sqrtf(x * x + y * y);
    }

    void normalize() {
        auto l = length();
        if (l == 0.0f) {
            return;
        }
        *this /= l;
    }

    float x, y;
};

Vec2 operator + (const Vec2& i, const Vec2& j) {
    return { i.x + j.x, i.y + j.y };
}

Vec2 operator / (const Vec2& v, float f) {
    return { v.x / f, v.y / f };
}

Vec2 operator * (const Vec2& i, float f) {
    return { i.x * f, i.y * f };
}

static const Vec2 zero2 = { 0.0f, 0.0f };

struct Vec3 {

    Vec3(const Vec2& xy, float z = 0.0f) {
        x = xy.x;
        y = xy.y;
        this->z = z;
    }

    Vec3(float x, float y, float z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Vec3() = default;

    Vec2 xy() const {

        return { x, y };
    }

    void operator += (const Vec3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
    }

    void operator *= (const float f) {
        x *= f;
        y *= f;
        z *= f;
    }

    void operator /= (const float f) {
        x /= f;
        y /= f;
        z /= f;
    }

    float length() {
        return sqrtf(x * x + y * y + z * z);
    }
    
    void normalize() {
        auto l = length();
        if (l == 0.0f) {
            return;
        }
        *this /= l;
    }

    float dot(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    void ortho_basis(Vec3& b1, Vec3& b2)
    {
        float sign = copysignf(1.0f, z);
        const float a = -1.0f / (sign + z);
        const float b = x * y * a;
        b1 = Vec3(1.0f + sign * x * x * a, sign * b, -sign * x);
        b2 = Vec3(b, sign + y * y * a, -y);
    }

    float x, y, z;
};

static const Vec3 x_axis = { 1.0f, 0.0f, 0.0f };
static const Vec3 y_axis = { 0.0f, 1.0f, 0.0f };
static const Vec3 z_axis = { 0.0f, 0.0f, 1.0f };
static const Vec3 zero3  = { 0.0f, 0.0f, 0.0f };

Vec3 operator - (const Vec3& i, const Vec3& j) {
    return { i.x - j.x, i.y - j.y, i.z - j.z };
}

Vec3 operator + (const Vec3& i, const Vec3& j) {
    return { i.x + j.x, i.y + j.y, i.z + j.z };
}

Vec3 operator / (const Vec3& v, float f) {
    return { v.x / f, v.y / f, v.z / f };
}

Vec3 operator * (const Vec3& v, float f) {
    return { v.x * f, v.y * f, v.z * f };
}

float dot(const Vec3& i, const Vec3& j) {
    return i.x * j.x + i.y * j.y + i.z * j.z;
}

Vec3 cross(const Vec3& i, const Vec3& j) {
    return { i.y * j.z - j.y * i.z, i.z * j.x - j.z * i.x, i.x * j.y - j.x * i.y };
}

