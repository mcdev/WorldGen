#pragma once

#include "helpers.h"

namespace RandGen {

    int32_t permute(uint32_t i, uint32_t l, uint32_t p) {
        uint32_t w = l - 1;
        w |= w >> 1;
        w |= w >> 2;
        w |= w >> 4;
        w |= w >> 8;
        w |= w >> 16;
        do {
            i ^= p;
            i *= 0xe170893d;
            i ^= p >> 16;
            i ^= (i & w) >> 4;
            i ^= p >> 8;
            i *= 0x0929eb3f;
            i ^= p >> 23;
            i ^= (i & w) >> 1;
            i *= 1 | p >> 27;
            i *= 0x6935fa69;
            i ^= (i & w) >> 11;
            i *= 0x74dcb303;
            i ^= (i & w) >> 2;
            i *= 0x9e501cc3;
            i ^= (i & w) >> 2;
            i *= 0xc860a3df;
            i &= w;
            i ^= i >> 5;
        } while(i >= l);
        return int32_t((i + p) % l);
    }

    float randfloat(uint32_t i, uint32_t p) {
        i ^= p; 
        i ^= i >> 17;
        i ^= i >> 10;
        i *= 0xb36534e5;
        i ^= i >> 12;
        i ^= i >> 21;
        i *= 0x93fc4795;
        i ^= 0xdf6e307f;
        i ^= i >> 17;
        i *= 1 | p >> 18;
        return i * (1.0f / 4294967808.0f);
    }

    Vec2 sample_2d(int s, int m, int n, int p) {
        s = permute(s, m * n, p * 0x51633e2d); 
        int sx = permute(s % m, m, p * 0x68bc21eb); 
        int sy = permute(s / m, n, p * 0x02e5be93); 
        float jx = randfloat(s, p * 0x967a889b); 
        float jy = randfloat(s, p * 0x368cc8b7); 
        Vec2 r = { (sx + (sy + jx) / n) / m, (s + jy) / (m * n) };
        return r; 
    }
};
