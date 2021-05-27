#ifndef COLORF_H
#define COLORF_H

#include <vector>
#include <cmath>

// Color class
class Colorf {
public:
    float r;
    float g;
    float b;
    Colorf() : r(0.0f), g(0.0f), b(0.0f) {};
    Colorf(const Colorf& C) : r(C.r), g(C.g), b(C.b) {};
    Colorf(float a, float b, float c) : r(a), g(b), b(c) {};
    Colorf(float a) : r(a), g(a), b(a) {};
    ~Colorf() {};
    Colorf operator+(const Colorf& other) const
    {
        return Colorf(this->r + other.r, this->g + other.g, this->b + other.b);
    };
    Colorf operator*(const float c) const
    {
        return Colorf(this->r * c, this->g * c, this->b * c);
    };
};
// color data to ramp between
const std::vector<Colorf> PortalW{
    Colorf(1.        , 0.36470588, 0.),
    Colorf(1.        , 0.60392157, 0.),
    Colorf(1.        , 1.        , 1.),
    Colorf(0.        , 0.63529412, 1.),
    Colorf(0.        , 0.39607843, 1.)
};
// clamp values
inline float clamp(const float x, const float x_min, const float x_max) {
    return std::min(std::max(x, x_min), x_max);
}
//  performs smooth Hermite interpolation between 0 and 1 when e0 < x < e1
float smoothstep(const float e0, const float e1, const float x) {
    float t = clamp((x - e0) / (e1 - e0), 0.0f, 1.0f);
    return t * t * (3.0f - 2.0f * t);
}
// return fractional part of number
float fract(const float t) {
    return t - std::floor(t);
}
// performs a linear interpolation between X and Y using a weight w between them. 
Colorf mix(const Colorf& X, const Colorf& Y, const float w) {
    return X * (1.0f - w) + Y * w;
}
// ramp between colors
Colorf ramp(const std::vector<Colorf>& Cols, const float t) {
    float t_temp = t * float(Cols.size() - 1);
    return mix(Cols[int(t_temp)], Cols[int(t_temp) + 1], smoothstep(0.0f, 1.0f, fract(t_temp)));
}
// clamp color components between low and high
Colorf clamp(const Colorf& color, const float low = 0.0f, const float up = 1.0f)
{
    return Colorf(std::max(std::min(color.r, up), low),
        std::max(std::min(color.g, up), low),
        std::max(std::min(color.b, up), low));
};

#endif
