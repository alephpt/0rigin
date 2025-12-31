// include/physics/PhysicsTypes.h
#pragma once
#include <array>

namespace physics {

struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
};

struct Vec4 {
    float t, x, y, z;
    Vec4() : t(0), x(0), y(0), z(0) {}
    Vec4(float t, float x, float y, float z) : t(t), x(x), y(y), z(z) {}
};

struct FieldTensor {
    // F^μν antisymmetric tensor (6 independent components)
    float Ex, Ey, Ez;  // Electric field
    float Bx, By, Bz;  // Magnetic field
};

} // namespace physics
