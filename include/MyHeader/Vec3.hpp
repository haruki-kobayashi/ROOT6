#pragma once

#include <cmath>

// 3次元ベクトルのクラス
class Vec3 {
public:
    double x, y, z;
    Vec3(double x = 0., double y = 0., double z = 0.) : x(x), y(y), z(z) {}

    // ベクトル演算
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; } // スカラー倍

    // 内積（dot product）
    double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    // 外積（cross product）
    Vec3 cross(const Vec3& o) const {
        return { y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x };
    }

    // ベクトルの大きさ
    double length() const { return std::sqrt(x * x + y * y + z * z); }

    // 単位ベクトルに正規化
    Vec3 normalized() const {
        double n = length();
        // ゼロ除算を避ける
        return (n < 1e-10) ? Vec3{0, 0, 0} : Vec3{x / n, y / n, z / n};
    }
    Vec3& normalize() {
        double n = length();
        if (n > 1e-10) {
            x /= n; y /= n; z /= n;
        } else {
            x = 0; y = 0; z = 0;
        }
        return *this;
    }
};
