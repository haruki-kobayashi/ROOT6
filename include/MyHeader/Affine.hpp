#pragma once
#include <cmath>

// Affine変換のクラス
class Affine {
public:
    double a, b, c, d, p, q;

    Affine() : a(1.0), b(0.0), c(0.0), d(1.0), p(0.0), q(0.0) {}
    Affine(double a_, double b_, double c_, double d_, double p_, double q_)
        : a(a_), b(b_), c(c_), d(d_), p(p_), q(q_) {}

    // 単位変換行列を返す
        static Affine identity() { return Affine(); }

    // 行列式を計算
    double determinant() const { return a * d - b * c; }

    // 逆行列が存在するかどうかを判定
    bool isInvertible(double eps = 1e-12) const { return std::fabs(determinant()) > eps; }

    // 逆行列を計算
    Affine inverse(double eps = 1e-12) const {
        double det = determinant();
        if (std::fabs(det) <= eps) {
            return Affine::identity();
        }
        Affine inv;
        inv.a =  d / det;
        inv.b = -b / det;
        inv.c = -c / det;
        inv.d =  a / det;
        inv.p = -(inv.a * p + inv.b * q);
        inv.q = -(inv.c * p + inv.d * q);
        return inv;
    }

    // 合成行列を計算
    // this ∘ pre （先に pre を適用し、次に this を適用）
    Affine compose(const Affine& pre) const {
        Affine out;
        out.a = a * pre.a + b * pre.c;
        out.b = a * pre.b + b * pre.d;
        out.c = c * pre.a + d * pre.c;
        out.d = c * pre.b + d * pre.d;
        out.p = p + a * pre.p + b * pre.q;
        out.q = q + c * pre.p + d * pre.q;
        return out;
    }

    // 座標変換
    void transform(double xin, double yin, double* xout, double* yout) const {
        *xout = a * xin + b * yin + p;
        *yout = c * xin + d * yin + q;
    }
    // 逆座標変換
    void transformInverse(double xin, double yin, double* xout, double* yout, double eps = 1e-12) const {
        Affine inv = inverse(eps);
        *xout = inv.a * xin + inv.b * yin + inv.p;
        *yout = inv.c * xin + inv.d * yin + inv.q;
    }

    // float版の座標変換
    void transform(float xin, float yin, float* xout, float* yout) const {
        *xout = static_cast<float>(a * xin + b * yin + p);
        *yout = static_cast<float>(c * xin + d * yin + q);
    }
    // float版の逆座標変換
    void transformInverse(float xin, float yin, float* xout, float* yout, double eps = 1e-12) const {
        Affine inv = inverse(eps);
        *xout = static_cast<float>(inv.a * xin + inv.b * yin + inv.p);
        *yout = static_cast<float>(inv.c * xin + inv.d * yin + inv.q);
    }

    Affine operator*(const Affine& pre) const { return compose(pre); }
    Affine& operator*=(const Affine& pre) { return (*this = compose(pre)); }

    Affine operator+(const Affine& rhs) const {
        Affine out(*this);
        out.a += rhs.a; out.b += rhs.b; out.c += rhs.c; out.d += rhs.d; out.p += rhs.p; out.q += rhs.q;
        return out;
    }
    Affine& operator+=(const Affine& rhs) {
        a += rhs.a; b += rhs.b; c += rhs.c; d += rhs.d; p += rhs.p; q += rhs.q;
        return *this;
    }

    Affine operator*(double k) const {
        Affine out(*this);
        out.a *= k; out.b *= k; out.c *= k; out.d *= k; out.p *= k; out.q *= k;
        return out;
    }
    Affine& operator*=(double k) {
        a *= k; b *= k; c *= k; d *= k; p *= k; q *= k;
        return *this;
    }
    Affine& operator/=(double k) { return (*this *= (1.0 / k)); }
};

inline Affine compose(const Affine& post, const Affine& pre) { return post.compose(pre); }
inline Affine invert(const Affine& af, double eps = 1e-12) { return af.inverse(eps); }
