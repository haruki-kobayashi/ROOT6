#pragma once

#include <algorithm>

#include <MyHeader/Vec3.hpp>

namespace {
    constexpr double EPS = 1e-10;
}

namespace CalcTrackPair
{
    // Radial方向の角度差を計算
    // tanθ < EPS のときは √(dax^2+day^2) を返す
    template <class Track1, class Track2>
    double RadialAngleDifference(const Track1& t1, const Track2& t2)
    {
        double denom = std::sqrt(t1.ax * t1.ax + t1.ay * t1.ay);
        double dax = t2.ax - t1.ax;
        double day = t2.ay - t1.ay;
        return (denom < EPS) ? std::hypot(dax, day) : (dax * t1.ax + day * t1.ay) / denom;
    }

    // Lateral方向の角度差を計算
    // tanθ < EPS のときは 0.0 を返す
    template <class Track1, class Track2>
    double LateralAngleDifference(const Track1& t1, const Track2& t2)
    {
        double denom = std::sqrt(t1.ax * t1.ax + t1.ay * t1.ay);
        double dax = t2.ax - t1.ax;
        double day = t2.ay - t1.ay;
        return (denom < EPS) ? 0.0 : (dax * t1.ay - day * t1.ax) / denom;
    }

    // 開き角を計算（ラジアン）
    template <class Track1, class Track2>
    double OpeningAngle(const Track1& t1, const Track2& t2)
    {
        // 各trackの方向ベクトル
        Vec3 dirV0 = {t1.ax, t1.ay, 1};
        Vec3 dirV1 = {t2.ax, t2.ay, 1};

        double cos_theta = dirV0.dot(dirV1) / (dirV0.length() * dirV1.length());

        // 浮動小数点誤差で値が範囲外になるのを防ぐ
        cos_theta = std::clamp(cos_theta, -1.0, 1.0);

        return std::acos(cos_theta);
    }

    // Radial方向の距離を計算
    template <class Track1, class Track2>
    double RadialDistance(const Track1& t1, const Track2& t2)
    {
        // 各trackの位置ベクトルと方向ベクトル
        Vec3 posV0 = {t1.x, t1.y, t1.z};
        Vec3 posV1 = {t2.x, t2.y, t2.z};
        Vec3 dirV0 = {t1.ax, t1.ay, 1};
        Vec3 dirV1 = {t2.ax, t2.ay, 1};

        // 中点と差分
        Vec3 mid = (posV0 + posV1) * 0.5;
        Vec3 diff = posV1 - posV0;
        if (diff.length() < EPS) return 0.0;

        // 中間面に位置ベクトルを外挿
        double coef0, coef1;
        double denom0 = dirV0.dot(diff);
        double denom1 = dirV1.dot(diff);
        if (std::abs(denom0) < EPS) { // 外挿不能->元の位置
            coef0 = 0.0;
        } else {
            coef0 = (mid - posV0).dot(diff) / denom0;
        }
        if (std::abs(denom1) < EPS) {
            coef1 = 0.0;
        } else {
            coef1 = (mid - posV1).dot(diff) / denom1;
        }
        Vec3 ext0 = posV0 + (dirV0 * coef0);
        Vec3 ext1 = posV1 + (dirV1 * coef1);

        // lateral方向の単位ベクトル
        Vec3 z_axis = {0, 0, 1};
        Vec3 lateral = z_axis.cross(diff);
        Vec3 unit_r = diff.cross(lateral).normalize();

        // radial方向の変位
        return (ext1 - ext0).dot(unit_r);
    }

    // Lateral方向の距離を計算
    template <class Track1, class Track2>
    double LateralDistance(const Track1& t1, const Track2& t2)
    {
        // 各trackの位置ベクトルと方向ベクトル
        Vec3 posV0 = {t1.x, t1.y, t1.z};
        Vec3 posV1 = {t2.x, t2.y, t2.z};
        Vec3 dirV0 = {t1.ax, t1.ay, 1};
        Vec3 dirV1 = {t2.ax, t2.ay, 1};

        // 中点と差分
        Vec3 mid = (posV0 + posV1) * 0.5;
        Vec3 diff = posV1 - posV0;
        if (diff.length() < EPS) return 0.0;

        // 中間面に位置ベクトルを外挿
        double coef0, coef1;
        double denom0 = dirV0.dot(diff);
        double denom1 = dirV1.dot(diff);
        if (std::abs(denom0) < EPS) { // 外挿不能->元の位置
            coef0 = 0.0;
        } else {
            coef0 = (mid - posV0).dot(diff) / denom0;
        }
        if (std::abs(denom1) < EPS) {
            coef1 = 0.0;
        } else {
            coef1 = (mid - posV1).dot(diff) / denom1;
        }
        Vec3 ext0 = posV0 + (dirV0 * coef0);
        Vec3 ext1 = posV1 + (dirV1 * coef1);

        // lateral方向の単位ベクトル
        Vec3 z_axis = {0, 0, 1};
        Vec3 unit_l = z_axis.cross(diff).normalize();

        // lateral方向の変位
        return (ext1 - ext0).dot(unit_l);
    }

    // 最近接点の距離を計算
    template <class Track1, class Track2>
    double MinimumDistance(const Track1& t1, const Track2& t2)
    {
        // 各trackの位置ベクトルと方向ベクトル
        Vec3 posV0 = {t1.x, t1.y, t1.z};
        Vec3 posV1 = {t2.x, t2.y, t2.z};
        Vec3 dirV0 = {t1.ax, t1.ay, 1};
        Vec3 dirV1 = {t2.ax, t2.ay, 1};

        // 差分
        Vec3 diff = posV1 - posV0;
        if (diff.length() < EPS) return 0.0;

        double coef0, coef1, denom;

        if (OpeningAngle(t1, t2) < EPS) { // 平行な場合
            // trackの方向ベクトルに垂直なdiffベクトルの成分を取り出す
            Vec3 perp_component = diff - dirV0 * diff.dot(dirV0);
            return perp_component.length();
        } else {
            denom = dirV0.dot(dirV0) * dirV1.dot(dirV1) - dirV0.dot(dirV1) * dirV0.dot(dirV1);
            if (std::abs(denom) < EPS) {
                // 平行な場合の処理と同じ
                Vec3 perp_component = diff - dirV0 * diff.dot(dirV0);
                return perp_component.length();
            } else {
                coef0 = (diff.dot(dirV0) * dirV1.dot(dirV1) - diff.dot(dirV1) * dirV0.dot(dirV1)) / denom;
                coef1 = (diff.dot(dirV0) * dirV0.dot(dirV1) - diff.dot(dirV1) * dirV0.dot(dirV0)) / denom;
            }
        }

        // 外挿した位置ベクトルが両端より外側に出ないように補正
        // 最近接点が両端より外側にある場合は両端に固定する
        double z_low = std::min(posV0.z, posV1.z);
        double z_high = std::max(posV0.z, posV1.z);
        coef0 = std::clamp(coef0, z_low - posV0.z, z_high - posV0.z);
        coef1 = std::clamp(coef1, z_low - posV1.z, z_high - posV1.z);

        // 外挿した位置ベクトル (最近接点)
        Vec3 intercept0 = posV0 + dirV0 * coef0;
        Vec3 intercept1 = posV1 + dirV1 * coef1;

        return (intercept1 - intercept0).length();
    }
}
