#pragma once

#include <algorithm>

#include <MyHeader/Vec3.hpp>

namespace {
    constexpr double PI = 3.14159265358979323846;
    constexpr double EPS = 1e-10;
}

namespace CalcTrackPair
{
    template <class Track1, class Track2>
    double OpeningAngle(const Track1& t1, const Track2& t2)
    {
        double oa;

        // 各trackの方向ベクトル
        Vec3 dirV0 = {t1.ax, t1.ay, 1};
        Vec3 dirV1 = {t2.ax, t2.ay, 1};

        double cos_theta = dirV0.dot(dirV1) / (dirV0.length() * dirV1.length());
        if (cos_theta > 1.0) oa = 0.0;
        else if (cos_theta < -1.0) oa = PI;
        else oa = std::acos(cos_theta);

        return oa;
    }

    template <class Track1, class Track2>
    double RadialDistance(const Track1& t1, const Track2& t2)
    {
        double dr;

        // 各trackの位置ベクトルと方向ベクトル
        Vec3 posV0 = {t1.x, t1.y, t1.z};
        Vec3 posV1 = {t2.x, t2.y, t2.z};
        Vec3 dirV0 = {t1.ax, t1.ay, 1};
        Vec3 dirV1 = {t2.ax, t2.ay, 1};

        // 中点と差分
        Vec3 mid = (posV0 + posV1) * 0.5;
        Vec3 diff = posV1 - posV0;

        // 中間面に外挿した位置ベクトル
        double coef0 = (mid - posV0).dot(diff) / dirV0.dot(diff);
        double coef1 = (mid - posV1).dot(diff) / dirV1.dot(diff);
        Vec3 ext0 = posV0 + (dirV0 * coef0);
        Vec3 ext1 = posV1 + (dirV1 * coef1);

        // lateral方向の単位ベクトル
        Vec3 z_axis = {0, 0, 1};
        Vec3 lateral = z_axis.cross(diff);
        Vec3 unit_r = dirV0.cross(lateral).normalize();

        // radial方向の変位
        dr = (ext1 - ext0).dot(unit_r);

        return dr;
    }

    template <class Track1, class Track2>
    double LateralDistance(const Track1& t1, const Track2& t2)
    {
        double dl;

        // 各trackの位置ベクトルと方向ベクトル
        Vec3 posV0 = {t1.x, t1.y, t1.z};
        Vec3 posV1 = {t2.x, t2.y, t2.z};
        Vec3 dirV0 = {t1.ax, t1.ay, 1};
        Vec3 dirV1 = {t2.ax, t2.ay, 1};

        // 中点と差分
        Vec3 mid = (posV0 + posV1) * 0.5;
        Vec3 diff = posV1 - posV0;

        // 中間面に外挿した位置ベクトル
        double coef0 = -1 * (mid - posV0).dot(diff) / dirV0.dot(diff);
        double coef1 = -1 * (mid - posV1).dot(diff) / dirV1.dot(diff);
        Vec3 ext0 = posV0 + (dirV0 * coef0);
        Vec3 ext1 = posV1 + (dirV1 * coef1);

        // lateral方向の単位ベクトル
        Vec3 z_axis = {0, 0, 1};
        Vec3 unit_l = z_axis.cross(diff).normalize();

        // lateral方向の変位
        dl = (ext1 - ext0).dot(unit_l);

        return dl;
    }

    template <class Track1, class Track2>
    double MinimumDistance(const Track1& t1, const Track2& t2)
    {
        double md, coef0, coef1, denom;

        // 各trackの位置ベクトルと方向ベクトル
        Vec3 posV0 = {t1.x, t1.y, t1.z};
        Vec3 posV1 = {t2.x, t2.y, t2.z};
        Vec3 dirV0 = {t1.ax, t1.ay, 1};
        Vec3 dirV1 = {t2.ax, t2.ay, 1};

        // 差分
        Vec3 diff = posV1 - posV0;

        if (OpeningAngle(t1, t2) < EPS) { // 平行な場合
            // trackの方向ベクトルに垂直なdiffベクトルの成分を取り出す
            Vec3 perp_component = diff - dirV0 * diff.dot(dirV0);
            return perp_component.length();
        } else {
            denom = dirV0.dot(dirV0) * dirV1.dot(dirV1) - dirV0.dot(dirV1) * dirV0.dot(dirV1);
            if (denom < EPS) {
                // 平行な場合の処理と同じ
                Vec3 perp_component = diff - dirV0 * diff.dot(dirV0);
                return perp_component.length();
            } else {
                coef0 = (diff.dot(dirV0) * dirV1.dot(dirV1) - diff.dot(dirV1) * dirV0.dot(dirV1)) / denom;
                coef1 = (diff.dot(dirV0) * dirV0.dot(dirV1) - diff.dot(dirV1) * dirV0.dot(dirV0)) / denom;
            }
        }

        double z_low = std::min(posV0.z, posV1.z);
        double z_high = std::max(posV0.z, posV1.z);

        if (posV0.z + coef0 < z_low || posV1.z + coef1 < z_low) {
            coef0 = z_low - posV0.z;
            coef1 = z_low - posV1.z;
        } else if (posV0.z + coef0 > z_high || posV1.z + coef1 > z_high) {
            coef0 = z_high - posV0.z;
            coef1 = z_high - posV1.z;
        }

        // 外挿した位置ベクトル (最近接点)
        Vec3 intercept0 = posV0 + dirV0 * coef0;
        Vec3 intercept1 = posV1 + dirV1 * coef1;

        md = (intercept1 - intercept0).length();

        return md;
    }
}
