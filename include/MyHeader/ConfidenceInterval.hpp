#pragma once

#include <utility>
#include <limits>
#include <cmath>

namespace Confidence {

// Wilson score interval (normal approximation method 1)
// z = 1.0 gives ~68% interval (1-sigma)
inline std::pair<double, double> wilson_interval(int n, int x, double z = 1.0)
{
    if (n <= 0) {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }

    const double nn = static_cast<double>(n);
    const double p = static_cast<double>(x) / nn;
    const double z2 = z * z;
    const double denom = 1.0 + z2 / nn;
    const double center = (p + z2 / (2.0 * nn)) / denom;
    const double half = (z / denom) * std::sqrt((p * (1.0 - p) / nn) + (z2 / (4.0 * nn * nn)));

    double lower = center - half;
    double upper = center + half;

    if (lower < 0.0) lower = 0.0;
    if (upper > 1.0) upper = 1.0;
    return {lower, upper};
}

} // namespace Confidence

