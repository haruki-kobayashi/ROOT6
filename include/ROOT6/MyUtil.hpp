#pragma once

#include <TH1.h>

namespace MyUtil {
    void PaintBin(TH1* h, int bin, int color, float alpha = 1.0f) noexcept;
    void PaintBins(TH1* h, double h_range_min, double h_range_max, float alpha = 1.0f) noexcept;
    void ShowProgress(int &page, double progress) noexcept;
    void ShowProgress(int &page, int total) noexcept;
}