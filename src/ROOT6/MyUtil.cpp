#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <TH1.h>
#include <TBox.h>
#include <TStyle.h>
#include <TPad.h>
#include <ROOT6/MyUtil.hpp>

namespace MyUtil {

    void PaintBin(TH1* h, int bin, int color, float alpha) noexcept {
        TBox* b = new TBox(
            h->GetBinLowEdge(bin),
            h->GetMinimum(),
            h->GetBinWidth(bin) + h->GetBinLowEdge(bin),
            h->GetBinContent(bin)
        );
        b->SetFillColorAlpha(color, alpha);
        b->Draw();
    }

    void PaintBins(TH1* h, double h_range_min, double h_range_max, float alpha) noexcept {
        int min = h->FindFixBin(h_range_min);
        int max = h->FindFixBin(h_range_max);

        for (int i = min; i <= max; ++i) {
            if (h->GetBinContent(i) == 0) continue;

            int nColors = gStyle->GetNumberOfColors();

            int ci = gStyle->GetColorPalette(
                nColors * (h->GetBinCenter(i) - h->GetBinLowEdge(min)) /
                (h->GetBinLowEdge(max + 1) - h->GetBinLowEdge(min))
            );

            if (h->GetBinLowEdge(i) < h_range_min) {
                TBox* b = new TBox(
                    h_range_min,
                    h->GetMinimum(),
                    h->GetBinWidth(i) + h->GetBinLowEdge(i),
                    std::min(h->GetBinContent(i), gPad->GetY2())
                );
                b->SetFillColorAlpha(ci, alpha);
                b->Draw();
            } else if (h->GetBinWidth(i) + h->GetBinLowEdge(i) > h_range_max) {
                TBox* b = new TBox(
                    h->GetBinLowEdge(i),
                    h->GetMinimum(),
                    h_range_max,
                    std::min(h->GetBinContent(i), gPad->GetY2())
                );
                b->SetFillColorAlpha(ci, alpha);
                b->Draw();
            } else {
                PaintBin(h, i, ci, alpha);
            }
        }
    }

    void ShowProgress(int &already_done, const double progress) noexcept
    {
        if (progress < 1.0) already_done += 1;
        constexpr uint8_t barWidth = 30;
        uint8_t bar_length = static_cast<uint8_t>(barWidth * progress);

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << (progress * 100);
        std::string percentage = oss.str() + " %";

        std::cout << " Progress [";

        for (uint8_t i = 0; i < barWidth; ++i)
        {
            if (i < bar_length) std::cout << "=";
            else if (i == bar_length) std::cout << ">";
            else std::cout << " ";
        }

        std::cout << "] " << percentage << "\r";
        std::cout.flush();
    }

    void ShowProgress(int &already_done, const int total) noexcept
    {
        double progress = static_cast<double>(already_done) / static_cast<double>(total);
        if (progress < 1.0) already_done += 1;
        constexpr uint8_t barWidth = 30;
        uint8_t bar_length = static_cast<uint8_t>(barWidth * progress);

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << (progress * 100);
        std::string percentage = oss.str() + " %";

        std::cout << " Progress [";

        for (uint8_t i = 0; i < barWidth; ++i)
        {
            if (i < bar_length) std::cout << "=";
            else if (i == bar_length) std::cout << ">";
            else std::cout << " ";
        }

        std::cout << "] " << percentage << "\r";
        std::cout.flush();
    }
}