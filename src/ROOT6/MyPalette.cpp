#include <iostream>
#include <TROOT.h>
#include <TColor.h>
#include <TStyle.h>
#include <magic_enum/magic_enum.hpp>
#include <ROOT6/MyPalette.hpp>

namespace MyPalette {

    namespace {
        PaletteData SetPaletteData(Palette palette) {
            PaletteData data;
            switch (palette) {
                case Palette::kBirdDark:
                    data = {
                        9,
                        {0.1652, 0.0367, 0.0524, 0.0138, 0.1427, 0.4800, 0.7400, 0.9418, 0.9270},
                        {0.1320, 0.3228, 0.4686, 0.5862, 0.6206, 0.6600, 0.6400, 0.7482, 0.9224},
                        {0.4684, 0.7826, 0.7547, 0.6757, 0.5211, 0.4300, 0.3200, 0.1582, 0.0457},
                        {0./8., 1./8., 2./8., 3./8., 4./8., 5./8., 6./8., 7./8., 8./8.}
                    };
                    break;
                case Palette::kBlueWhiteRed:
                    data = {
                        10,
                        { 27./255.,  37./255.,   7./255.,  54./255., 255./255., 248./255., 223./255., 194./255., 143./255., 116./255.},
                        { 48./255.,  59./255.,  99./255., 205./255., 255./255., 194./255.,  76./255.,  37./255.,  35./255.,   0./255.},
                        {155./255., 163./255., 255./255., 255./255., 255./255.,  90./255.,  24./255.,  10./255.,   3./255.,   0./255.},
                        {0.00, 0.14, 0.29, 0.46, 0.50, 0.54, 0.73, 0.82, 0.95, 1.00}
                    };
                    break;
                case Palette::kBlueBlackRed:
                    data = {
                        9,
                        { 54./255.,  16./255.,  26./255.,  27./255.,   0./255., 137./255., 210./255., 229./255., 248./255.},
                        {205./255., 120./255.,  74./255.,  45./255.,   0./255.,  35./255.,  58./255., 106./255., 194./255.},
                        {255./255., 255./255., 198./255., 130./255.,   0./255.,   5./255.,  18./255.,  41./255.,  90./255.},
                        {0./8., 1./8., 2./8., 3.3/8., 4./8., 4.7/8., 6./8., 7./8., 8./8.}
                    };
                    break;
                case Palette::kGreenWhiteMagenta:
                    data = {
                        9,
                        {  0./255.,  13./255.,  62./255., 135./255., 255./255., 238./255., 225./255., 199./255., 155./255.},
                        {121./255., 179./255., 212./255., 228./255., 255./255., 143./255.,  53./255.,   9./255.,   0./255.},
                        {  0./255.,  13./255.,  62./255., 135./255., 255./255., 238./255., 225./255., 199./255., 155./255.},
                        {0./8., 1./8., 2./8., 3./8., 4./8., 5./8., 6./8., 7./8., 8./8.}
                    };
                    break;
                case Palette::kGreenBlackMagenta:
                    data = {
                        9,
                        {109./255.,  33./255.,   0./255.,   0./255.,   0./255., 100./255., 155./255., 212./255., 255./255.},
                        {255./255., 192./255., 170./255., 112./255.,   0./255.,   0./255.,   0./255.,  65./255., 119./255.},
                        {109./255.,  33./255.,   0./255.,   0./255.,   0./255., 100./255., 155./255., 212./255., 255./255.},
                        {0./8., 1./8., 2./8., 3./8., 4./8., 5./8., 6./8., 7./8., 8./8.}
                    };
                    break;
                case Palette::kLegacy:
                    data = {
                        20,
                        { 98./255.,  51./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,  34./255., 105./255., 153./255., 225./255., 255./255., 255./255., 255./255., 255./255., 255./255.},
                        {  0./255.,   0./255.,  20./255.,  68./255., 139./255., 187./255., 255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 238./255., 167./255., 119./255.,  47./255.,   0./255.},
                        {255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 252./255., 204./255., 133./255.,  85./255.,  13./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.},
                        {0./19., 1./19., 2./19., 3./19., 4./19., 5./19., 6./19., 7./19., 8./19., 9./19., 10./19., 11./19., 12./19., 13./19., 14./19., 15./19., 16./19., 17./19., 18./19., 19./19.}
                    };
                    break;
                default:
                    std::cerr << "Palette not found." << std::endl;
            }
            return data;
        }
    }

    int CreateCustomPalette(PaletteData& data) {
        if (data.size <= 0) {
            std::cerr << "Error: Invalid vector length specified in PaletteData." << std::endl;
            return 1;
        }
        if (data.Red.empty() || data.Green.empty() || data.Blue.empty() || data.Stops.empty()) {
            std::cerr << "Error: One or more input vectors of custom palette are null." << std::endl;
            return 1;
        }

        // 配列長が一致しているか確認
        if (data.Red.size() != data.size || data.Green.size() != data.size || 
            data.Blue.size() != data.size || data.Stops.size() != data.size) {
            std::cerr << "Error: Vector lengths do not match the specified size."
                      << " Expected size: " << data.size
                      << ", Red: " << data.Red.size()
                      << ", Green: " << data.Green.size()
                      << ", Blue: " << data.Blue.size()
                      << ", Stops: " << data.Stops.size() << std::endl;
            return 1;
        }

        constexpr uint32_t NColors = 256;
        int Palette[NColors];
        int FI = TColor::CreateGradientColorTable(
            data.size, data.Stops.data(), data.Red.data(), data.Green.data(), data.Blue.data(), NColors
        );
        for (int i = 0; i < NColors; i++) {
            Palette[i] = FI + i;
        }
        gStyle->SetPalette(NColors, Palette);

        return 0;
    }

    void InvertPalette() {
        Int_t nColors = gStyle->GetNumberOfColors();
        if (nColors <= 0) {
            std::cerr << "Error: Invalid palette size or no palette is currently set." << std::endl;
        }

        // 現在のパレットを取得
        const TArrayI& currentPaletteArray = TColor::GetPalette();
        const int* currentPalette = currentPaletteArray.GetArray();
        if (!currentPalette) {
            std::cerr << "Error: Failed to retrieve the current palette." << std::endl;
        }

        // 新しい上下反転パレットを作成
        std::vector<int> invertedPalette(nColors);
        for (int i = 0; i < nColors; ++i) {
            invertedPalette[i] = currentPalette[nColors - 1 - i];
        }

        // 反転パレットを設定
        gStyle->SetPalette(nColors, invertedPalette.data());
    }

    void NegatePalette() {
        Int_t nColors = gStyle->GetNumberOfColors();
        if (nColors <= 0) {
            std::cerr << "Error: Invalid palette size or no palette is currently set." << std::endl;
        }

        // 現在のパレットを取得
        const TArrayI& currentPaletteArray = TColor::GetPalette();
        const int* currentPalette = currentPaletteArray.GetArray();
        if (!currentPalette) {
            std::cerr << "Error: Failed to retrieve the current palette." << std::endl;
        }

        // 新しいネガティブパレットを作成
        std::vector<int> negatedPalette(nColors);
        for (int i = 0; i < nColors; ++i) {
            float ri, gi, bi;
            gROOT->GetColor(currentPalette[i])->GetRGB(ri, gi, bi);
            int ci = TColor::GetFreeColorIndex();
            auto color = new TColor(ci, 1.0 - ri, 1.0 - gi, 1.0 - bi);
            negatedPalette[i] = ci;
        }

        // ネガティブパレットを設定
        gStyle->SetPalette(nColors, negatedPalette.data());
    }

    void SetPalette(int palette, uint32_t NContours) {
        if (palette >= 114) {
            PaletteData data = SetPaletteData(static_cast<Palette>(palette));
            CreateCustomPalette(data);
        } else {
            gStyle->SetPalette(palette);
        }
        gStyle->SetNumberContours(NContours);
    }

    void SetPalette(std::string palette, uint32_t NContours) {
        int paletteIndex = -1;

        if (auto PaletteOpt = magic_enum::enum_cast<Palette>(palette)) {
            paletteIndex = static_cast<int>(*PaletteOpt);
        } else if (auto PaletteOpt = magic_enum::enum_cast<EColorPalette>(palette)) {
            paletteIndex = static_cast<int>(*PaletteOpt);
        } else {
            std::cerr << "Unknown palette name: " << palette << std::endl;
            gStyle->SetNumberContours(NContours);
            return;
        }

        if (paletteIndex >= 114) {
            PaletteData data = SetPaletteData(static_cast<Palette>(paletteIndex));
            CreateCustomPalette(data);
        } else {
            gStyle->SetPalette(paletteIndex);
        }
        gStyle->SetNumberContours(NContours);
    }

    void SetPalette(Palette palette, uint32_t NContours) {
        PaletteData data = SetPaletteData(palette);
        CreateCustomPalette(data);
        gStyle->SetNumberContours(NContours);
    }

    void SetPalette(EColorPalette palette, uint32_t NContours) {
        gStyle->SetPalette(palette);
        gStyle->SetNumberContours(NContours);
    }
}
