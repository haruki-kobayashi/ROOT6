#include <ROOT6/MyPalette.hpp>
#include <iostream>
#include <magic_enum.hpp>
#include <TColor.h>
#include <TStyle.h>
#include <variant>

namespace MyPalette {

    // 各パレット設定関数
    static void RedWhiteBlue(uint32_t NContours) {
        constexpr uint32_t Number = 10;
        static double Red[Number]   = { 27./255.,  37./255.,   7./255.,  54./255., 255./255., 248./255., 223./255., 194./255., 143./255., 116./255.};
        static double Green[Number] = { 48./255.,  59./255.,  99./255., 205./255., 255./255., 194./255.,  76./255.,  37./255.,  35./255.,   0./255.};
        static double Blue[Number]  = {155./255., 163./255., 255./255., 255./255., 255./255.,  90./255.,  24./255.,  10./255.,   3./255.,   0./255.};
        static double Stops[Number] = {0.00, 0.14, 0.29, 0.46, 0.50, 0.54, 0.73, 0.82, 0.95, 1.00};

        PaletteData data = {Number, Red, Green, Blue, Stops, NContours};
        CreateCustomPalette(data);
    }

    static void RedBlackBlue(uint32_t NContours) {
        constexpr uint32_t Number = 9;
        static double Red[Number]   = { 54./255.,  16./255.,  26./255.,  27./255.,   0./255., 137./255., 210./255., 229./255., 248./255.};
        static double Green[Number] = {205./255., 120./255.,  74./255.,  45./255.,   0./255.,  35./255.,  58./255., 106./255., 194./255.};
        static double Blue[Number]  = {255./255., 255./255., 198./255., 130./255.,   0./255.,   5./255.,  18./255.,  41./255.,  90./255.};
        static double Stops[Number] = {0./8., 1./8., 2./8., 3.3/8., 4./8., 4.7/8., 6./8., 7./8., 8./8.};

        PaletteData data = {Number, Red, Green, Blue, Stops, NContours};
        CreateCustomPalette(data);
    }

    static void MagentaWhiteGreen(uint32_t NContours) {
        constexpr uint32_t Number = 9;
        static double Red[Number]   = {  0./255.,  13./255.,  62./255., 135./255., 255./255., 238./255., 225./255., 199./255., 155./255.};
        static double Green[Number] = {121./255., 179./255., 212./255., 228./255., 255./255., 143./255.,  53./255.,   9./255.,   0./255.};
        static double Blue[Number]  = {  0./255.,  13./255.,  62./255., 135./255., 255./255., 238./255., 225./255., 199./255., 155./255.};
        static double Stops[Number] = {0./8., 1./8., 2./8., 3./8., 4./8., 5./8., 6./8., 7./8., 8./8.};

        PaletteData data = {Number, Red, Green, Blue, Stops, NContours};
        CreateCustomPalette(data);
    }

    static void MagentaBlackGreen(uint32_t NContours) {
        constexpr uint32_t Number = 9;
        static double Red[Number]   = {109./255.,  33./255.,   0./255.,   0./255.,   0./255., 100./255., 155./255., 212./255., 255./255.};
        static double Green[Number] = {255./255., 192./255., 170./255., 112./255.,   0./255.,   0./255.,   0./255.,  65./255., 119./255.};
        static double Blue[Number]  = {109./255.,  33./255.,   0./255.,   0./255.,   0./255., 100./255., 155./255., 212./255., 255./255.};
        static double Stops[Number] = {0./8., 1./8., 2./8., 3./8., 4./8., 5./8., 6./8., 7./8., 8./8.};

        PaletteData data = {Number, Red, Green, Blue, Stops, NContours};
        CreateCustomPalette(data);
    }

    static void BirdDark(uint32_t NContours) {
        constexpr uint32_t Number = 8;
        static double Red[Number]   = {0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956};
        static double Green[Number] = {0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862};
        static double Blue[Number]  = {0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968};
        static double Stops[Number] = {0./7., 1./7., 2./7., 3./7., 4./7., 5./7., 6./7., 7./7.};

        PaletteData data = {Number, Red, Green, Blue, Stops, NContours};
        CreateCustomPalette(data);
    }

    static void Legacy() {
        constexpr uint32_t Number = 20;
        static double Red[Number]   = { 98./255.,  51./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,  34./255., 105./255., 153./255., 225./255., 255./255., 255./255., 255./255., 255./255., 255./255.};
        static double Green[Number] = {  0./255.,   0./255.,  20./255.,  68./255., 139./255., 187./255., 255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 238./255., 167./255., 119./255.,  47./255.,   0./255.};
        static double Blue[Number]  = {255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 252./255., 204./255., 133./255.,  85./255.,  13./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.,   0./255.};
        static double Stops[Number] = {0./19., 1./19., 2./19., 3./19., 4./19., 5./19., 6./19., 7./19., 8./19., 9./19., 10./19., 11./19., 12./19., 13./19., 14./19., 15./19., 16./19., 17./19., 18./19., 19./19.};

        PaletteData data = {Number, Red, Green, Blue, Stops, Number};
        CreateCustomPalette(data);
    }

    void SetPalette(std::variant<int, std::string, Palette, EColorPalette> palette, uint32_t NContours) {
        int paletteIndex = -1;

        // palette が int 型の場合
        if (std::holds_alternative<int>(palette)) {
            paletteIndex = std::get<int>(palette);
        }
        // palette が string 型の場合
        else if (std::holds_alternative<std::string>(palette)) {
            const auto& paletteName = std::get<std::string>(palette);
            if (auto PaletteOpt = magic_enum::enum_cast<Palette>(paletteName)) {
                paletteIndex = static_cast<int>(*PaletteOpt);
            } else if (auto PaletteOpt = magic_enum::enum_cast<EColorPalette>(paletteName)) {
                paletteIndex = static_cast<int>(*PaletteOpt);
            } else {
                std::cerr << "\nUnknown palette name: " << paletteName << std::endl;
                return;
            }
        }
        // palette が Palette 型の場合
        else if (std::holds_alternative<Palette>(palette)) {
            paletteIndex = static_cast<int>(std::get<Palette>(palette));
        }
        // palette が EColorPalette 型の場合
        else if (std::holds_alternative<EColorPalette>(palette)) {
            paletteIndex = static_cast<int>(std::get<EColorPalette>(palette));
        } else {
            std::cerr << "\nUnknown palette type!" << std::endl;
            return;
        }

        // パレットを設定
        if (paletteIndex >= 114) {
            switch (static_cast<Palette>(paletteIndex)) {
                case Palette::kRedWhiteBlue:
                    RedWhiteBlue(NContours);
                    break;
                case Palette::kRedBlackBlue:
                RedBlackBlue(NContours);
                    break;
                case Palette::kMagentaWhiteGreen:
                    MagentaWhiteGreen(NContours);
                    break;
                case Palette::kMagentaBlackGreen:
                    MagentaBlackGreen(NContours);
                    break;
                case Palette::kBirdDark:
                    BirdDark(NContours);
                    break;
                case Palette::kLegacy:
                Legacy();
                    break;
                default:
                    std::cerr << "\nUnknown custom palette!" << std::endl;
                    break;
            }
        } else {
            gStyle->SetPalette(paletteIndex);
            gStyle->SetNumberContours(NContours);
        }
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

        // 新しい反転パレットを作成
        std::vector<int> invertedPalette(nColors);
        for (int i = 0; i < nColors; ++i) {
            invertedPalette[i] = currentPalette[nColors - 1 - i];
        }

        // 反転パレットを設定
        gStyle->SetPalette(nColors, invertedPalette.data());
    }

    int CreateCustomPalette(PaletteData& data) {
        if (data.Number <= 0) {
            std::cerr << "Error: Invalid array length specified in PaletteData." << std::endl;
            return 1;
        }
        if (!data.Red || !data.Green || !data.Blue || !data.Stops) {
            std::cerr << "Error: One or more input arrays of custom palette are null." << std::endl;
            return 1;
        }

        // 配列長が一致しているか確認
        for (uint32_t i = 0; i < data.Number; ++i) {
            if (!data.Red || !data.Green || !data.Blue || !data.Stops) {
                std::cerr << "Error: Array lengths of Red, Green, Blue, or Stops do not match the specified Number." << std::endl;
                return 1;
            }
        }

        const uint32_t NColors = 256;
        int Palette[NColors];
        int FI = TColor::CreateGradientColorTable(data.Number, data.Stops, data.Red, data.Green, data.Blue, NColors);
        for (int i = 0; i < NColors; i++) {
            Palette[i] = FI + i;
        }
        gStyle->SetPalette(NColors, Palette);
        gStyle->SetNumberContours(data.NContours);

        return 0;
    }
}
