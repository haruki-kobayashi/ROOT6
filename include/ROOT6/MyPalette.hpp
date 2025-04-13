#pragma once

#include <TColor.h>
#include <variant>

namespace MyPalette {
    enum class Palette {
        kBirdDark = 114,
        kRedWhiteBlue = 115,
        kRedBlackBlue = 116,
        kMagentaWhiteGreen = 117,
        kMagentaBlackGreen = 118,
        kLegacy = 119,
    };

    struct PaletteData {
        uint32_t Number = 0;
        double* Red = nullptr;
        double* Green = nullptr;
        double* Blue = nullptr;
        double* Stops = nullptr;
        uint32_t NContours = 256U;
    };

    void SetPalette(std::variant<int, std::string, Palette, EColorPalette> palette, uint32_t NContours = 256U);
    void InvertPalette();
    int CreateCustomPalette(PaletteData& data);
}