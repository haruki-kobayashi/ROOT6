#pragma once

#include <TColor.h>

namespace MyPalette {
    enum class Palette {
        kBirdDark = 114,
        kBlueWhiteRed = 115,
        kBlueBlackRed = 116,
        kGreenWhiteMagenta = 117,
        kGreenBlackMagenta = 118,
        kLegacy = 119,
    };

    struct PaletteData {
        size_t size = 0;
        std::vector<double> Red;
        std::vector<double> Green;
        std::vector<double> Blue;
        std::vector<double> Stops;
    };

    int CreateCustomPalette(PaletteData& data);
    void InvertPalette();
    void NegatePalette();
    void SetPalette(int palette, uint32_t NContours = 256U);
    void SetPalette(std::string palette, uint32_t NContours = 256U);
    void SetPalette(Palette palette, uint32_t NContours = 256U);
    void SetPalette(EColorPalette palette, uint32_t NContours = 256U);
}