#pragma once

#include <TColor.h>

namespace MyPalette {
    enum class Palette {
        kRedWhiteBlue = 114,
        kBirdDark = 115,
    };

    void SetPalette(int palette, UInt_t NContours = 256, Float_t alpha = 1.0);
    void RedWhiteBlue(UInt_t NContours, Float_t alpha);
    void BirdDark(UInt_t NContours, Float_t alpha);
}