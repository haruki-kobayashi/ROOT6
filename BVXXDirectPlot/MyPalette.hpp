#pragma once

#include <TColor.h>
#include <variant>

namespace MyPalette {
    enum class Palette {
        kRedWhiteBlue = 114,
        kBirdDark = 115,
    };

    void SetPalette(std::variant<int, std::string, Palette, EColorPalette> palette, uint32_t NContours = 256U);
    void RedWhiteBlue(uint32_t NContours);
    void BirdDark(uint32_t NContours);
}