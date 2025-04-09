#include "MyPalette.hpp"
#include <iostream>
#include <magic_enum.hpp>
#include <TColor.h>
#include <TStyle.h>
#include <variant>

namespace MyPalette {

    static void RedWhiteBlue(uint32_t NContours) {
        const uint32_t NColors = 256;
        const uint32_t Number = 9;
        int Palette[NColors];
        double Red[Number]    = { 31./255.,  71./255., 123./255., 160./255., 255./255., 222./255., 214./255., 199./255., 163./255.};
        double Green[Number]  = { 40./255., 117./255., 171./255., 211./255., 255./255., 190./255., 150./255.,  92./255.,  20./255.};
        double Blue[Number]   = {234./255., 214./255., 228./255., 222./255., 255./255., 160./255.,  95./255.,  60./255.,  34./255.};
        double Stops[Number]  = {0./8., 1./8., 2./8., 3./8., 4./8., 5./8., 6./8., 7./8., 8./8.};
        int FI = TColor::CreateGradientColorTable(Number, Stops, Red, Green, Blue, NColors);
        for (int i = 0; i < NColors; i++) {
            Palette[i] = FI + i;
        }
        gStyle->SetPalette(NColors, Palette);
        gStyle->SetNumberContours(NContours);
    }

    static void BirdDark(uint32_t NContours) {
        const uint32_t NColors = 256;
        const uint32_t Number = 8;
        int Palette[NColors];
        double Red[Number]    = {0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956};
        double Green[Number]  = {0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862};
        double Blue[Number]   = {0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968};
        double Stops[Number]  = {0./7., 1./7., 2./7., 3./7., 4./7., 5./7., 6./7., 7./7.};
        int FI = TColor::CreateGradientColorTable(Number, Stops, Red, Green, Blue, NColors);
        for (int i = 0; i < NColors; i++) {
            Palette[i] = FI + i;
        }
        gStyle->SetPalette(NColors, Palette);
        gStyle->SetNumberContours(NContours);
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
                std::cerr << "Unknown palette name: " << paletteName << std::endl;
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
            std::cerr << "Unknown palette type!" << std::endl;
            return;
        }

        // パレットを設定
        if (paletteIndex >= 114) {
            switch (static_cast<Palette>(paletteIndex)) {
                case Palette::kRedWhiteBlue:
                    RedWhiteBlue(NContours);
                    break;
                case Palette::kBirdDark:
                    BirdDark(NContours);
                    break;
                default:
                    std::cerr << "Unknown custom palette!" << std::endl;
                    break;
            }
        } else {
            gStyle->SetPalette(paletteIndex);
            gStyle->SetNumberContours(NContours);
        }
    }
}
