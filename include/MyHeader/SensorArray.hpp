#pragma once

#include <array>

namespace SensorArray{

    constexpr int calcValue(int x, int y) {
        // x: (0-7), y: (0-8)
        // xが偶数と奇数で共通のオフセット
        int row_offset = (x % 2 == 0) ? x / 2 : (x - 1) / 2;
        int g = y / 3; // グループ番号
        int r = y % 3; // グループ内の位置

        if (x % 2 == 0) { // xが偶数の場合
            if (r == 0) return 12 + 4 * g + row_offset;
            else if (r == 1) return 44 - 4 * g + row_offset;
            else return 68 - 4 * g + row_offset;
        } else { // xが奇数の場合
            if (r == 0) return 8 - 4 * g + row_offset;
            else if (r == 1) return 24 + 4 * g + row_offset;
            else return 48 + 4 * g + row_offset;
        }
    }

    constexpr std::array<std::array<int, 9>, 8> generateSensorArray() {
        std::array<std::array<int, 9>, 8> arr{};

        for (int x = 0; x < 8; ++x)
            for (int y = 0; y < 9; ++y)
                arr[x][y] = calcValue(x, y);

        return arr;
    }

    constexpr auto SensorArray = generateSensorArray();

    constexpr int kMaxValue = 71;
    constexpr int kInvalid = -1;

    struct Index {
        int x = kInvalid;
        int y = kInvalid;

        constexpr Index() = default;
        constexpr Index(int x_, int y_) : x(x_), y(y_) {}
    };

    // 逆引きテーブル
    constexpr auto generateReverseMap() {
        std::array<Index, kMaxValue + 1> map{};

        for (int x = 0; x < 8; ++x)
            for (int y = 0; y < 9; ++y) {
                int val = SensorArray[x][y];
                map[val] = Index{x, y};
            }
        return map;
    }

    constexpr auto ReverseMap = generateReverseMap();

    // センサー番号からx, yのインデックスを取得する関数
    inline Index GetIndex(int value) {
        if (value < 0 || value > kMaxValue) return Index{};
        return ReverseMap[value];
    }
}