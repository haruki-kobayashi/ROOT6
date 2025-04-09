#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "ShowProgress.hpp"

void ShowProgress(int &page, double progress)
{
    if (progress < 1.0) page += 1;
    uint8_t barWidth = 30;
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
