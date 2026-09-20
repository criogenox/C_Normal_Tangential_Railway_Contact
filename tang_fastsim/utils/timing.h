#pragma once

#include <chrono>
#include <iostream>

namespace fastsim::timing {
    inline auto tick_now() {
        return std::chrono::high_resolution_clock::now();
    }

    inline void bench(auto startTime, auto endTime) {
        const std::chrono::duration<double, std::micro> elapsed{endTime - startTime};
        std::cout << "\033[0;38;5;144mComputation time elapsed: "
                << elapsed.count() << " microseconds" << std::endl;
    }
}
