#pragma once

#include <filesystem>

#include "Types.h"

namespace simrail::geometry {
    /**
     * Load a 2D profile from file (Y,Z coordinate pairs).
     *
     * File format expectations:
     * - Each line: "y,z" (comma-separated doubles)
     * - Lines starting with '#' or whitespace+# are comments (skipped)
     * - Empty lines are ignored
     * - UTF-8 BOM is handled automatically
     *
     * @param filepath Path to profile file
     * @return         Profile with validated Y, Z data
     * @throws std::runtime_error if file cannot be opened, parsed, or validation fails
     */
    [[nodiscard]] Profile loadProfileFromFile(const std::filesystem::path &filepath);

    /**
     * Validate a profile for use with PCHIP interpolation.
     *
     * Requirements:
     * - Minimum 2 points
     * - Y values must be strictly increasing
     * - No NaN or Inf values
     *
     * @param profile   Profile to validate
     * @throws std::invalid_argument if validation fails
     */
    void validateProfile(const Profile &profile);
} // namespace simrail::geometry
