#include "ProfileLoader.h"

#include <cmath>
#include <fstream>

#include "Constants.h"

namespace simrail::geometry {
    Profile loadProfileFromFile(const std::filesystem::path &filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filepath.string());
        }

        Profile p;
        std::string line;
        int line_number = 0;
        int valid_points = 0;

        // Skip UTF-8 BOM if present (EF BB BF)
        char bom_check[3];
        if (file.read(bom_check, 3).gcount() == 3) {
            bool has_bom = (static_cast<unsigned char>(bom_check[0]) == 0xEF &&
                            static_cast<unsigned char>(bom_check[1]) == 0xBB &&
                            static_cast<unsigned char>(bom_check[2]) == 0xBF);
            if (!has_bom) {
                file.seekg(0, std::ios::beg);
            }
        } else {
            file.seekg(0, std::ios::beg);
        }

        while (std::getline(file, line)) {
            ++line_number;

            // Skip empty lines
            if (line.empty()) {
                continue;
            }

            // Trim leading whitespace
            auto start = line.find_first_not_of(" \t\r\n");
            if (start == std::string::npos) {
                continue; // Whitespace-only line
            }

            // Skip comment lines (after trimming)
            if (line[start] == '#') {
                continue;
            }

            // Parse Y,Z values
            std::stringstream ss(line.substr(start));
            char comma;
            double y, z;

            if (!(ss >> y >> comma >> z)) {
                throw std::runtime_error("Failed to parse line " + std::to_string(line_number) +
                                         ": '" + line + "'");
            }

            // Verify comma separator
            if (comma != ',') {
                throw std::runtime_error("Expected comma separator at line " +
                                         std::to_string(line_number) + ": '" + line + "'");
            }

            // Check for extra data after Z
            char extra;
            if (ss >> extra) {
                throw std::runtime_error("Extra data after Z value at line " +
                                         std::to_string(line_number) + ": '" + line + "'");
            }

            // Validate numeric values (no NaN or Inf)
            if (std::isnan(y) || std::isinf(y) || std::isnan(z) || std::isinf(z)) {
                throw std::runtime_error("Invalid numeric value (NaN or Inf) at line " +
                                         std::to_string(line_number));
            }

            p.y.push_back(y);
            p.z.push_back(z);
            ++valid_points;
        }

        file.close();

        // Validate minimum points
        if (valid_points < 2) {
            throw std::runtime_error("Profile must have at least 2 points (found " +
                                     std::to_string(valid_points) + ") in: " + filepath.string());
        }

        // Validate Y values are strictly increasing
        for (size_t i = 1; i < p.y.size(); ++i) {
            const double dy = p.y[i] - p.y[i - 1];
            if (dy <= constants::FLOATING_POINT_TOL) {
                throw std::runtime_error("Y values must be strictly increasing at index " +
                                         std::to_string(i) + " (duplicate or decreasing value: " +
                                         std::to_string(p.y[i]) + ") in: " + filepath.string());
            }
        }

        return p;
    }

    void validateProfile(const Profile &profile) {
        // Minimum points check
        if (profile.size() < 2) {
            throw std::invalid_argument("Profile must have at least 2 points");
        }

        // Check for NaN/Inf in both Y and Z
        for (size_t i = 0; i < profile.size(); ++i) {
            if (std::isnan(profile.y[i]) || std::isinf(profile.y[i])) {
                throw std::invalid_argument("Invalid Y value at index " + std::to_string(i));
            }
            if (std::isnan(profile.z[i]) || std::isinf(profile.z[i])) {
                throw std::invalid_argument("Invalid Z value at index " + std::to_string(i));
            }
        }

        // Check strictly increasing Y
        for (size_t i = 1; i < profile.size(); ++i) {
            const double dy = profile.y[i] - profile.y[i - 1];
            if (dy <= constants::FLOATING_POINT_TOL) {
                throw std::invalid_argument("Y values must be strictly increasing at index " +
                                            std::to_string(i));
            }
        }

        // Check Y and Z size match
        if (profile.y.size() != profile.z.size()) {
            throw std::invalid_argument("Y and Z vectors must have equal size");
        }
    }
} // namespace simrail::geometry
