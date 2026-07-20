#ifndef GINIBRE_Q3_FRONTIER_RUNTIME_SCOPE_HPP
#define GINIBRE_Q3_FRONTIER_RUNTIME_SCOPE_HPP

#include <cstdlib>
#include <stdexcept>
#include <string>

inline bool frontier_environment_flag(const char* name, bool default_value) {
    const char* raw = std::getenv(name);
    if (raw == nullptr) return default_value;
    const std::string value(raw);
    if (value == "1") return true;
    if (value == "0") return false;
    throw std::runtime_error(
        std::string("environment flag ") + name + " must be 0 or 1"
    );
}

inline bool frontier_run_b() {
    return frontier_environment_flag("RUN_B", true);
}

#endif
