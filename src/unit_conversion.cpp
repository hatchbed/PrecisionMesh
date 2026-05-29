#include <precision_mesh/unit_conversion.h>

#include <unordered_map>

#include <boost/algorithm/string/case_conv.hpp>

static const std::unordered_map<std::string, float> to_meters = {
    {"mm", 0.001f},
    {"cm", 0.01f},
    {"m",  1.0f},
    {"in", 0.0254f},
    {"ft", 0.3048f},
    {"yd", 0.9144f}
};

std::string normalizeUnit(const std::string& unit) {
    if (unit.empty()) {
        return "mesh units";
    }

    auto lower = boost::algorithm::to_lower_copy(unit);
    if (lower == "millimetre" || lower == "millimetres" || lower == "millimeter" ||
        lower == "millimeters" || lower == "mm") {
        return "mm";
    }
    if (lower == "centimetre" || lower == "centimetres" || lower == "centimeter" ||
        lower == "centimeters" || lower == "cm") {
        return "cm";
    }
    if (lower == "metre" || lower == "metres" || lower == "meter" || lower == "meters" ||
        lower == "m") {
        return "m";
    }
    if (lower == "inch" || lower == "inches" || lower == "in") {
        return "in";
    }
    if (lower == "foot" || lower == "feet" || lower == "ft") {
        return "ft";
    }
    if (lower == "yard" || lower == "yards" || lower == "yd") {
        return "yd";
    }

    return lower;
}

bool isKnownUnit(const std::string& unit) {
    return to_meters.count(unit) > 0;
}

float getUnitConversionScale(const std::string& from, const std::string& to) {
    auto from_it = to_meters.find(from);
    auto to_it   = to_meters.find(to);

    if (from_it == to_meters.end() || to_it == to_meters.end()) {
        return 1.0f;
    }

    return from_it->second / to_it->second;
}
