#pragma once

#include <string>

// Normalize a long-form unit name (e.g. "millimetre", "Millimeters") to its
// canonical short form ("mm", "cm", "m", "in", "ft", "yd").  Returns the
// lower-cased input unchanged if unrecognized, or "mesh units" if empty.
std::string normalizeUnit(const std::string& unit);

// Returns true if unit is a known, convertible unit (one of mm/cm/m/in/ft/yd).
bool isKnownUnit(const std::string& unit);

// Returns the scale factor to convert a length from `from` to `to`.
// Returns 1.0 if either unit is unknown.
float getUnitConversionScale(const std::string& from, const std::string& to);
