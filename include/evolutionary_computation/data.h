#pragma once

#include <optional>
#include <string>
#include <vector>

struct DataEntry {
    int x;
    int y;
    int cost;
    std::optional<std::string> filename;
};

struct Data {
    std::vector<DataEntry> entries;
    std::optional<std::string> filename;
};
