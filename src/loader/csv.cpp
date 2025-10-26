#include "evolutionary_computation/data.h"
#include <sstream>
#include <string>
#include <evolutionary_computation/loader/csv.h>
#include <vector>

Data CSVLoader::load() {
    auto result = std::vector<DataEntry>();

    std::string line;
    while (std::getline(this->file, line)) {
        auto lineStream = std::stringstream(line);
        DataEntry data;
        std::string cell;

        std::getline(lineStream, cell, ';');
        data.x = std::stoi(cell);

        std::getline(lineStream, cell, ';');
        data.y = std::stoi(cell);

        std::getline(lineStream, cell, ';');
        data.cost = std::stoi(cell);

        result.push_back(data);
    }

    return Data { result };
}
