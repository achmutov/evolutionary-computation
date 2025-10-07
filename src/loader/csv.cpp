#include <sstream>
#include <string>
#include <evolutionary_computation/loader/csv.h>
#include <vector>

std::vector<Data> CSVLoader::load() {
    auto result = std::vector<Data>();

    std::string line;
    while (std::getline(this->file, line)) {
        auto lineStream = std::stringstream(line);
        Data data;
        std::string cell;

        std::getline(lineStream, cell, ';');
        data.x = std::stoi(cell);

        std::getline(lineStream, cell, ';');
        data.y = std::stoi(cell);

        std::getline(lineStream, cell, ';');
        data.cost = std::stoi(cell);

        result.push_back(data);
    }
    return result;
}
