#pragma once

#include <fstream>
#include <vector>
#include <evolutionary_computation/data.h>

class FileLoader {
public:
    FileLoader(std::ifstream& file);

    virtual std::vector<Data> load() = 0;
protected:
    std::ifstream& file;
};
