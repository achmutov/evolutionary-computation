#pragma once

#include <evolutionary_computation/loader/base.h>
#include <optional>

class CSVLoader : public FileLoader {
public:
    using FileLoader::FileLoader;
    std::optional<std::string> filename;

    CSVLoader(std::ifstream& file, std::string filename) : FileLoader(file), filename{filename} {}

    Data load() override;
};
