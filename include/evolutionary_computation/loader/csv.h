#pragma once

#include <evolutionary_computation/loader/base.h>

class CSVLoader : public FileLoader {
public:
    using FileLoader::FileLoader;

    std::vector<Data> load() override;
};
