#include "yaml-cpp/yaml.h"
#include <iostream>

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "usage: " << argv[0] << " path/to/file.yaml" << std::endl;
        return -1;
    }

    // load test config
    YAML::Node config = YAML::LoadFile(argv[1]);
    std::cout << config << std::endl;

    std::cout << "name " << config["name"].as<std::string>() << std::endl;
    std::cout << "mass " << config["inertia"]["mass"].as<double>() << std::endl;

    auto inertia = config["inertia"];
    std::cout << "inertia " << inertia << std::endl;
}
