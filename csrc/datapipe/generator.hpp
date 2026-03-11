#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <random>

#include "../chemspace/chemspace.hpp"

namespace prexsyn::datapipe {

class Generator {
public:
    struct Config {
        unsigned int max_building_blocks;
        unsigned int heavy_atom_limit;
    };
    constexpr const static Config default_config = {
        .max_building_blocks = 5,
        .heavy_atom_limit = 50,
    };

private:
    std::shared_ptr<chemspace::ChemicalSpace> cs_;
    Config config_;

    std::unique_ptr<chemspace::Synthesis> synthesis_;
    std::mt19937 rng_;

    bool is_limit_exceeded() const;

    void clear_synthesis();
    void init_synthesis();
    void grow_synthesis();
    std::optional<std::shared_ptr<chemspace::Synthesis>> try_next();

public:
    Generator(std::shared_ptr<chemspace::ChemicalSpace> cs, const Config &config = default_config,
              size_t random_seed = std::random_device{}());

    std::shared_ptr<chemspace::Synthesis> next();
};

} // namespace prexsyn::datapipe
