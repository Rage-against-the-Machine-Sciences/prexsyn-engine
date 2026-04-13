#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <random>
#include <utility>
#include <vector>

#include "../chemistry/chemistry.hpp"
#include "../chemspace/chemspace.hpp"
#include "base.hpp"

namespace prexsyn::enumerator {

class GuidedEnumerator {
public:
    using Config = EnumeratorConfig;

private:
    std::shared_ptr<chemspace::ChemicalSpace> cs_;
    Config config_;
    std::vector<double> bb_weights_;             // smoothed, length == bb_lib().size()
    std::discrete_distribution<size_t> bb_dist_; // for init_synthesis()
    std::unique_ptr<chemspace::Synthesis> synthesis_;
    std::mt19937 rng_;

    bool not_growable() const;

    void clear_synthesis();
    void init_synthesis();
    void grow_synthesis();
    std::optional<std::shared_ptr<chemspace::Synthesis>> try_next();

    // weighted sample from a subset of BB indices
    size_t weighted_sample_from(const std::vector<size_t> &candidates);

public:
    // raw_bb_weights: length == bb_lib().size(), value == corpus count (0 = unseen)
    // smoothing_alpha added to every entry before normalisation
    GuidedEnumerator(std::shared_ptr<chemspace::ChemicalSpace> cs,
                     std::vector<double> raw_bb_weights, double smoothing_alpha = 1.0,
                     const Config &config = kDefaultEnumeratorConfig,
                     std::optional<size_t> random_seed = std::nullopt);

    std::shared_ptr<chemspace::Synthesis> next();
    std::pair<std::shared_ptr<chemspace::Synthesis>, std::shared_ptr<Molecule>> next_with_product();
};

} // namespace prexsyn::enumerator
