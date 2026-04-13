#include "guided.hpp"

#include <cstddef>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../chemistry/chemistry.hpp"
#include "../chemspace/chemspace.hpp"
#include "../utility/random.hpp"

namespace prexsyn::enumerator {

GuidedEnumerator::GuidedEnumerator(std::shared_ptr<chemspace::ChemicalSpace> cs,
                                   std::vector<double> raw_bb_weights, double smoothing_alpha,
                                   const Config &config, std::optional<size_t> random_seed)
    : cs_(std::move(cs)), config_(config), synthesis_(nullptr),
      rng_(random_seed.value_or(std::random_device{}())) {

    if (cs_ == nullptr)
        throw std::invalid_argument("null pointer for chemical space");

    const size_t n = cs_->bb_lib().size();
    if (n == 0)
        throw std::invalid_argument("empty building block library");
    if (cs_->rxn_lib().size() == 0)
        throw std::invalid_argument("empty reaction library");
    if (!raw_bb_weights.empty() && raw_bb_weights.size() != n)
        throw std::invalid_argument("raw_bb_weights size must equal bb_lib size");

    bb_weights_.resize(n);
    for (size_t i = 0; i < n; ++i) {
        double raw = raw_bb_weights.empty() ? 0.0 : raw_bb_weights[i];
        bb_weights_[i] = raw + smoothing_alpha;
    }

    bb_dist_ = std::discrete_distribution<size_t>(bb_weights_.begin(), bb_weights_.end());
}

size_t GuidedEnumerator::weighted_sample_from(const std::vector<size_t> &candidates) {
    std::vector<double> w;
    w.reserve(candidates.size());
    for (size_t idx : candidates)
        w.push_back(bb_weights_[idx]);
    std::discrete_distribution<size_t> dist(w.begin(), w.end());
    return candidates[dist(rng_)];
}

bool GuidedEnumerator::not_growable() const {
    if (synthesis_ == nullptr)
        return false;
    auto products = synthesis_->products();
    if (products.empty())
        return false;
    const auto &product = products.front();
    return product->num_heavy_atoms() >= config_.heavy_atom_limit ||
           synthesis_->count_building_blocks() >= config_.max_building_blocks;
}

void GuidedEnumerator::clear_synthesis() { synthesis_.reset(); }

void GuidedEnumerator::init_synthesis() {
    synthesis_ = cs_->new_synthesis();
    size_t index = bb_dist_(rng_); // ← weighted instead of uniform
    synthesis_->add_building_block(index);
}

void GuidedEnumerator::grow_synthesis() {
    auto products = synthesis_->products();
    if (products.empty()) {
        clear_synthesis();
        return;
    }
    auto product = random_choice(products, rng_);

    auto matches = cs_->rxn_lib().match_reactants(*product);
    std::erase_if(matches,
                  [&](const auto &match) { return match.count > config_.selectivity_cutoff; });
    if (matches.empty()) {
        clear_synthesis();
        return;
    }
    auto match = random_choice(matches, rng_);

    auto rxn = cs_->rxn_lib().get(match.reaction_index);
    chemspace::Synthesis::Result result;

    for (size_t i = 0; i < rxn.reaction->num_reactants(); ++i) {
        if (i == match.reactant_index)
            continue;

        const auto &rlist_bb = cs_->building_block_reactant_lists().get(match.reaction_index, i);
        const auto &rlist_int = cs_->intermediate_reactant_lists().get(match.reaction_index, i);

        if (rlist_bb.empty() && rlist_int.empty()) {
            clear_synthesis();
            return;
        }

        if (rlist_bb.empty()) {
            // no BB candidates, pick intermediate uniformly (unchanged from random)
            std::uniform_int_distribution<size_t> u(0, rlist_int.size() - 1);
            result =
                synthesis_->add_intermediate(rlist_int[u(rng_)], config_.max_outcomes_per_reaction);
        } else if (rlist_int.empty()) {
            // only BBs available — weighted pick
            result = synthesis_->add_building_block(weighted_sample_from(rlist_bb));
        } else {
            // both available: weight total_bb vs total_int by same scale
            double total_bb = 0.0;
            for (size_t idx : rlist_bb)
                total_bb += bb_weights_[idx];
            double avg_w = total_bb / rlist_bb.size();
            double total_int = avg_w * rlist_int.size();

            std::uniform_real_distribution<double> coin(0.0, total_bb + total_int);
            if (coin(rng_) < total_bb) {
                result = synthesis_->add_building_block(weighted_sample_from(rlist_bb));
            } else {
                std::uniform_int_distribution<size_t> u(0, rlist_int.size() - 1);
                result = synthesis_->add_intermediate(rlist_int[u(rng_)],
                                                      config_.max_outcomes_per_reaction);
            }
        }

        if (!result) {
            clear_synthesis();
            return;
        }
    }

    result = synthesis_->add_reaction(match.reaction_index, config_.max_outcomes_per_reaction);
    if (!result) {
        clear_synthesis();
        return;
    }
}

std::optional<std::shared_ptr<chemspace::Synthesis>> GuidedEnumerator::try_next() {
    if (synthesis_ == nullptr || not_growable())
        init_synthesis();
    else
        grow_synthesis();

    if (synthesis_ == nullptr)
        return std::nullopt;
    return std::make_shared<chemspace::Synthesis>(*synthesis_);
}

std::shared_ptr<chemspace::Synthesis> GuidedEnumerator::next() {
    constexpr int max_attempts = 100;
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        auto result = try_next();
        if (result.has_value())
            return result.value();
    }
    throw std::runtime_error("too many failed attempts in GuidedEnumerator");
}

std::pair<std::shared_ptr<chemspace::Synthesis>, std::shared_ptr<Molecule>>
GuidedEnumerator::next_with_product() {
    auto syn = next();
    auto product = random_choice(syn->products(), rng_);
    return {syn, product};
}

} // namespace prexsyn::enumerator
