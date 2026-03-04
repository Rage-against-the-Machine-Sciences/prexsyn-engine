#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "molecule.hpp"
#include "reaction.hpp"

namespace prexsyn::synthesis {

class SynthesisError : public std::runtime_error {
public:
    explicit SynthesisError(const std::string &message) : std::runtime_error(message) {}
};

class SynthesisNode {
private:
    struct Item {
        std::variant<std::shared_ptr<Molecule>, ReactionOutcomeWithReactantAssignment> m;
        std::vector<size_t> precursor_indices;
    };

    std::vector<Item> items_;
    std::shared_ptr<Reaction> reaction_;
    std::vector<std::pair<std::shared_ptr<SynthesisNode>, size_t>> precursors_;

public:
    SynthesisNode(std::shared_ptr<Molecule> mol) { items_.push_back({mol, {}}); };

    size_t size() const;
    const Item &at(size_t) const;
};

} // namespace prexsyn::synthesis
