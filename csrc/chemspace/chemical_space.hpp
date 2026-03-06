#pragma once

#include <memory>

#include "building_block.hpp"
#include "reaction.hpp"

namespace prexsyn::chemspace {

class ChemicalSpace {
private:
    std::unique_ptr<BuildingBlockLibrary> bb_lib_;
    std::unique_ptr<ReactionLibrary> rxn_lib_;

public:
    const BuildingBlockLibrary &bb_lib() const { return *bb_lib_; }
    const ReactionLibrary &rxn_lib() const { return *rxn_lib_; }
};

} // namespace prexsyn::chemspace
