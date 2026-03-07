#include "synthesis.hpp"

#include <exception>
#include <string>

#include "bb_lib.hpp"
#include "chemical_space.hpp"
#include "rxn_lib.hpp"

namespace prexsyn::chemspace {

ChemicalSpaceSynthesis::Result
ChemicalSpaceSynthesis::add_building_block(BuildingBlockLibrary::Index index) noexcept {
    try {
        const auto &bb_item = cs_.bb_lib().get(index);
        synthesis_->push(bb_item.molecule);
        postfix_notation_.append(bb_item.index, PostfixNotation::Token::Type::BuildingBlock);
        return Result::ok();
    } catch (const std::exception &e) {
        return Result::error(e.what());
    }
}

ChemicalSpaceSynthesis::Result
ChemicalSpaceSynthesis::add_building_block(const std::string &index) noexcept {
    try {
        const auto &bb_item = cs_.bb_lib().get(index);
        synthesis_->push(bb_item.molecule);
        postfix_notation_.append(bb_item.index, PostfixNotation::Token::Type::BuildingBlock);
        return Result::ok();
    } catch (const std::exception &e) {
        return Result::error(e.what());
    }
}

ChemicalSpaceSynthesis::Result
ChemicalSpaceSynthesis::add_reaction(ReactionLibrary::Index index) noexcept {
    try {
        const auto &rxn_item = cs_.rxn_lib().get(index);
        synthesis_->push(rxn_item.reaction);
        postfix_notation_.append(rxn_item.index, PostfixNotation::Token::Type::Reaction);
        return Result::ok();
    } catch (const std::exception &e) {
        return Result::error(e.what());
    }
}

ChemicalSpaceSynthesis::Result
ChemicalSpaceSynthesis::add_reaction(const std::string &index) noexcept {
    try {
        const auto &rxn_item = cs_.rxn_lib().get(index);
        synthesis_->push(rxn_item.reaction);
        postfix_notation_.append(rxn_item.index, PostfixNotation::Token::Type::Reaction);
        return Result::ok();
    } catch (const std::exception &e) {
        return Result::error(e.what());
    }
}

} // namespace prexsyn::chemspace
