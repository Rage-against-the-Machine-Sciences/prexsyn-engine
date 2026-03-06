#include "reaction.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace prexsyn::chemspace {

const ReactionItem &ReactionLibrary::get(Index index) const {
    if (index >= reactions_.size()) {
        throw std::out_of_range("Reaction index out of range");
    }
    return reactions_[index];
}

const ReactionItem &ReactionLibrary::get(const std::string &name) const {
    auto it = name_to_index_.find(name);
    if (it == name_to_index_.end()) {
        throw std::out_of_range("Reaction name not found: " + name);
    }
    return reactions_[it->second];
}

ReactionLibrary::Index ReactionLibrary::add(const ReactionEntry &entry) {
    if (name_to_index_.contains(entry.name)) {
        throw std::invalid_argument("Reaction with the same name already exists: " + entry.name);
    }
    auto new_index = reactions_.size();
    reactions_.push_back(ReactionItem{entry, new_index});
    name_to_index_[entry.name] = new_index;
    return new_index;
}

} // namespace prexsyn::chemspace
