#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "../chemistry/chemistry.hpp"
#include "bb_lib.hpp"
#include "rxn_lib.hpp"

namespace prexsyn::chemspace {

class PostfixNotation {
public:
    struct Token {
        size_t index;
        enum Type { BuildingBlock = 1, Reaction = 2 } type; // NOLINT(performance-enum-size)

        template <typename Archive> void serialize(Archive &ar, const unsigned int /* version */) {
            ar & index;
            ar & type;
        }
    };

private:
    std::vector<Token> tokens_;

public:
    template <typename Archive> void serialize(Archive &ar, const unsigned int /* version */) {
        ar & tokens_;
    }

    const std::vector<Token> &tokens() const { return tokens_; }

    void append(size_t index, Token::Type type) {
        tokens_.push_back(Token{.index = index, .type = type});
    }

    void extend(const std::vector<Token> &new_tokens) {
        tokens_.insert(tokens_.end(), new_tokens.begin(), new_tokens.end());
    }
};

class ChemicalSpace;

class ChemicalSpaceSynthesis {
private:
    const ChemicalSpace &cs_;
    PostfixNotation postfix_notation_;
    std::shared_ptr<Synthesis> synthesis_;

public:
    ChemicalSpaceSynthesis(const ChemicalSpace &cs)
        : cs_(cs), postfix_notation_(), synthesis_(std::make_shared<Synthesis>()) {}

    const PostfixNotation &postfix_notation() const { return postfix_notation_; }
    const Synthesis &synthesis() const { return *synthesis_; }
    std::shared_ptr<Synthesis> &get_synthesis() { return synthesis_; }

    struct Result {
        bool is_ok;
        std::string message;

        static Result ok() { return {.is_ok = true, .message = {}}; }
        static Result error(const std::string &msg) { return {.is_ok = false, .message = msg}; }

        operator bool() const { return is_ok; }
    };

    Result add_building_block(BuildingBlockLibrary::Index) noexcept;
    Result add_building_block(const std::string &) noexcept;
    Result add_reaction(ReactionLibrary::Index) noexcept;
    Result add_reaction(const std::string &) noexcept;
};

} // namespace prexsyn::chemspace
