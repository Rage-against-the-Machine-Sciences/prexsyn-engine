#pragma once

#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace prexsyn::chemspace {

class PostfixNotation {
public:
    struct Token {
        enum Type { BuildingBlock = 1, Reaction = 2 }; // NOLINT(performance-enum-size)
        size_t index;
        Type type;

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

    size_t size() const { return tokens_.size(); }
    const std::vector<Token> &tokens() const { return tokens_; }

    void append(size_t index, Token::Type type) {
        tokens_.push_back(Token{.index = index, .type = type});
    }

    void extend(const std::vector<Token> &new_tokens) {
        tokens_.insert(tokens_.end(), new_tokens.begin(), new_tokens.end());
    }

    void pop_back() {
        if (tokens_.empty()) {
            throw std::out_of_range("Cannot pop_back from an empty PostfixNotation");
        }
        tokens_.pop_back();
    }
};

inline std::ostream &operator<<(std::ostream &os, const PostfixNotation &pn) {
    os << "PostfixNotation(";
    for (const auto &token : pn.tokens()) {
        os << (token.type == PostfixNotation::Token::Type::BuildingBlock ? "b" : "r") << token.index
           << " ";
    }
    os << ")";
    return os;
}

} // namespace prexsyn::chemspace
