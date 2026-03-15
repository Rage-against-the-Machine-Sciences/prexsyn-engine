#include <map>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "chemistry.hpp"

namespace {

using prexsyn::Molecule;
using prexsyn::Reaction;
using prexsyn::ReactionError;

constexpr const char *kReactionSmarts = "[NH2:1][c:2][c:3][C:4](=[O:5])[O][C].[C:6][NH2:7]>>[NH:1]"
                                        "1[c:2][c:3][C:4](=[O:5])[N:7]([C:6])C1=O";
constexpr const char *kExpectedProductSmiles = "CN(C)CCCn1c(=O)[nH]c2csnc2c1=O";

std::unique_ptr<Reaction> make_test_reaction() {
    return Reaction::from_smarts(kReactionSmarts, {"A", "B"});
}

std::shared_ptr<Molecule> make_reactant_a() { return Molecule::from_smiles("COC(=O)c1nscc1N"); }

std::shared_ptr<Molecule> make_reactant_b() { return Molecule::from_smiles("CN(C)CCCN"); }

} // namespace

TEST(ReactionTest, ApplyNamedReactantsProducesExpectedOutcome) {
    auto reaction = make_test_reaction();

    std::map<std::string, std::shared_ptr<Molecule>> reactants;
    reactants["A"] = make_reactant_a();
    reactants["B"] = make_reactant_b();

    const auto outcomes = reaction->apply(reactants);

    ASSERT_EQ(outcomes.size(), 1);
    EXPECT_FALSE(outcomes.front().empty());
    EXPECT_EQ(outcomes.front().num_products(), 1);
    EXPECT_EQ(outcomes.front().main_product()->smiles(), kExpectedProductSmiles);
}

TEST(ReactionTest, ApplyVectorReactantsTracksAssignmentAndProduct) {
    auto reaction = make_test_reaction();

    const auto outcomes = reaction->apply(std::vector{make_reactant_a(), make_reactant_b()});

    ASSERT_EQ(outcomes.size(), 1);
    ASSERT_EQ(outcomes.front().reactant_names.size(), 2);
    EXPECT_EQ(outcomes.front().reactant_names.at(0), "A");
    EXPECT_EQ(outcomes.front().reactant_names.at(1), "B");
    EXPECT_EQ(outcomes.front().main_product()->smiles(), kExpectedProductSmiles);
}

TEST(ReactionTest, ApplyNamedReactantsThrowsOnMismatchedReactantNames) {
    auto reaction = make_test_reaction();

    std::map<std::string, std::shared_ptr<Molecule>> reactants;
    reactants["A"] = make_reactant_a();
    reactants["C"] = make_reactant_b();

    EXPECT_THROW(
        {
            const auto outcomes = reaction->apply(reactants);
            (void)outcomes;
        },
        ReactionError);
}
