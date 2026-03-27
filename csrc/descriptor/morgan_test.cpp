#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "../chemistry/chemistry.hpp"
#include "morgan.hpp"

namespace {

using prexsyn::Molecule;
using prexsyn::descriptor::MorganFingerprint;

struct FingerprintReference {
    std::string smiles;
    std::vector<size_t> ecfp4_nonzero_indices;
    std::vector<size_t> fcfp4_nonzero_indices;
};

std::filesystem::path reference_csv_path() {
    return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() /
           "resources/test/fp_ref.csv";
}

std::vector<size_t> parse_indices(const std::string &field) {
    std::vector<size_t> indices;
    if (field.empty()) {
        return indices;
    }

    std::istringstream stream(field);
    std::string token;
    while (std::getline(stream, token, ';')) {
        if (!token.empty()) {
            indices.push_back(static_cast<size_t>(std::stoul(token)));
        }
    }
    return indices;
}

FingerprintReference parse_reference_line(const std::string &line) {
    const size_t first_comma = line.find(',');
    const size_t second_comma = line.find(',', first_comma + 1);
    const size_t third_comma = line.find(',', second_comma + 1);

    if (first_comma == std::string::npos || second_comma == std::string::npos ||
        third_comma == std::string::npos) {
        throw std::runtime_error("Malformed fingerprint reference CSV row");
    }

    FingerprintReference row;
    row.smiles = line.substr(0, first_comma);
    row.ecfp4_nonzero_indices =
        parse_indices(line.substr(first_comma + 1, second_comma - first_comma - 1));
    row.fcfp4_nonzero_indices =
        parse_indices(line.substr(second_comma + 1, third_comma - second_comma - 1));
    return row;
}

std::vector<FingerprintReference> load_reference_rows() {
    const auto csv_path = reference_csv_path();
    std::ifstream input(csv_path);
    if (!input.is_open()) {
        throw std::runtime_error("Cannot open fingerprint reference file: " + csv_path.string());
    }

    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("Fingerprint reference CSV is empty");
    }

    std::vector<FingerprintReference> rows;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (!line.empty()) {
            rows.push_back(parse_reference_line(line));
        }
    }
    return rows;
}

std::vector<size_t> compute_nonzero_indices(const MorganFingerprint &descriptor,
                                            const Molecule &molecule) {
    std::vector<std::byte> output(descriptor.size_in_bytes());
    std::span<std::byte> out_span(output.data(), output.size());
    descriptor(molecule, out_span);

    std::vector<size_t> indices;
    indices.reserve(output.size());
    for (size_t i = 0; i < output.size(); ++i) {
        if (std::to_integer<unsigned char>(output[i]) != 0U) {
            indices.push_back(i);
        }
    }
    return indices;
}

} // namespace

TEST(MorganFingerprintTest, ECFP4MatchesReferenceBitsFromCsv) {
    const auto references = load_reference_rows();
    ASSERT_FALSE(references.empty());

    auto descriptor = MorganFingerprint::ecfp4();
    ASSERT_TRUE(descriptor);

    for (const auto &row : references) {
        SCOPED_TRACE("SMILES: " + row.smiles);

        const auto molecule = Molecule::from_smiles(row.smiles);
        ASSERT_TRUE(molecule);

        const auto actual = compute_nonzero_indices(*descriptor, *molecule);
        EXPECT_EQ(actual, row.ecfp4_nonzero_indices);
    }
}

TEST(MorganFingerprintTest, FCFP4MatchesReferenceBitsFromCsv) {
    const auto references = load_reference_rows();
    ASSERT_FALSE(references.empty());

    auto descriptor = MorganFingerprint::fcfp4();
    ASSERT_TRUE(descriptor);

    for (const auto &row : references) {
        SCOPED_TRACE("SMILES: " + row.smiles);

        const auto molecule = Molecule::from_smiles(row.smiles);
        ASSERT_TRUE(molecule);

        const auto actual = compute_nonzero_indices(*descriptor, *molecule);
        EXPECT_EQ(actual, row.fcfp4_nonzero_indices);
    }
}
