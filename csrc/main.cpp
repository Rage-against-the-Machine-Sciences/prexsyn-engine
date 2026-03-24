#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "chemspace/chemspace.hpp"
#include "datapipe/datapipe.hpp"
#include "descriptor/descriptor.hpp"

using namespace prexsyn;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <chemical_space_file>\n";
        return 1;
    }
    std::string cs_path = argv[1];
    std::ifstream ifs(cs_path, std::ios::binary);
    if (!ifs) {
        std::cerr << "Failed to open chemical space file: " << cs_path << "\n";
        return 1;
    }
    std::shared_ptr<chemspace::ChemicalSpace> cs = chemspace::ChemicalSpace::deserialize(ifs);

    auto dp =
        datapipe::DataPipeline(cs,
                               {{"ecfp4", descriptor::MorganFingerprint::ecfp4()},
                                {"fcfp4", descriptor::MorganFingerprint::fcfp4()}},
                               {{"synthesis", descriptor::SynthesisPostfixNotation::create(16)}});

    dp.start_workers({1, 2, 3, 4});

    for (int i = 0; i < 1000; ++i) {
        auto start_time = std::chrono::high_resolution_clock::now();

        auto batch = dp.get(1024);

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start_time);
        std::cout << "Batch " << i << " retrieved (" << duration.count() << " ms)" << "\n";
    }

    dp.stop_workers();
    return 0;
}
