#include <pybind11/pybind11.h>

#include "chemistry/binding.hpp"
#include "chemspace/binding.hpp"

namespace py = pybind11;

PYBIND11_MODULE(prexsyn_engine, m) {
    m.doc() = "PrexSyn Engine";

    auto m_chemistry = m.def_submodule("chemistry");
    def_module_chemistry(m_chemistry);

    auto m_chemspace = m.def_submodule("chemspace");
    def_module_chemspace(m_chemspace);
}
