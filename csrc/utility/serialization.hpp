#pragma once

#include <istream>
#include <ostream>

// IWYU pragma: begin_exports
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
// IWYU pragma: end_exports

namespace prexsyn {

struct SerializationVersionTag {
    int version = 0;

    static SerializationVersionTag read(std::istream &is) {
        boost::archive::binary_iarchive ia(is);
        SerializationVersionTag tag;
        ia >> tag.version;
        return tag;
    }

    void write(std::ostream &os) const {
        boost::archive::binary_oarchive oa(os);
        oa << version;
    }

    operator int() const { return version; }
    bool operator==(const SerializationVersionTag &other) const { return version == other.version; }
    bool operator<(const SerializationVersionTag &other) const { return version < other.version; }
};

} // namespace prexsyn
