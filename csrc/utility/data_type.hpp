#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace prexsyn {

struct DataType {
    enum T : std::uint8_t {
        float32,
        int64,
    };

    template <typename U> static constexpr T get_dtype() {
        if constexpr (std::is_same_v<U, float>) {
            return T::float32;
        } else if constexpr (std::is_same_v<U, std::int64_t>) {
            return T::int64;
        } else {
            static_assert(!std::is_same_v<U, U>, "Unsupported data type");
        }
    }

    static size_t get_size(const T &t) {
        switch (t) {
        case T::float32:
            return 4;
        case T::int64:
            return 8;
        }
        return 0;
    }
};

} // namespace prexsyn
