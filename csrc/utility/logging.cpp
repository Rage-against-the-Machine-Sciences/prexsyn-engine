#include "logging.hpp"

#include <memory>
#include <string>

#include <spdlog/logger.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace prexsyn {

std::shared_ptr<spdlog::logger> logger() {
    static auto logger = spdlog::stderr_color_mt("prexsyn");
    return logger;
}

std::shared_ptr<spdlog::logger> create_logger(const std::string &module) {
    return spdlog::stderr_color_mt("prexsyn." + module);
}

} // namespace prexsyn
