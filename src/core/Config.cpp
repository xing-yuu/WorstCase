#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <wcso/core/Config.h>

namespace wcso {
// clang-format off
#define WCSO_FIELDS(X) \
    X(lame_lambda) X(lame_mu) X(gamma) X(nullspace_multiplier) \
    X(eigen_tolerance) X(linear_tolerance) X(eigen_max_iterations) \
    X(linear_max_iterations) X(gcr_restart) X(linear_solver) \
    X(w_ss) X(w_si) X(w_ef) X(w_hd) X(w_lds) X(d0) X(d1) X(d2) \
    X(epsilon_d) X(quality_threshold) X(tet_edge_length) X(eta_initial) X(eta_stop) \
    X(outer_iterations) X(inner_iterations) X(max_phase1_steps) X(max_failed_steps) \
    X(max_backtracks) X(armijo) X(step_fraction) X(step_tolerance) X(sampling_density) \
    X(threads) X(save_every) X(checkpoint_every) X(seed) X(remesh) X(debug)
// clang-format on

namespace {
template <typename T> void read(const nlohmann::json &j, const char *key, T &value) {
    if (!j.contains(key))
        return;
    const auto &v = j.at(key);
    bool valid = false;
    if constexpr (std::is_same_v<T, bool>)
        valid = v.is_boolean();
    else if constexpr (std::is_same_v<T, std::string>)
        valid = v.is_string();
    else if constexpr (std::is_integral_v<T>)
        valid = v.is_number_integer();
    else
        valid = v.is_number();
    if (!valid)
        throw std::invalid_argument(std::string("Wrong type for config key: ") + key);
    if constexpr (std::is_integral_v<T> && !std::is_same_v<T, bool>) {
        const long double number = v.get<long double>();
        if (number < static_cast<long double>(std::numeric_limits<T>::lowest()) ||
            number > static_cast<long double>(std::numeric_limits<T>::max()))
            throw std::invalid_argument(std::string("Integer out of range: ") + key);
    }
    value = v.get<T>();
}
} // namespace
nlohmann::json Config::to_json() const {
    nlohmann::json j = nlohmann::json::object();
#define PUT(name) j[#name] = name;
    WCSO_FIELDS(PUT)
#undef PUT
    return j;
}
Config Config::from_json(const nlohmann::json &j) {
    if (!j.is_object())
        throw std::invalid_argument("Configuration must be a JSON object");
    Config config;
    const auto known = config.to_json();
    for (const auto &item : j.items())
        if (!known.contains(item.key()))
            throw std::invalid_argument("Unknown config key: " + item.key());
#define GET(name) read(j, #name, config.name);
    WCSO_FIELDS(GET)
#undef GET
    config.validate();
    return config;
}
void Config::validate() const {
    const auto serialized = to_json();
    for (const auto &item : serialized.items()) {
        if (item.value().is_number_float() && !std::isfinite(item.value().get<double>()))
            throw std::invalid_argument("Non-finite config value: " + item.key());
    }
    const auto require = [](bool condition, const char *message) {
        if (!condition)
            throw std::invalid_argument(message);
    };
    require(lame_mu > 0 && 3 * lame_lambda + 2 * lame_mu > 0,
            "Lame parameters require mu>0 and 3*lambda+2*mu>0");
    require(gamma > 0 && nullspace_multiplier > 0, "Stress parameters must be positive");
    require(eigen_tolerance > 0 && eigen_tolerance < 1 && linear_tolerance > 0 && linear_tolerance < 1,
            "Solver tolerances must be in (0,1)");
    require(linear_solver == "cg" || linear_solver == "gcr", "linear_solver must be cg or gcr");
    require(eigen_max_iterations > 0 && linear_max_iterations > 0 && gcr_restart > 0,
            "Solver iteration limits must be positive");
    require(w_ss > 0 && w_si > 0 && w_ef > 0 && w_hd >= 0 && w_lds >= 0,
            "Constraint weights must be positive; shape weights must be nonnegative");
    require(d0 > 0 && d1 > 0 && d2 > 0 && epsilon_d > 0 && tet_edge_length > 0,
            "Geometric thresholds must be positive");
    require(quality_threshold >= 0 && quality_threshold < 2.01, "Invalid quality threshold");
    require(eta_initial > 0 && eta_initial < eta_stop && eta_stop < 1, "Require 0<eta_initial<eta_stop<1");
    require(outer_iterations > 0 && inner_iterations > 0 && max_phase1_steps > 0 && max_failed_steps > 0 &&
                max_backtracks > 0,
            "Iteration limits must be positive");
    require(armijo > 0 && armijo < 1 && step_fraction > 0 && step_fraction <= .5 && step_tolerance > 0,
            "Invalid line search parameters");
    require(sampling_density > 0 && sampling_density <= 1e8, "Invalid sampling_density");
    require(threads > 0 && threads <= 256 && save_every >= 0 && checkpoint_every >= 0,
            "Invalid thread count or output interval");
}
#undef WCSO_FIELDS
} // namespace wcso
