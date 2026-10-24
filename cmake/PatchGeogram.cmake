# Geogram 1.9.6 bypasses Numeric::random_reset() in its shuffle helper.
# Preserve the upstream BSD header and make the fixed dependency reproducible.
FetchContent_GetProperties(geogram SOURCE_DIR _wcso_geogram_source_dir)
set(_wcso_shuffle_header "${_wcso_geogram_source_dir}/src/lib/geogram/basic/algorithm.h")
file(READ "${_wcso_shuffle_header}" _wcso_shuffle_source)
if(NOT _wcso_shuffle_source MATCHES "WCSO deterministic mesher shuffle")
  string(REGEX REPLACE
    "std::random_device rng;[\r\n \t]*std::mt19937 urng\\(rng\\(\\)\\);"
    "// WCSO deterministic mesher shuffle; upstream BSD license retained.\n        std::mt19937 urng(42);"
    _wcso_shuffle_patched "${_wcso_shuffle_source}")
  if(_wcso_shuffle_patched STREQUAL _wcso_shuffle_source)
    message(FATAL_ERROR "Pinned Geogram shuffle patch no longer matches its source")
  endif()
  file(WRITE "${_wcso_shuffle_header}" "${_wcso_shuffle_patched}")
endif()
