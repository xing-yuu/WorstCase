include(FetchContent)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

FetchContent_Declare(eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0 GIT_SHALLOW TRUE)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
  FetchContent_Populate(eigen)
endif()
add_library(wcso_eigen INTERFACE)
add_library(Eigen3::Eigen ALIAS wcso_eigen)
target_include_directories(wcso_eigen SYSTEM INTERFACE "${eigen_SOURCE_DIR}")

FetchContent_Declare(spectra
  URL https://codeload.github.com/yixuan/spectra/tar.gz/fef8a396048040a411ed0d880ac7b82d63fb98e3)
FetchContent_GetProperties(spectra)
if(NOT spectra_POPULATED)
  FetchContent_Populate(spectra)
endif()
add_library(wcso_spectra INTERFACE)
target_include_directories(wcso_spectra SYSTEM INTERFACE "${spectra_SOURCE_DIR}/include")

FetchContent_Declare(amgcl
  URL https://codeload.github.com/ddemidov/amgcl/tar.gz/73a0b8dd8494ff5d4253685b20effcc736359c5f)
FetchContent_GetProperties(amgcl)
if(NOT amgcl_POPULATED)
  FetchContent_Populate(amgcl)
endif()
add_library(wcso_amgcl INTERFACE)
target_include_directories(wcso_amgcl SYSTEM INTERFACE "${amgcl_SOURCE_DIR}")
set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
FetchContent_Declare(json
  URL https://codeload.github.com/nlohmann/json/tar.gz/9cca280a4d0ccf0c08f47a99aa71d1b0e52f8d03)
FetchContent_MakeAvailable(json)
# Reuse this target when the optional mesher requests its older JSON dependency.
add_library(json INTERFACE)
target_link_libraries(json INTERFACE nlohmann_json::nlohmann_json)

set(IPC_TOOLKIT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(IPC_TOOLKIT_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
set(IPC_TOOLKIT_WITH_ABSEIL OFF CACHE BOOL "" FORCE)
set(IPC_TOOLKIT_WITH_ROBIN_MAP OFF CACHE BOOL "" FORCE)
set(IPC_TOOLKIT_WITH_FILIB OFF CACHE BOOL "" FORCE)
set(IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION OFF CACHE BOOL "" FORCE)
set(IPC_TOOLKIT_WITH_CUDA OFF CACHE BOOL "" FORCE)
set(LIBIGL_PREDICATES ON CACHE BOOL "" FORCE)
if(WCSO_WITH_TETWILD)
  # Share the same fmt ABI between IPC's spdlog and fTetWild.
  set(FMT_TEST OFF CACHE BOOL "" FORCE)
  set(FMT_INSTALL ON CACHE BOOL "" FORCE)
  FetchContent_Declare(fmt URL https://codeload.github.com/fmtlib/fmt/tar.gz/9.1.0)
  FetchContent_MakeAvailable(fmt)
  set(SPDLOG_FMT_EXTERNAL ON CACHE BOOL "" FORCE)
endif()
FetchContent_Declare(ipc_toolkit
  URL https://codeload.github.com/ipc-sim/ipc-toolkit/tar.gz/1885ce847e10498d239fbaed341526e62c4366c1)
FetchContent_MakeAvailable(ipc_toolkit)

set(BUILD MMG3D CACHE STRING "" FORCE)
set(LIBMMG3D_STATIC OFF CACHE BOOL "" FORCE)
set(LIBMMG3D_SHARED ON CACHE BOOL "" FORCE)
set(USE_SCOTCH OFF CACHE BOOL "" FORCE)
set(USE_ELAS OFF CACHE BOOL "" FORCE)
FetchContent_Declare(mmg
  URL https://codeload.github.com/MmgTools/mmg/tar.gz/588182d6497961435bb08221832f1358a957da64)
set(_wcso_previous_shared "${BUILD_SHARED_LIBS}")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries" FORCE)
FetchContent_MakeAvailable(mmg)
set(BUILD_SHARED_LIBS "${_wcso_previous_shared}" CACHE BOOL "Build shared libraries" FORCE)
if(NOT TARGET libmmg3d_so)
  message(FATAL_ERROR "MMG3D shared-library target was not created")
endif()
add_library(wcso_mmg INTERFACE)
target_link_libraries(wcso_mmg INTERFACE libmmg3d_so)
target_include_directories(wcso_mmg SYSTEM INTERFACE "${mmg_BINARY_DIR}/include"
  "${mmg_SOURCE_DIR}/src/mmg3d" "${mmg_SOURCE_DIR}/src/common")

if(WCSO_WITH_TETWILD)
  # fTetWild itself requires GMP. Pass GMP_INCLUDE_DIRS/GMP_LIBRARIES or
  # use a package-manager toolchain; the optimizer-only build needs neither.
  set(FLOAT_TETWILD_ENABLE_TBB ON CACHE BOOL "" FORCE)
  set(FLOAT_TETWILD_WITH_EXACT_ENVELOPE OFF CACHE BOOL "" FORCE)
  set(GEOGRAM_WITH_TETGEN OFF CACHE BOOL "" FORCE)
  set(GEOGRAM_WITH_HLBFGS OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(ftetwild
    URL https://codeload.github.com/wildmeshing/fTetWild/tar.gz/8118f810478e0e65a7bf2d8cecdc5e203a876e97)
  FetchContent_MakeAvailable(ftetwild)
  include("${CMAKE_CURRENT_LIST_DIR}/PatchGeogram.cmake")
endif()
