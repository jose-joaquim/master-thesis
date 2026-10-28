cmake_minimum_required(VERSION 3.31)
set(CMAKE_CXX_STANDARD 20)

project(WirelessScheduling)

option(USE_VNS_MATH_SOLVER "VNS using Math Solvers in Local Search" OFF)
option(USE_VNS_PURE "VNS with DP in Local Search" OFF)
option(USE_VRBSP_IP "VRBSP integer program" OFF)
option(USE_MDVRBSP_IP "MDVRBSP integer program" OFF)
option(USE_MVRBSP_COL_GEN "MDVRBSP with column geneneration" OFF)
option(USE_BRKGA "MDVRBSP with column geneneration" OFF)

if(NOT USE_VNS_MATH_SOLVER
   AND NOT USE_VNS_PURE
   AND NOT USE_VRBSP_IP
   AND NOT USE_MDVRBSP_IP
   AND NOT USE_MDVRBSP_COL_GEN
   AND NOT USE_VRBSP_COL_GEN
  AND NOT USE_BRKGA)
  message(FATAL_ERROR "You need to select an option!")
endif()

option(CXX "enable C++ compilation" ON)

if(CXX)
  enable_language(CXX)
endif()

include(FetchContent)

FetchContent_Declare(absl
                     GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git)

FetchContent_Declare(
  json
  URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz)

FetchContent_MakeAvailable(json absl)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
add_library(basic STATIC "./src/basic.cc")

if(USE_BRKGA)
  if(APPLE)
    execute_process(
      COMMAND brew --prefix libomp
      OUTPUT_VARIABLE HOMEBREW_OMP_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(HOMEBREW_OMP_PREFIX)
      set(OpenMP_ROOT "${HOMEBREW_OMP_PREFIX}")
    endif()
  endif()

  add_library(mtrand STATIC "./src/MTRand.cc")

  find_package(OpenMP REQUIRED)
  add_library(brkga STATIC "./src/brkga.cc" "./src/population.cc")
  target_link_libraries(
    brkga PUBLIC OpenMP::OpenMP_CXX absl::log absl::log_initialize
                 nlohmann_json::nlohmann_json mtrand)
endif()

include_directories("./include")

if(USE_VNS_MATH_SOLVER
   OR USE_VRBSP_IP
   OR USE_MDVRBSP_IP
   OR USE_MDVRBSP_COL_GEN
   OR USE_VRBSP_COL_GEN)
  target_compile_definitions(basic PUBLIC USE_VRBSP_IP)

  find_package(GUROBI REQUIRED)
  target_include_directories(basic PUBLIC "${GUROBI_INCLUDE_DIRS}")
  target_link_libraries(basic ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
  target_link_libraries(basic nlohmann_json::nlohmann_json absl::log
                        absl::log_initialize absl::flags absl::flags_parse)
else()
  target_link_libraries(
    basic
    nlohmann_json::nlohmann_json
    absl::log
    absl::flags
    absl::flags_parse
    absl::log_globals
    absl::log_initialize)
endif()

if(USE_VNS_MATH_SOLVER)
  target_compile_definitions(basic PUBLIC USE_VNS_MATH_SOLVER)
  add_executable(vns_math_solver "./src/vns.cc")
  target_link_libraries(vns_math_solver basic)
elseif(USE_VNS_PURE)
  target_compile_definitions(basic PUBLIC USE_VNS_PURE)
  add_executable(vns_pure "./src/vns.cc")
  target_link_libraries(vns_pure basic nlohmann_json::nlohmann_json)
elseif(USE_VRBSP_IP)
  target_compile_definitions(basic PUBLIC USE_VRBSP_IP)
  add_executable(vrbsp_ip "./src/vrbsp.cc")
  target_link_libraries(vrbsp_ip basic)
elseif(USE_MDVRBSP_IP)
  target_compile_definitions(basic PUBLIC USE_MDVRBSP_IP)
  add_executable(mdvrbsp_ip "./src/mdvrbsp.cc")
  target_link_libraries(mdvrbsp_ip basic)
elseif(USE_VRBSP_COL_GEN)
  target_compile_definitions(basic PUBLIC USE_VRBSP_IP)
  add_executable(vrbsp_col_gen "./src/vrbsp_col_gen.cc")
  target_link_libraries(vrbsp_col_gen basic)
elseif(USE_MDVRBSP_COL_GEN)
  target_compile_definitions(basic PUBLIC USE_MDVRBSP_COL_GEN)
  add_executable(mdvrbsp_col_gen "./src/mdvrbsp_col_gen.cc")
  target_link_libraries(mdvrbsp_col_gen basic nlohmann_json::nlohmann_json)
elseif(USE_BRKGA)
  target_compile_definitions(basic PUBLIC USE_BRKGA)

  add_executable(brkga_decoder "./src/brkga_decoder.cc")
  target_link_libraries(brkga_decoder brkga basic mtrand
                        nlohmann_json::nlohmann_json ${OMP_LIBRARY})
elseif(USE_BRKGA)
  target_compile_definitions(basic PUBLIC USE_BRKGA)

  add_executable(brkga_decoder "./src/brkga_decoder.cc")
  target_link_libraries(brkga_decoder brkga basic mtrand
                        nlohmann_json::nlohmann_json ${OMP_LIBRARY})
endif()
