cmake_minimum_required(VERSION 3.31)
project(WirelessScheduling)

option(USE_VNS_MATH_SOLVER "VNS using Math Solvers in Local Search" OFF)
option(USE_VNS_PURE "VNS with DP in Local Search" OFF)
option(USE_VRBSP_IP "VRBSP integer program" OFF)
option(USE_MDVRBSP_IP "MDVRBSP integer program" OFF)
option(USE_MVRBSP_COL_GEN "MDVRBSP with column geneneration" OFF)
option(USE_DECODER "VRBSP with brkga" OFF)

if(NOT USE_VNS_MATH_SOLVER
   AND NOT USE_VNS_PURE
   AND NOT USE_VRBSP_IP
   AND NOT USE_MDVRBSP_IP
   AND NOT USE_MRVRBSP_COL_GEN
   AND NOT USE_DECODER)
  message(FATAL_ERROR "You need to select an option!")
endif()

option(CXX "enable C++ compilation" ON)

if(CXX)
  enable_language(CXX)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_STANDARD 20)
add_library(basic STATIC "./src/basic.cc")

include_directories("./include")
# target_include_directories(basic PUBLIC "./include")

if(USE_VNS_MATH_SOLVER
   OR USE_VRBSP_IP
   OR USE_MDVRBSP_IP
   OR USE_MDVRBSP_COL_GEN)
  target_compile_definitions(basic PUBLIC USE_MATH_SOLVER)

  find_package(GUROBI REQUIRED)
  target_include_directories(basic PUBLIC "${GUROBI_INCLUDE_DIRS}")
  target_link_libraries(basic ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
endif()

if(USE_VNS_MATH_SOLVER)
  target_compile_definitions(basic PUBLIC USE_VNS_MATH_SOLVER)
  find_package(GUROBI REQUIRED)
  add_executable(vns_math_solver "./src/vns.cc")

  target_link_libraries(vns_math_solver basic)
elseif(USE_VNS_PURE)
  target_compile_definitions(basic PUBLIC USE_VNS_PURE)
  add_executable(vns_pure "./src/vns.cc")
  # target_include_directories(vns_pure PUBLIC "./include")
  target_link_libraries(vns_pure basic)
elseif(USE_VRBSP_IP)
  find_package(GUROBI REQUIRED)
  add_executable(vrbsp_ip "./src/vrbsp.cc")
  # target_include_directories(vrbsp_ip PUBLIC "./include")
  target_link_libraries(vrbsp_ip basic)
elseif(USE_MDVRBSP_IP)
  target_compile_definitions(basic PUBLIC USE_MDVRBSP_IP)
  add_executable(mdvrbsp_ip "./src/mdvrbsp.cc")
  # target_include_directories(mdvrbsp_ip PUBLIC "./include")
  target_link_libraries(mdvrbsp_ip basic)
elseif(USE_MDVRBSP_COL_GEN)

elseif(USE_DECODER)

endif()
