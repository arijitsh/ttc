cmake_minimum_required(VERSION 3.10)
project(ttc)

include(CTest)

# Option to build a static binary (default is OFF)
option(BUILD_STATIC "Build static binary" OFF)

# Optionally build with debugging symbols
option(ENABLE_DEBUG "Enable debug symbols" OFF)

# Optional DDNNF/D4 counting support (default OFF)
option(ENABLE_DDNNF "Enable DDNNF/D4 counting backends" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(ENABLE_DEBUG)
    message(STATUS "Debugging symbols enabled")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
    add_compile_options(-g)
else()
    # Default to an optimized build.  The numerical core (Eigen + volesti) is
    # orders of magnitude slower without optimization, so always enable -O3.
    if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
    endif()
    add_compile_options(-O3 -DNDEBUG)
endif()

if(BUILD_STATIC)
    message(STATUS "Building a static binary")
    # Linux-specific flag for static linking; adjust for your environment.
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif()

# Enable feature macro when DDNNF/D4 support is requested
if(ENABLE_DDNNF)
    add_compile_definitions(TTC_ENABLE_DDNNF)
endif()

# Locate cvc5 headers and libraries manually
find_path(CVC5_INCLUDE_DIR NAMES cvc5/cvc5.h
          HINTS $ENV{HOME}/solvers/cvc5_solver/include)
find_library(CVC5_LIB NAMES cvc5
            HINTS $ENV{HOME}/solvers/cvc5_solver/build/src)
find_library(CVC5PARSER_LIB NAMES cvc5parser
            HINTS $ENV{HOME}/solvers/cvc5_solver/build/src/parser)
if(NOT CVC5_INCLUDE_DIR OR NOT CVC5_LIB OR NOT CVC5PARSER_LIB)
  message(FATAL_ERROR "cvc5 library not found. Please install libcvc5-dev and libcvc5parser1")
endif()
include_directories(${CVC5_INCLUDE_DIR})
# Print the full paths of cvc5 headers and libraries
message(STATUS "CVC5 include directory: ${CVC5_INCLUDE_DIR}")
message(STATUS "CVC5 library: ${CVC5_LIB}")
message(STATUS "CVC5Parser library: ${CVC5PARSER_LIB}")
# Find ISL via pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(ISL REQUIRED IMPORTED_TARGET isl)
# Find Boost (we need program_options)
find_package(Boost REQUIRED COMPONENTS program_options)
if(Boost_FOUND)
    message(STATUS "Found Boost: ${Boost_INCLUDE_DIRS}")
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Boost not found!")
endif()

# cvc5 static builds depend on the poly, CaDiCaL and CryptoMiniSat libraries.
# Prefer the versions shipped with the local cvc5 build (e.g. the XOR-patched
# CaDiCaL), falling back to system-wide installs.
set(CVC5_DEPS_LIB_DIR "$ENV{HOME}/solvers/cvc5_solver/build/deps/lib")
find_library(PICPOLY_LIB picpoly HINTS ${CVC5_DEPS_LIB_DIR})
if(NOT PICPOLY_LIB)
  message(WARNING "picpoly library not found; continuing without explicit linkage")
  set(PICPOLY_LIB "")
endif()
find_library(PICPOLYXX_LIB picpolyxx HINTS ${CVC5_DEPS_LIB_DIR})
if(NOT PICPOLYXX_LIB)
  message(WARNING "picpolyxx library not found; continuing without explicit linkage")
  set(PICPOLYXX_LIB "")
endif()
find_library(CADICAL_LIB cadical HINTS ${CVC5_DEPS_LIB_DIR})
if(NOT CADICAL_LIB)
  message(WARNING "CaDiCaL library not found; continuing without explicit linkage")
  set(CADICAL_LIB "")
endif()
# CryptoMiniSat is only required when cvc5 was built statically with CMS support.
find_library(CMS_LIB cryptominisat5 HINTS ${CVC5_DEPS_LIB_DIR})
if(NOT CMS_LIB)
  message(STATUS "CryptoMiniSat library not found; continuing without explicit linkage")
  set(CMS_LIB "")
endif()

# cddlib (double precision) for polytope canonicalization in the LRA volume path.
find_path(CDD_INCLUDE_DIR cddlib/cdd.h)
find_library(CDD_LIB NAMES cdd)
if(NOT CDD_INCLUDE_DIR OR NOT CDD_LIB)
  message(FATAL_ERROR "cddlib not found. Install libcdd-dev (provides cddlib/cdd.h and libcdd).")
endif()
message(STATUS "cddlib include: ${CDD_INCLUDE_DIR}")
message(STATUS "cddlib library: ${CDD_LIB}")
include_directories(${CDD_INCLUDE_DIR})
find_library(COLAMD_LIB colamd)
if(NOT COLAMD_LIB)
  message(WARNING "COLAMD library not found; continuing without explicit linkage")
  set(COLAMD_LIB "")
endif()
find_library(AMD_LIB amd)
if(NOT AMD_LIB)
  message(WARNING "AMD library not found; continuing without explicit linkage")
  set(AMD_LIB "")
endif()

# Add our project source directories
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/ddnnf)
include_directories(${PROJECT_SOURCE_DIR}/src/csb)
include_directories(${PROJECT_SOURCE_DIR}/src/pact)
include_directories(${PROJECT_SOURCE_DIR}/src/tocnf)
include_directories(${PROJECT_SOURCE_DIR}/src/volesti/include)
include_directories(${PROJECT_SOURCE_DIR}/src/volesti/external)
include_directories(${PROJECT_SOURCE_DIR}/src/allsat-circuits/src)
include_directories(${PROJECT_SOURCE_DIR}/src/allsat-circuits)

find_library(LPSOLVE_LIB lpsolve55)
find_path(LPSOLVE_INCLUDE_DIR lp_lib.h PATH_SUFFIXES lpsolve)
if(LPSOLVE_LIB AND LPSOLVE_INCLUDE_DIR)
    message(STATUS "Found lpsolve library: ${LPSOLVE_LIB}")
    include_directories(${LPSOLVE_INCLUDE_DIR})
else()
    if(NOT LPSOLVE_LIB)
        message(WARNING "lpsolve library not found; disabling lpsolve integration")
    elseif(NOT LPSOLVE_INCLUDE_DIR)
        message(WARNING "lpsolve headers not found; disabling lpsolve integration")
    endif()
    add_compile_definitions(DISABLE_LPSOLVE)
    set(LPSOLVE_LIB "")
endif()

find_path(EIGEN3_INCLUDE_DIR Eigen/Dense
          HINTS $ENV{EIGEN3_INCLUDE_DIR} /usr/include/eigen3)
if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Eigen3 headers not found. Please install libeigen3-dev")
endif()
include_directories(${EIGEN3_INCLUDE_DIR})

set(SOURCES
    src/main.cpp
    src/parser.cpp
    src/pact/pact.cpp
    src/pact/saturatingcount.cpp
    src/pact/hash_constraint.cpp
    src/pact/bool_hash.cpp
    src/pact/bv_hash.cpp
    src/projected_enumerator.cpp
    src/csb/stp_sat.cpp
    src/tocnf/tocnf.cpp
    src/eager/bvcnf.cpp
    src/weighted/tounweighted.cpp
    src/volume/allsat_volume.cpp
    src/volume/lra_volume.cpp
    src/allsat-circuits/src/AllSatAlgo/AllSatAlgoBase.cpp
    src/allsat-circuits/src/AllSatAlgo/Blocking/AllSatAlgoBlockingBase.cpp
    src/allsat-circuits/src/AllSatAlgo/Blocking/TseitinEnc/AllSatAlgoTseitinEnc.cpp
    src/allsat-circuits/src/AllSatSolver/AllSatSolverBase.cpp
    src/allsat-circuits/src/AllSatSolver/Topor/AllSatSolverTopor.cpp
    src/allsat-circuits/src/AllSatSolver/Ipasir/AllSatSolverIpasir.cpp
    src/allsat-circuits/src/CirSimulation/CirSim.cpp
    src/allsat-circuits/src/Globals/ipasir_stub.cpp
)

if(ENABLE_DDNNF)
    include_directories(${PROJECT_SOURCE_DIR}/src/ddnnf/flowcutter/src)
    include_directories(${PROJECT_SOURCE_DIR}/src/d4)
    include_directories(${PROJECT_SOURCE_DIR}/src/d4/src)

    list(APPEND SOURCES
        src/ddnnf/proj_ddnnf.cpp
        src/ddnnf/cache.cpp
        src/ddnnf/var_order.cpp
        src/ddnnf/iflowcutter.cpp
        src/ddnnf/flowcutter/src/greedy_order.cpp
        src/ddnnf/TreeDecomposition.cpp
        src/ddnnf/arjun.cpp
        src/d4_interface.cpp
    )

    # GMP for big integer arithmetic used by D4
    find_library(GMP_LIB gmp)
    find_library(GMPXX_LIB gmpxx)

    # D4 sources (trimmed to avoid heavy optional dependencies)
    set(D4_SOURCES
        src/d4/src/problem/ProblemManager.cpp
        src/d4/src/problem/ProblemTypes.cpp
        src/d4/src/problem/cnf/ProblemManagerCnf.cpp
        src/d4/src/problem/cnf/ParserDimacs.cpp
        src/d4/src/solvers/ActivityManager.cpp
        src/d4/src/solvers/WrapperSolver.cpp
        src/d4/src/solvers/cnf/WrapperMinisat.cpp
        src/d4/src/solvers/cnf/minisat/Solver.cpp
        src/d4/src/formulaManager/FormulaManager.cpp
        src/d4/src/formulaManager/cnf/CnfManager.cpp
        src/d4/src/formulaManager/cnf/CnfManagerDyn.cpp
        src/d4/src/formulaManager/cnf/CnfManagerDynBlockedCl.cpp
        src/d4/src/formulaManager/cnf/CnfManagerDynPure.cpp
        src/d4/src/formulaManager/cnf/ResidualSimplifier.cpp
        src/d4/src/caching/DataInfo.cpp
        src/d4/src/caching/bucket/BucketAllocator.cpp
        src/d4/src/caching/bucket/cnf/BucketInConstruction.cpp
        src/d4/src/heuristics/branchingHeuristic/BranchingHeuristic.cpp
        src/d4/src/heuristics/branchingHeuristic/BranchingHeuristicClassic.cpp
        src/d4/src/heuristics/branchingHeuristic/BranchingHeuristicLargeArity.cpp
        src/d4/src/heuristics/branchingHeuristic/BranchingHeuristicHybridPartialClassic.cpp
        src/d4/src/heuristics/partialOrder/PartialOrderHeuristic.cpp
        src/d4/src/heuristics/partialOrder/PartialOrderHeuristicGiven.cpp
        src/d4/src/heuristics/partialOrder/PartialOrderHeuristicNone.cpp
        src/d4/src/heuristics/phaseSelection/PhaseHeuristic.cpp
        src/d4/src/heuristics/phaseSelection/PhaseHeuristicFalse.cpp
        src/d4/src/heuristics/phaseSelection/PhaseHeuristicPolarity.cpp
        src/d4/src/heuristics/phaseSelection/PhaseHeuristicTrue.cpp
        src/d4/src/heuristics/phaseSelection/cnf/PhaseHeuristicOccurrence.cpp
        src/d4/src/heuristics/scoringVariable/ScoringMethod.cpp
        src/d4/src/heuristics/scoringVariable/cnf/ScoringMethodDlcs.cpp
        src/d4/src/heuristics/scoringVariable/cnf/ScoringMethodJwts.cpp
        src/d4/src/heuristics/scoringVariable/cnf/ScoringMethodMom.cpp
        src/d4/src/heuristics/scoringVariable/cnf/ScoringMethodVsads.cpp
        src/d4/src/heuristics/scoringVariable/cnf/ScoringMethodVsids.cpp
        src/d4/src/methods/MethodManager.cpp
        src/d4/src/methods/QueryManager.cpp
        src/d4/src/options/cache/OptionCacheManager.cpp
        src/d4/src/options/branchingHeuristic/OptionBranchingHeuristic.cpp
        src/d4/src/options/branchingHeuristic/OptionPartialOrderHeuristic.cpp
        src/d4/src/options/methods/OptionDpllStyleMethod.cpp
        src/d4/src/options/methods/OptionProjMcMethod.cpp
        src/d4/src/utils/AtMost1Extractor.cpp
        src/d4/src/utils/EquivExtractor.cpp
        src/d4/src/utils/Parsing.cpp
    )

    add_library(d4lib ${D4_SOURCES})

    # D4 depends on GMP for big integer support
    target_link_libraries(d4lib ${GMPXX_LIB} ${GMP_LIB})
endif()


# Attempt to find STP for optional satisfiability checking
find_path(STP_INCLUDE_DIR NAMES stp/cpp_interface.h)
find_library(STP_LIB NAMES stp)
if(STP_INCLUDE_DIR AND STP_LIB)
    message(STATUS "Found STP: ${STP_INCLUDE_DIR}")
    include_directories(${STP_INCLUDE_DIR})
    add_definitions(-DTTC_HAS_STP)
else()
    message(STATUS "STP library not found; STP support disabled")
    set(STP_LIB "")
endif()

add_executable(ttc ${SOURCES})

# Link against the cvc5 libraries, Boost, optional STP and optional D4 support
target_link_libraries(ttc
    ${CVC5PARSER_LIB}
    ${CVC5_LIB}
    ${PICPOLYXX_LIB}
    ${PICPOLY_LIB}
    ${CADICAL_LIB}
    ${CMS_LIB}
    ${CDD_LIB}
    ${Boost_LIBRARIES}
    ${STP_LIB}
    ${LPSOLVE_LIB}
    ${COLAMD_LIB}
    ${AMD_LIB}
    PkgConfig::ISL
    m)

find_library(GMP_LIB gmp REQUIRED)
find_library(GMPXX_LIB gmpxx REQUIRED)
target_link_libraries(ttc ${GMP_LIB} ${GMPXX_LIB})

if(ENABLE_DDNNF)
    add_dependencies(ttc d4lib)
    target_link_libraries(ttc d4lib ${GMPXX_LIB} ${GMP_LIB})
endif()

if(BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    set(LRA_VOLUME_TEST_SCRIPT ${CMAKE_SOURCE_DIR}/tests/test_lra_volume.py)

    add_test(
        NAME lra_volume_union
        COMMAND
            ${Python3_EXECUTABLE}
            ${LRA_VOLUME_TEST_SCRIPT}
            $<TARGET_FILE:ttc>
            ${CMAKE_SOURCE_DIR}/tests/lra_union.smt2
            --expected-polytopes 2
            --expected-rows 2
            --expected-volumes 400
            --volume-tolerance 120
            --final-range 650 820)

    add_test(
        NAME lra_volume_single_box
        COMMAND
            ${Python3_EXECUTABLE}
            ${LRA_VOLUME_TEST_SCRIPT}
            $<TARGET_FILE:ttc>
            ${CMAKE_SOURCE_DIR}/tests/lra_single_box.smt2
            --expected-polytopes 1
            --expected-rows 1
            --expected-volumes 400
            --expected-deleted 0
            --volume-tolerance 80
            --final-range 360 440)

    add_test(
        NAME lra_volume_disjoint
        COMMAND
            ${Python3_EXECUTABLE}
            ${LRA_VOLUME_TEST_SCRIPT}
            $<TARGET_FILE:ttc>
            ${CMAKE_SOURCE_DIR}/tests/lra_disjoint.smt2
            --expected-polytopes 2
            --expected-rows 2
            --expected-volumes 400
            --expected-deleted 0 0
            --volume-tolerance 120
            --final-range 720 880)

    add_test(
        NAME lra_volume_redundant
        COMMAND
            ${Python3_EXECUTABLE}
            ${LRA_VOLUME_TEST_SCRIPT}
            $<TARGET_FILE:ttc>
            ${CMAKE_SOURCE_DIR}/tests/lra_redundant.smt2
            --expected-polytopes 3
            --expected-rows 3
            --expected-volumes 100 40 0
            --expected-deleted 0 44 0
            --volume-tolerance 60
            --final-range 70 140)
endif()
