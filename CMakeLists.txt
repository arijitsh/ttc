cmake_minimum_required(VERSION 3.10)
project(ttc)

include(CTest)

# Directory holding the vendored solver dependencies (cvc5, cadical,
# cryptominisat, skolemfc). configure.sh clones+builds them here and passes
# -DTTC_DEPS_DIR. All solver paths below resolve under it, so ttc links ONLY
# against deps/ and never against a host-specific ~/solvers layout. System
# libraries (Boost, GMP, ISL, cddlib, lpsolve, SuiteSparse) stay as system pkgs.
if(NOT DEFINED TTC_DEPS_DIR)
  set(TTC_DEPS_DIR "${CMAKE_SOURCE_DIR}/deps")
endif()
set(TTC_DEPS_DIR "${TTC_DEPS_DIR}" CACHE PATH "Directory holding vendored dependencies")
message(STATUS "TTC dependency directory: ${TTC_DEPS_DIR}")

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
    if(APPLE)
        # Apple's toolchain has no static libc, so a fully `-static` binary is
        # not possible. On macOS BUILD_STATIC keeps the static-archive wiring
        # (build_static deps, single CMS, cadiback recompile) but links the .a
        # archives into an otherwise dynamic binary -- no `-static`, no GNU-only
        # objcopy/--exclude-libs (all guarded by `if(NOT APPLE)` below).
        message(STATUS "  (macOS: archives linked into a dynamic binary; not fully static)")
    else()
        # Linux: a genuinely static binary.
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
        # Prefer the static archive (.a) over the shared object so nothing sneaks
        # a .so into the -static link, and pull Boost from static archives too.
        # On macOS the binary is dynamic, so let system libs (cddlib, lpsolve,
        # isl, gmp, boost...) resolve to .dylib -- many Homebrew kegs ship no .a;
        # only the vendored solver deps (found by explicit path) stay .a.
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        set(Boost_USE_STATIC_LIBS ON)
    endif()
endif()

# Enable feature macro when DDNNF/D4 support is requested
if(ENABLE_DDNNF)
    add_compile_definitions(TTC_ENABLE_DDNNF)
endif()

# Locate cvc5 headers and libraries manually (vendored at deps/cvc5).
# NO_DEFAULT_PATH pins resolution to the vendored tree -- otherwise a stray
# system/${HOME} cvc5 install can satisfy the lookup and mask the build-tree
# include below (which only the vendored build has).
find_path(CVC5_INCLUDE_DIR NAMES cvc5/cvc5.h
          HINTS ${TTC_DEPS_DIR}/cvc5/include NO_DEFAULT_PATH)
find_library(CVC5_LIB NAMES cvc5
            HINTS ${TTC_DEPS_DIR}/cvc5/build/src NO_DEFAULT_PATH)
find_library(CVC5PARSER_LIB NAMES cvc5parser
            HINTS ${TTC_DEPS_DIR}/cvc5/build/src/parser NO_DEFAULT_PATH)
if(NOT CVC5_INCLUDE_DIR OR NOT CVC5_LIB OR NOT CVC5PARSER_LIB)
  message(FATAL_ERROR "cvc5 not found under ${TTC_DEPS_DIR}/cvc5. Run ./configure.sh to build it.")
endif()
include_directories(${CVC5_INCLUDE_DIR})
# cvc5 generates some public headers (cvc5_export.h, cvc5_proof_rule.h) into its
# build tree's include/, NOT the source include/. cvc5.h #includes
# <cvc5/cvc5_export.h>, so add the build-tree include dir too.
set(CVC5_BUILD_INCLUDE_DIR "${TTC_DEPS_DIR}/cvc5/build/include")
if(EXISTS "${CVC5_BUILD_INCLUDE_DIR}/cvc5/cvc5_export.h")
  include_directories(${CVC5_BUILD_INCLUDE_DIR})
else()
  message(WARNING "cvc5_export.h not found in ${CVC5_BUILD_INCLUDE_DIR}; cvc5 may not be built yet")
endif()
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
set(CVC5_DEPS_LIB_DIR "${TTC_DEPS_DIR}/cvc5/build/deps/lib")
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
# CaDiCaL: the vendored `xor` branch (deps/cadical) is the SINGLE cadical for
# the whole link. It exports Solver::add_xor(int) (cvc5's static libcvc5.a needs
# it) AND cadiback's extras Solver::traverse_red_clauses()/get_eqiv_lits(), so
# the same archive satisfies cvc5 AND SkolemFC's cadiback -- no second cadical in
# the link means no "multiple definition" collision.
set(NEW_CADICAL_ROOT "${TTC_DEPS_DIR}/cadical")
find_library(CADICAL_LIB cadical
             HINTS "${NEW_CADICAL_ROOT}/build"
             NO_DEFAULT_PATH)
if(NOT CADICAL_LIB)
  message(FATAL_ERROR "sc2025 CaDiCaL not found at ${NEW_CADICAL_ROOT}/build/libcadical.a")
endif()
message(STATUS "CaDiCaL library: ${CADICAL_LIB}")
set(CADICAL_SRC_INCLUDE_DIR "${NEW_CADICAL_ROOT}/src")
if(EXISTS "${CADICAL_SRC_INCLUDE_DIR}/internal.hpp")
  include_directories(${CADICAL_SRC_INCLUDE_DIR})
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

# SkolemFC: uninterpreted-function (Skolem-function) counting engine. Links the
# locally built libskolemfc and its CryptoMiniSat / ApproxMC / UniGen / Arjun
# dependencies. The dynamic build uses build_test/ (shared .so); a fully static
# ttc (BUILD_STATIC) uses build_static/ which holds the whole stack as .a
# archives (configure.sh builds it: cd deps/skolemfc && ./configure.sh --static
# --build-dir build_static && cmake --build build_static -j).
set(SKOLEMFC_ROOT "${TTC_DEPS_DIR}/skolemfc")
if(BUILD_STATIC)
  set(SKOLEMFC_BUILD "${SKOLEMFC_ROOT}/build_static")
else()
  set(SKOLEMFC_BUILD "${SKOLEMFC_ROOT}/build_test")
endif()
find_library(SKOLEMFC_LIB skolemfc HINTS "${SKOLEMFC_BUILD}/src" NO_DEFAULT_PATH)
set(SKOLEMFC_DEPS_LIB_DIR "${SKOLEMFC_BUILD}/deps_install/lib")
set(SKOLEMFC_DEPS_INCLUDE_DIR "${SKOLEMFC_BUILD}/deps_install/include")
if(SKOLEMFC_LIB AND EXISTS "${SKOLEMFC_ROOT}/src/skolemfc.h")
  message(STATUS "Found SkolemFC library: ${SKOLEMFC_LIB}")
  find_library(SKOLEMFC_CMS_LIB cryptominisat5 HINTS ${SKOLEMFC_DEPS_LIB_DIR} NO_DEFAULT_PATH)
  find_library(SKOLEMFC_APPROXMC_LIB approxmc HINTS ${SKOLEMFC_DEPS_LIB_DIR} NO_DEFAULT_PATH)
  find_library(SKOLEMFC_UNIGEN_LIB unigen HINTS ${SKOLEMFC_DEPS_LIB_DIR} NO_DEFAULT_PATH)
  find_library(SKOLEMFC_ARJUN_LIB arjun HINTS ${SKOLEMFC_DEPS_LIB_DIR} NO_DEFAULT_PATH)
  if(BUILD_STATIC)
    # Static cryptominisat pulls in louvain-community, sbva, cadical and cadiback
    # symbols; the shared build hid these inside libcryptominisat5.so, but a
    # -static link must name every archive explicitly.
    find_library(SKOLEMFC_LOUVAIN_LIB louvain_communities HINTS ${SKOLEMFC_DEPS_LIB_DIR} NO_DEFAULT_PATH)
    find_library(SKOLEMFC_SBVA_LIB sbva HINTS ${SKOLEMFC_DEPS_LIB_DIR} NO_DEFAULT_PATH)
    # Use the SAME (xor-branch) CaDiCaL as cvc5 (CADICAL_LIB) for cadiback too.
    # The bundled SkolemFC cadical is v2.1.3; linking it alongside ours pulls
    # both archives' elim.o and triggers "multiple definition". The xor-branch
    # cadical exports cadiback's traverse_red_clauses()/get_eqiv_lits() extras,
    # so the single archive resolves cadiback's CaDiCaL refs with no second
    # cadical in the link -- mirroring the single-CMS fix below.
    set(SKOLEMFC_CADICAL_LIB "${CADICAL_LIB}")
    # The shipped libcadiback.a was compiled against the v2.1.3 cadical headers;
    # linking it against the xor-branch lib is an ABI mismatch (heap corruption
    # in CadiBack::doit). Recompile cadiback.cpp against the xor-branch headers
    # so its inlines/struct layouts match the single cadical we link.
    set(CADIBACK_SRC "${SKOLEMFC_ROOT}/deps/cadiback/cadiback.cpp")
    set(CADIBACK_XOR_LIB "${CMAKE_BINARY_DIR}/libcadiback_xor.a")
    add_custom_command(
        OUTPUT "${CADIBACK_XOR_LIB}"
        COMMAND ${CMAKE_CXX_COMPILER} -fPIC -O3 -std=c++17 -DNDEBUG
                -I${NEW_CADICAL_ROOT}/src -c "${CADIBACK_SRC}"
                -o "${CMAKE_BINARY_DIR}/cadiback_xor.o"
        COMMAND ${CMAKE_AR} rcs "${CADIBACK_XOR_LIB}" "${CMAKE_BINARY_DIR}/cadiback_xor.o"
        DEPENDS "${CADIBACK_SRC}" "${NEW_CADICAL_ROOT}/src/cadical.hpp"
        VERBATIM
        COMMENT "Recompiling cadiback against xor-branch CaDiCaL headers")
    add_custom_target(cadiback_xor DEPENDS "${CADIBACK_XOR_LIB}")
    set(SKOLEMFC_CADIBACK_LIB "${CADIBACK_XOR_LIB}")
    # CryptoMiniSat collision: cvc5 bundles its own libcryptominisat5.a and so
    # does SkolemFC; a single -static link would see "multiple definition" on
    # the shared CCNR/CMSat objects. Resolve by linking exactly ONE copy and
    # dropping cvc5's bundled one. We use SkolemFC's own built CMS
    # (deps_install/lib/libcryptominisat5.a, pinned commit b09bd6b): it is the
    # version arjun/approxmc/unigen and cadiback were all built against (6-arg
    # CadiBack::doit, translate_sampl_set/update_var_mapping present), so the
    # whole SkolemFC stack is ABI-consistent with it. cvc5 needs only the stable
    # CMSat::SATSolver facade methods, which this 5.13 archive also provides.
    set(SHARED_CMS_LIB "${SKOLEMFC_DEPS_LIB_DIR}/libcryptominisat5.a")
    if(NOT EXISTS "${SHARED_CMS_LIB}")
      message(FATAL_ERROR "SkolemFC CryptoMiniSat static lib not found at ${SHARED_CMS_LIB}; "
                          "run ./configure.sh to build the SkolemFC stack.")
    endif()
    set(CMS_LIB "")                       # drop cvc5's bundled copy
    set(SKOLEMFC_CMS_LIB "${SHARED_CMS_LIB}")  # single CMS for the whole link
  endif()
  set(HAVE_SKOLEMFC ON)
  add_compile_definitions(TTC_ENABLE_SKOLEMFC)
  include_directories(${PROJECT_SOURCE_DIR}/src/skolemfc)
  include_directories(${SKOLEMFC_ROOT}/src)
  include_directories(${SKOLEMFC_DEPS_INCLUDE_DIR})
else()
  message(WARNING "SkolemFC library not found at ${SKOLEMFC_BUILD}/src; "
                  "--skolemfc engine disabled")
  set(HAVE_SKOLEMFC OFF)
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
    if(BUILD_STATIC AND NOT APPLE)
        # liblpsolve55.a and libcdd.a both export the C symbol `set_int` (two
        # unrelated functions). The shared build tolerated it (separate .so);
        # a -static link reports "multiple definition". Nothing outside lpsolve
        # references lpsolve's set_int, so rewrite that symbol in a patched copy
        # of the archive, leaving cddlib's set_int as the sole definition.
        # (GNU objcopy only; macOS lacks it -- the dynamic macOS link tolerates
        # the duplicate the same way a shared Linux build did.)
        set(LPSOLVE_PATCHED "${CMAKE_BINARY_DIR}/liblpsolve55_patched.a")
        add_custom_command(
            OUTPUT "${LPSOLVE_PATCHED}"
            COMMAND ${CMAKE_OBJCOPY} --redefine-sym set_int=lpsolve_set_int
                    "${LPSOLVE_LIB}" "${LPSOLVE_PATCHED}"
            DEPENDS "${LPSOLVE_LIB}"
            VERBATIM
            COMMENT "Patching liblpsolve55.a: rename set_int -> lpsolve_set_int (avoids cddlib clash)")
        add_custom_target(lpsolve_patched DEPENDS "${LPSOLVE_PATCHED}")
        set(LPSOLVE_LIB "${LPSOLVE_PATCHED}")
    endif()
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
          HINTS $ENV{EIGEN3_INCLUDE_DIR} /usr/include/eigen3
          PATH_SUFFIXES eigen3)   # Homebrew/others install under <prefix>/include/eigen3
if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Eigen3 headers not found. Please install libeigen3-dev")
endif()
include_directories(${EIGEN3_INCLUDE_DIR})

set(SOURCES
    src/main.cpp
    src/parser.cpp
    src/ddnnf/arjun.cpp
    src/pact/pact.cpp
    src/pact/saturatingcount.cpp
    src/pact/hash_constraint.cpp
    src/pact/bool_hash.cpp
    src/pact/bv_hash.cpp
    src/projected_enumerator.cpp
    src/tdd/tdd.cpp
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

if(HAVE_SKOLEMFC)
    list(APPEND SOURCES src/skolemfc/skolemfc_engine.cpp)
endif()

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

if(NOT BUILD_STATIC AND UNIX AND NOT APPLE)
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.so.6
        OUTPUT_VARIABLE TTC_LIBSTDCXX
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(TTC_LIBSTDCXX AND IS_ABSOLUTE "${TTC_LIBSTDCXX}")
        get_filename_component(TTC_LIBSTDCXX_DIR "${TTC_LIBSTDCXX}" DIRECTORY)
        get_filename_component(TTC_LIBSTDCXX_DIR "${TTC_LIBSTDCXX_DIR}" REALPATH)
        message(STATUS "C++ runtime RPATH: ${TTC_LIBSTDCXX_DIR}")
        set_target_properties(ttc PROPERTIES
            BUILD_RPATH "${TTC_LIBSTDCXX_DIR}"
            INSTALL_RPATH "${TTC_LIBSTDCXX_DIR}")
        # Use DT_RPATH instead of DT_RUNPATH so an older libstdc++ injected by
        # LD_LIBRARY_PATH (for example Anaconda) cannot override the compiler
        # runtime needed by this executable and its solver libraries.
        target_link_options(ttc PRIVATE "LINKER:--disable-new-dtags")
    endif()
endif()

# The static cryptominisat/louvain stack calls OpenMP runtime routines
# (omp_get_wtime, etc.); a -static link needs libgomp named explicitly.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(ttc OpenMP::OpenMP_CXX)
endif()

# The static lpsolve archive is patched (set_int renamed) by a custom command;
# make ttc wait for it.
if(BUILD_STATIC AND TARGET lpsolve_patched)
    add_dependencies(ttc lpsolve_patched)
endif()
if(BUILD_STATIC AND TARGET cadiback_xor)
    add_dependencies(ttc cadiback_xor)
endif()

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

if(HAVE_SKOLEMFC)
    target_link_libraries(ttc
        ${SKOLEMFC_LIB}
        ${SKOLEMFC_APPROXMC_LIB}
        ${SKOLEMFC_UNIGEN_LIB}
        ${SKOLEMFC_ARJUN_LIB}
        ${SKOLEMFC_CMS_LIB})
    if(BUILD_STATIC)
        # Static link: name every transitive archive (cryptominisat -> cadical,
        # cadiback, sbva, louvain-community) explicitly, in dependency order.
        # gmpxx/gmp are repeated here (after libskolemfc) so the GMP C++ stream
        # operators it references resolve, and libsuitesparseconfig satisfies
        # colamd's SuiteSparse_config_* symbols in a -static link.
        find_library(SUITESPARSECONFIG_LIB suitesparseconfig)
        target_link_libraries(ttc
            ${SKOLEMFC_CADIBACK_LIB}
            ${SKOLEMFC_SBVA_LIB}
            ${SKOLEMFC_CADICAL_LIB}
            ${SKOLEMFC_LOUVAIN_LIB}
            ${SUITESPARSECONFIG_LIB}
            ${GMPXX_LIB}
            ${GMP_LIB}
            ${CMAKE_DL_LIBS})   # `dl` on Linux, empty on macOS (in libSystem)
    else()
        # Resolve libskolemfc and its CryptoMiniSat/ApproxMC/UniGen deps at runtime
        # without polluting LD_LIBRARY_PATH.
        set_property(TARGET ttc APPEND PROPERTY
            BUILD_RPATH "${SKOLEMFC_BUILD}/src;${SKOLEMFC_DEPS_LIB_DIR}")
        set_property(TARGET ttc APPEND PROPERTY
            INSTALL_RPATH "${SKOLEMFC_BUILD}/src;${SKOLEMFC_DEPS_LIB_DIR}")
    endif()
    # cvc5 is statically linked and bundles an OLDER CryptoMiniSat. Without this,
    # those CMSat symbols are exported in ttc's dynamic symbol table and the
    # dynamic linker interposes them into SkolemFC's libapproxmc /
    # libcryptominisat5.so.5.13 (a different, incompatible ABI), crashing inside
    # ApproxMC. Hiding all static-archive symbols from the executable's dynamic
    # symbol table lets the SkolemFC shared libraries resolve CryptoMiniSat
    # among themselves. GNU-ld only; macOS ld64 has no equivalent (and the
    # single-static-CMS design has no shared SkolemFC libs to interpose into).
    if(NOT APPLE)
        target_link_options(ttc PRIVATE "LINKER:--exclude-libs,ALL")
    endif()
endif()

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
