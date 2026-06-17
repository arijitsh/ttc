#!/usr/bin/env bash
#
# configure.sh -- fetch ttc's solver dependencies into ./deps and configure the
# build to link ONLY against them (cvc5-style vendoring). System libraries
# (Boost, GMP, ISL, cddlib, lpsolve, SuiteSparse, OpenMP) are NOT vendored; they
# are resolved as ordinary system packages, exactly as cvc5 does.
#
# Vendored dependencies (cloned + built under deps/):
#   deps/cvc5          arijitsh/cvc5        @ ttc-cpp   (static libcvc5 + parser)
#   deps/cadical       arijitsh/cadical     @ xor       (XOR-reasoning CaDiCaL)
#   deps/cryptominisat msoos/cryptominisat  @ master    (single shared CMS 5.x)
#   deps/skolemfc      meelgroup/skolemfc   @ devel     (UF counting + ApproxMC/Arjun)
#
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEPS_DIR="$SCRIPT_DIR/deps"
BUILD_DIR="$SCRIPT_DIR/build"
JOBS="$(nproc)"

# ----------------------------------------------------------------------------
# Dependency pins (repo + branch + exact commit). Keep in sync with CMakeLists
# TTC_DEPS_DIR. Commits are pinned for reproducibility: floating branch HEADs
# drift and break the cross-dependency ABI (e.g. msoos/cryptominisat master
# moved to an 8-arg CadiBack::doit incompatible with skolemfc's cadiback).
# ----------------------------------------------------------------------------
CVC5_REPO="https://github.com/arijitsh/cvc5.git";          CVC5_BRANCH="ttc-cpp"
CVC5_COMMIT="445cb7a0d9b0e8c01eecf01815ba958c7c42bfc3"
CADICAL_REPO="https://github.com/arijitsh/cadical.git";    CADICAL_BRANCH="xor"
CADICAL_COMMIT="52592df3a73a5cb2ce3387bd3bee6c1a78dda12f"
# CryptoMiniSat is NOT cloned separately: ttc links the single CMS that the
# SkolemFC build produces (deps/skolemfc/.../libcryptominisat5.a, SkolemFC's
# pinned deps/cryptominisat submodule). That copy is ABI-consistent with the
# rest of the SkolemFC stack (arjun/approxmc/cadiback) and also serves cvc5's
# CMSat::SATSolver facade -- see CMakeLists. A standalone msoos/cryptominisat
# HEAD drifts (8-arg CadiBack::doit) and cannot even build alone (needs cadiback).
SKOLEMFC_REPO="https://github.com/meelgroup/skolemfc.git"; SKOLEMFC_BRANCH="devel"
SKOLEMFC_COMMIT="4b1889f6ee23d9309ac87b59ea4340e3c9313715"

# ----------------------------------------------------------------------------
# Options
# ----------------------------------------------------------------------------
DEBUG=false
ENABLE_DDNNF=false
STATIC=true          # vendored deps are built static; default to a static ttc
SKIP_DEPS=false      # --skip-deps: don't touch deps/, just (re)configure cmake
REBUILD_DEPS=false   # --rebuild-deps: force rebuild even if outputs exist
CMAKE_ARGS=()

usage() {
  echo "Usage: $0 [--debug] [--ddnnf] [--shared] [--skip-deps] [--rebuild-deps]" >&2
  echo "  (deps are fetched into ./deps and the build links only against them)" >&2
  exit 1
}

for arg in "$@"; do
  case $arg in
    --debug)        DEBUG=true ;;
    --ddnnf)        ENABLE_DDNNF=true ;;
    --static)       STATIC=true ;;     # kept for compatibility (default)
    --shared)       STATIC=false ;;
    --skip-deps)    SKIP_DEPS=true ;;
    --rebuild-deps) REBUILD_DEPS=true ;;
    -h|--help)      usage ;;
    *)              echo "Unknown option: $arg" >&2; usage ;;
  esac
done

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
log()  { echo -e "\033[1;34m[deps]\033[0m $*"; }
warn() { echo -e "\033[1;33m[deps]\033[0m $*" >&2; }

# clone_dep <name> <repo> <branch> <commit>
clone_dep() {
  local name="$1" repo="$2" branch="$3" commit="$4"
  local dir="$DEPS_DIR/$name"
  if [ -d "$dir/.git" ]; then
    log "$name already cloned"
  else
    log "cloning $name ($repo @ $branch)"
    git clone --branch "$branch" "$repo" "$dir"
  fi
  # Pin to the exact reproducible commit.
  if [ -n "$commit" ] && [ "$(git -C "$dir" rev-parse HEAD)" != "$commit" ]; then
    log "checking out $name @ $commit"
    git -C "$dir" fetch --quiet origin "$commit" 2>/dev/null || true
    git -C "$dir" checkout --quiet "$commit"
  fi
}

# ----------------------------------------------------------------------------
# Fetch + build dependencies
# ----------------------------------------------------------------------------
build_cadical() {
  local dir="$DEPS_DIR/cadical"
  if [ -f "$dir/build/libcadical.a" ] && ! $REBUILD_DEPS; then
    log "cadical already built (deps/cadical/build/libcadical.a)"; return
  fi
  log "building cadical (xor branch, static)"
  cmake -S "$dir" -B "$dir/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF
  cmake --build "$dir/build" -j"$JOBS"
}

build_skolemfc() {
  local dir="$DEPS_DIR/skolemfc"
  if [ -f "$dir/build_static/src/libskolemfc.a" ] && ! $REBUILD_DEPS; then
    log "skolemfc already built (deps/skolemfc/build_static/src/libskolemfc.a)"; return
  fi
  log "building skolemfc + its deps (approxmc/arjun/unigen/sbva/louvain/cadiback, static)"
  # Initialise only the dependency submodules skolemfc needs into
  # deps/skolemfc/deps/* (CMakeLists expects cadiback's source at
  # SKOLEMFC_ROOT/deps/cadiback). A fresh clone has them uninitialised, so
  # skolemfc's configure errors without this. We list them explicitly to avoid
  # the unrelated benchmark-baseline submodules (e.g. GPMC) and the nested
  # test-only submodules a recursive init would pull.
  local skolemfc_submodules=(
    deps/louvain-community deps/cryptominisat deps/arjun deps/approxmc
    deps/unigen deps/cadical deps/cadiback deps/sbva)
  ( cd "$dir" && git submodule update --init "${skolemfc_submodules[@]}" )
  # Apply ttc's library-interface patch (apply_default_config / get_count_log2 /
  # count_succeeded) -- these live in a ttc-local patch because they are not yet
  # upstream in meelgroup/skolemfc. Idempotent: skip if already applied.
  local patch="$SCRIPT_DIR/patches/skolemfc-library-interface.patch"
  if [ -f "$patch" ]; then
    if git -C "$dir" apply --reverse --check "$patch" 2>/dev/null; then
      log "skolemfc library-interface patch already applied"
    else
      log "applying skolemfc library-interface patch"
      git -C "$dir" apply "$patch"
    fi
  fi
  ( cd "$dir" && ./configure.sh --static --build-dir build_static )
  cmake --build "$dir/build_static" -j"$JOBS"
}

build_cvc5() {
  local dir="$DEPS_DIR/cvc5"
  if [ -f "$dir/build/src/libcvc5.a" ] && ! $REBUILD_DEPS; then
    log "cvc5 already built (deps/cvc5/build/src/libcvc5.a)"; return
  fi
  log "building cvc5 (static, auto-download deps) -- this is the long pole"
  # --cryptominisat: cvc5 needs the CryptoMiniSat bv-sat-solver backend (the
  # SkolemFC determinism check drives cvc5's BV SAT through it). cvc5's bundled
  # CMS .a is dropped at ttc link time; only its CMSat::SATSolver facade is used,
  # resolved against the single SkolemFC CMS (see CMakeLists).
  ( cd "$dir" && ./configure.sh --static --auto-download --cryptominisat --name=build )
  cmake --build "$dir/build" -j"$JOBS"
}

if ! $SKIP_DEPS; then
  mkdir -p "$DEPS_DIR"
  clone_dep cvc5          "$CVC5_REPO"     "$CVC5_BRANCH"     "$CVC5_COMMIT"
  clone_dep cadical       "$CADICAL_REPO"  "$CADICAL_BRANCH"  "$CADICAL_COMMIT"
  clone_dep skolemfc      "$SKOLEMFC_REPO" "$SKOLEMFC_BRANCH" "$SKOLEMFC_COMMIT"

  # Build lighter deps first so a cvc5 failure doesn't waste the others.
  # (CryptoMiniSat is built as part of the SkolemFC stack -- see build_skolemfc.)
  build_cadical
  build_skolemfc
  build_cvc5
  log "all dependencies built under $DEPS_DIR"
else
  log "--skip-deps: leaving deps/ untouched"
fi

# ----------------------------------------------------------------------------
# Configure ttc to link only against deps/
# ----------------------------------------------------------------------------
CMAKE_ARGS+=("-DTTC_DEPS_DIR=$DEPS_DIR")
CMAKE_ARGS+=("-DCMAKE_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu")
CMAKE_ARGS+=("-DCMAKE_IGNORE_PATH=$HOME/anaconda3")
# Force system Boost 1.83 static (not anaconda's), GMP archives for -static.
CMAKE_ARGS+=("-DBoost_NO_BOOST_CMAKE=ON" "-DBoost_NO_SYSTEM_PATHS=ON")
CMAKE_ARGS+=("-DBOOST_ROOT=/usr" "-DBOOST_LIBRARYDIR=/usr/lib/x86_64-linux-gnu")
CMAKE_ARGS+=("-DBoost_DIR=Boost_DIR-NOTFOUND")

$DEBUG        && CMAKE_ARGS+=("-DENABLE_DEBUG=ON")  || CMAKE_ARGS+=("-DENABLE_DEBUG=OFF")
$ENABLE_DDNNF && CMAKE_ARGS+=("-DENABLE_DDNNF=ON")  || CMAKE_ARGS+=("-DENABLE_DDNNF=OFF")
if $STATIC; then
  CMAKE_ARGS+=("-DBUILD_STATIC=ON")
  CMAKE_ARGS+=("-DGMP_LIB=/usr/lib/x86_64-linux-gnu/libgmp.a")
  CMAKE_ARGS+=("-DGMPXX_LIB=/usr/lib/x86_64-linux-gnu/libgmpxx.a")
else
  CMAKE_ARGS+=("-DBUILD_STATIC=OFF")
fi

mkdir -p "$BUILD_DIR"
# Wipe cached find_library results so paths resolve against deps/ cleanly.
rm -f "$BUILD_DIR/CMakeCache.txt"
cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" "${CMAKE_ARGS[@]}"

echo "-- Configured ttc to link only against $DEPS_DIR"
$STATIC && echo "-- Static build" || echo "-- Shared build"
$ENABLE_DDNNF && echo "-- DDNNF/D4 support enabled" || echo "-- DDNNF/D4 support disabled"
echo "-- Build configuration complete. To compile run:"
echo "  cmake --build $BUILD_DIR -j$JOBS"
