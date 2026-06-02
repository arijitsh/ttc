#!/usr/bin/env bash
set -e

BUILD_DIR="build"
USE_LOCAL=false
DEBUG=false
ENABLE_DDNNF=false
STATIC=false
CMAKE_ARGS=()

# Parse arguments
for arg in "$@"; do
  case $arg in
    --local)
      USE_LOCAL=true
      ;;
    --debug)
      DEBUG=true
      ;;
    --ddnnf)
      ENABLE_DDNNF=true
      ;;
    --static)
      STATIC=true
      ;;
    *)
      echo "Usage: $0 [--local] [--debug] [--ddnnf] [--static]" >&2
      exit 1
      ;;
  esac
done

CMAKE_ARGS+=("-DCMAKE_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu")
CMAKE_ARGS+=("-DCMAKE_IGNORE_PATH=$HOME/anaconda3")

if $DEBUG; then
  CMAKE_ARGS+=("-DENABLE_DEBUG=ON")
fi

if $STATIC; then
  CMAKE_ARGS+=("-DBUILD_STATIC=ON")
else
  CMAKE_ARGS+=("-DBUILD_STATIC=OFF")
fi

if $ENABLE_DDNNF; then
  CMAKE_ARGS+=("-DENABLE_DDNNF=ON")
else
  CMAKE_ARGS+=("-DENABLE_DDNNF=OFF")
fi

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

if $USE_LOCAL; then
  CVC5_ROOT="$HOME/solvers/cvc5_solver"
  CVC5_LIB="$CVC5_ROOT/build/src/libcvc5.so.1"
  CVC5_PARSER_LIB="$CVC5_ROOT/build/src/libcvc5parser.so.1"
  CVC5_INCLUDE="$CVC5_ROOT/include"

  cmake .. "${CMAKE_ARGS[@]}" \
    -DCVC5_LIB="$CVC5_LIB" \
    -DCVC5PARSER_LIB="$CVC5_PARSER_LIB" \
    -DCVC5_INCLUDE_DIR="$CVC5_INCLUDE"

  echo "-- Configured with local cvc5 from $CVC5_ROOT"
else
  cmake .. "${CMAKE_ARGS[@]}"
  echo "-- Configured with system-installed cvc5"
fi

if $ENABLE_DDNNF; then
  echo "-- DDNNF/D4 support enabled"
else
  echo "-- DDNNF/D4 support disabled (PACT default)"
fi

if $STATIC; then
  echo "-- Building static library"
else
  echo "-- Building shared library"
fi

if $DEBUG; then
  echo "-- Debug symbols enabled"
else
  echo "-- Debug symbols disabled"
fi

echo "-- Build configuration complete. To compile run:"
echo "  cd $BUILD_DIR && make"
