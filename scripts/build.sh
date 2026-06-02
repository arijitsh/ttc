#!/usr/bin/env bash
set -e

clean=0
while getopts "c" opt; do
  case $opt in
    c) clean=1 ;;
  esac
done

if [ $clean -eq 1 ]; then
  rm -rf CMakeFiles
  rm -rf .cmake
  rm -f CMakeCache.txt Makefile cmake_install.cmake libd4lib.a ttc
  echo "Cleaned build artifacts."
fi

cmake -DENABLE_DEBUG=ON ..
make -j"$(nproc)"
