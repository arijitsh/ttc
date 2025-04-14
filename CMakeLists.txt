cmake_minimum_required(VERSION 3.10)
project(ttc)

# Option to build a static binary (default is OFF)
option(BUILD_STATIC "Build static binary" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(BUILD_STATIC)
    message(STATUS "Building a static binary")
    # Linux-specific flag for static linking; adjust for your environment.
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif()

# Find cvc5. Note: Depending on your cvc5 installation, find_package(cvc5)
# may set variables for both the core and parser libraries.
find_package(cvc5 REQUIRED)
if(cvc5_FOUND)
  message(STATUS "Found cvc5 in ${cvc5_DIR}, include dir: ${cvc5_INCLUDE_DIRS}, lib dir: ${cvc5_LIBRARY_DIRS}")
  include_directories(${cvc5_INCLUDE_DIRS})
else()
  message(WARNING "cvc5 not found via find_package. Please set include directories manually.")
  # Example:
  # include_directories(/path/to/cvc5/include)
  # link_directories(/path/to/cvc5/lib)
endif()

# Find Boost (we need program_options)
find_package(Boost REQUIRED COMPONENTS program_options)
if(Boost_FOUND)
    message(STATUS "Found Boost: ${Boost_INCLUDE_DIRS}")
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Boost not found!")
endif()

# Add our project source directory
include_directories(${PROJECT_SOURCE_DIR}/src)

set(SOURCES
    src/main.cpp
    src/parser.cpp
)

add_executable(ttc ${SOURCES})

# Link against both the cvc5 core library and its parser library along with Boost.
# Adjust the library names if your installation uses different identifiers.
target_link_libraries(ttc cvc5parser cvc5 ${Boost_LIBRARIES})
