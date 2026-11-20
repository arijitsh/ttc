# Standalone volesti unit tests.  Every test_*.cpp becomes bin/<name> (doctest).
#   make -C tests/volesti -j8          # build all
#   make -C tests/volesti run          # build + run all
#   make -C tests/volesti bin/test_smoke && tests/volesti/bin/test_smoke
# OVERLAY=<dir> puts <dir> before src/volesti/include, to try patched headers.
ROOT     := $(abspath ../..)
CXX      ?= g++
OVERLAY  ?=
CXXFLAGS ?= -O2 -std=gnu++17 -DNDEBUG -w
INCS     := $(if $(OVERLAY),-I$(OVERLAY)) -I$(ROOT)/src/volesti/include -I$(ROOT)/src -I$(ROOT)/src/ddnnf \
            -I$(ROOT)/src/volesti/external -I$(ROOT)/src/volesti/include/generators \
            -I/usr/include/eigen3 -I/usr/include/lpsolve -I$(ROOT)/src/volesti/test
LIBS     := /usr/lib/liblpsolve55.a -lcolamd -lgmpxx -lgmp -ldl -lm -lpthread
SRCS     := $(wildcard test_*.cpp)
BINS     := $(patsubst %.cpp,bin/%,$(SRCS))

# The cvc5-free parts of ttc's LRA volume engine (src/volume/README.md).
TTC_NAMES := polytope_balls polytope_volume polytope_sampler union_parameters get_precision
TTC_OBJS := $(patsubst %,bin/%.o,$(TTC_NAMES))
TTC_INCS := -I$(ROOT)/src/ddnnf
# Shared objects of ttc code: keep them (pattern-rule intermediates are
# otherwise deleted after a build, breaking concurrent links).
.SECONDARY: $(TTC_OBJS) bin/doctest_main.o

all: $(BINS)
bin/doctest_main.o: doctest_main.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -I$(ROOT)/src/volesti/test -c $< -o $@
$(TTC_OBJS): bin/%.o: $(ROOT)/src/volume/%.cpp $(ROOT)/src/volume/*.hpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCS) $(TTC_INCS) -c $< -o $@
bin/%: %.cpp common.hpp bin/doctest_main.o $(TTC_OBJS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCS) $< bin/doctest_main.o $(TTC_OBJS) -o $@ $(LIBS)
run: $(BINS)
	@fail=0; for t in $(BINS); do echo "== $$t"; ./$$t || fail=1; done; exit $$fail
clean:
	rm -rf bin
.PHONY: all run clean
