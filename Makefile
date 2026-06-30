CXX=clang++-11
LD=clang++-11
CXXFLAGS=-std=c++14 -O2 -fPIC -fno-rtti

GSLCFG=gsl-config
GSLCFLAGS=$(shell $(GSLCFG) --cflags)
GSLLDFLAGS=$(shell $(GSLCFG) --libs)

CLANG=clang-11
CLANGFLAGS=-fPIC -fno-rtti

OMPFLAGS=-fopenmp
OS_NAME=$(shell uname -s)
ifeq ($(OS_NAME), Darwin)
	# For MacOS: https://iscinumpy.gitlab.io/post/omp-on-high-sierra/
	OMPFLAGS:=-Xpreprocessor $(OMPFLAGS)
endif

GREEN=$(shell tput setaf 2)
NOCOLOR=$(shell tput sgr0)

define MSG
	@tput setaf 2
	@echo [DONE]$1
	@tput sgr0
endef

BASIC_DEPS = src/fpInterface.h src/fpUtil.h
OBS_DEPS = src/local-sampler.h src/interval-sampler.h \
	src/observation.h src/local-observation.h \
	src/serializer.h src/polynomial.h
# DEP-LIBS = external-libs

.PHONY: clean utility targets gsl-fit testground oracles

all: utility targets gsl-fit testground oracles

# The floating-point utility start
utility: build/fpUtil.o
	$(call MSG, "The FP Utility Compiled (but not linked): $<")
build/fpUtil.o: src/fpUtil.cpp src/fpUtil.h
	@mkdir -p build
	$(CXX) -c $(CXXFLAGS) -o $@ $<
# The floating-point utility finish

# Compile the target under analysis start
targets: build/targets.o
	$(call MSG, "The Targets Compiled: $^")
build/targets.o: inputs/targets.cpp
	@mkdir -p build
	$(CLANG) -c $(CLANGFLAGS) $(GSLCFLAGS) -o $@ $<
# Compile the target under analysis finish

oracles: build/oracles.o
	$(call MSG, "The Oracles Compiled: $^")
build/oracles.o: inputs/oracles.cpp
	@mkdir -p build
	$(CLANG) -c $(CLANGFLAGS) $(GSLCFLAGS) -o $@ $<

testground: bin/testground.out
	$(call MSG, "The Test Ground")
bin/testground.out: build/testground.o build/fpUtil.o build/targets.o build/oracles.o build/herbie_optimized.o autornp_patch.o
	@mkdir -p bin
	$(LD) -o $@ $^ $(GSLLDFLAGS) -lmpfr -lgmp -L/usr/lib/llvm-10/lib -lomp -lalglib
build/testground.o: src/testground.cpp src/gsl-fit.h $(BASIC_DEPS) $(OBS_DEPS)
	@mkdir -p build
	$(CXX) -c $(CXXFLAGS) $(OMPFLAGS) $(GSLCFLAGS) -o $@ $<
build/herbie_optimized.o: herbie_optimized.c
	@mkdir -p build
	$(CLANG) -c $(CLANGFLAGS) -o $@ $<

gsl-fit: bin/gsl-fit.out
	$(call MSG, "The GSL Fit Binary")
bin/gsl-fit.out: build/gsl-fit.o build/fpUtil.o build/targets.o
	@mkdir -p bin
	@mkdir -p data
	$(LD) -o $@ $^ $(GSLLDFLAGS) -lmpfr -lgmp -L/usr/lib/llvm-10/lib -lomp -lalglib
build/gsl-fit.o: src/gsl-fit.cpp src/gsl-fit.h $(BASIC_DEPS) $(OBS_DEPS)
	@mkdir -p build
	$(CXX) -c $(CXXFLAGS) $(OMPFLAGS) $(GSLCFLAGS) -o $@ $<

runtime-cost: bin/runtime-cost.out
	$(call MSG, "Measure the Runtime of Patches")
build/runtime.o: src/measure-runtime.cpp patches/patches.h
	@mkdir -p build
	$(CXX) -c $(CXXFLAGS) $< -o $@
bin/runtime-cost.out: build/runtime.o build/targets.o
	@mkdir -p bin
	$(LD) -o $@ $^

sanity-check: bin/sanity-check.out
	$(call MSG, "Sanity Check")
bin/sanity-check.out: build/sanity-check.o
	@mkdir -p bin
	$(LD) -o $@ $^ -lmpfr -lgsl
build/sanity-check.o: src/sanity-check.cpp
	@mkdir -p build
	$(CXX) -c $(CXXFLAGS) $< -o $@

build/measure-runtime.o: src/measure-runtime.cpp patches/patches.h
	@mkdir -p build
	$(CXX) -c $(CXXFLAGS) $< -o $@
bin/measure-rumtime: build/measure-runtime.o build/targets.o
	@mkdir -p bin
	$(LD) -o $@ $^ -lgsl -lalglib

clean:
	find . -maxdepth 1 -name '*.o' ! -name 'autornp_patch.o' -delete
	rm -f *.so *.out
	rm -f *.ll *.bc *.s
	rm -f *.dwo
	rm -rf *.dSYM
	rm -rf lib/
	rm -rf build/
	rm -rf bin/
	rm -rf pdfs/
	rm -rf data/
	rm -rf tmp/
	rm -rf outputs/
	rm -rf figures/
	rm -rf patches/
