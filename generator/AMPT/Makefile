# Hybrid High-Performance Makefile for AMPT v1.26t9b-v2.26t9b + CERN ROOT
# Native CPU optimizations and LTO with IEEE-exact arithmetic: no -ffast-math,
# and -ffp-contract=off so that -march=native does not fuse a*b+c into FMA
# instructions (which round differently and change every event). With it the
# output is the same as from the original code built with gfortran -O.
FC       = gfortran
# compiler driver for the final link (e.g. LINK=gfortran-14 if the default
# gfortran installation cannot link)
LINK    ?= $(FC)
CXX      = g++

# LTO= (empty) when the compiler has no LTO linker plugin (liblto_plugin.so)
LTO     ?= -flto
OPT      = -O3 -march=native -mtune=native -funroll-loops $(LTO) -ffp-contract=off
FFLAGS   = $(OPT) -fno-automatic -finit-local-zero \
           -std=legacy -Wno-legacy-extensions -Wno-argument-mismatch
CXXFLAGS = $(OPT) $(shell root-config --cflags)
# extra link options, e.g. LDEXTRA=-L/usr/lib/gcc/x86_64-linux-gnu/15 (libstdc++ of g++)
LDEXTRA ?=
LDFLAGS  = $(OPT) $(LDEXTRA) $(shell root-config --libs) -lEG -lstdc++

FSRCS    = main.f amptsub.f art1f.f hijing1.383_ampt.f hipyset1.35.f zpc.f linana.f
FOBJS    = $(FSRCS:.f=.o)
CXXOBJS  = AmptRootIO.o

all: ampt

ampt: $(FOBJS) $(CXXOBJS)
	$(LINK) -o $@ $(FOBJS) $(CXXOBJS) $(LDFLAGS)

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

AmptRootIO.o: AmptRootIO.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Profile-guided build (~10% faster): instrumented build, PGO_EVENTS events of
# the system in input.ampt as training, rebuild. Same arithmetic (no fast-math).
PGO_EVENTS ?= 5
pgo:
	$(MAKE) clean
	$(MAKE) ampt OPT="$(OPT) -fprofile-generate -fprofile-update=single"
	rm -rf pgo-train && mkdir -p pgo-train/ana
	sed -e 's/^[0-9]*\(\s*! NEVNT\)/$(PGO_EVENTS)\1/' -e 's/^[0-9]*\(\s*! ihjsed\)/11\1/' \
	    input.ampt > pgo-train/input.ampt
	cd pgo-train && printf '123456789\n8\n' | AMPT_ROOT_OUT=ana/train.root ../ampt > train.log 2>&1
	rm -f *.o ampt
	$(MAKE) ampt OPT="$(OPT) -fprofile-use -fprofile-partial-training -Wno-missing-profile"
	rm -rf pgo-train *.gcda

clean:
	rm -f *.o *.gcda ampt core

.PHONY: all pgo clean
