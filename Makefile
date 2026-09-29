CXX       = g++
CXXFLAGS  = -O2 -Wall

# ---- PYTHIA (from environment variable) ----
ifeq ($(filter clean,$(MAKECMDGOALS)),)
ifndef PYTHIA8
$(error PYTHIA8 environment variable is not set)
endif
endif

PYINC     = -I$(PYTHIA8)/include
PYLIB     = -L$(PYTHIA8)/lib -lpythia8 -ldl
PYRPATH   = -Wl,-rpath,$(PYTHIA8)/lib

# ---- ROOT (auto-detected via root-config; its cflags already set the C++ standard) ----
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
ROOTRPATH  := -Wl,-rpath,$(shell root-config --libdir)

# ---- Sources: every .cc in the repo, at any depth ----
SOURCES := $(shell find . -name '*.cc')
TARGETS := $(SOURCES:.cc=)

.PHONY: all clean

all: $(TARGETS)

%: %.cc
	$(CXX) $< $(CXXFLAGS) \
	$(PYINC) $(ROOTCFLAGS) \
	$(PYLIB) $(ROOTLIBS) \
	$(PYRPATH) $(ROOTRPATH) \
	-o $@

clean:
	rm -f $(TARGETS)