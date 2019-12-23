# $Id: GNUmakefile 68058 2013-03-13 14:47:43Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -g
CXXFLAGS += -g
CXXFLAGS += $(shell root-config --cflags)
CPPFLAGS += $(shell root-config --cflags) -Wl,--no-as-needed
LDFLAGS += $(shell root-config --libs --glibs)

name := FCP_Simulation
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

