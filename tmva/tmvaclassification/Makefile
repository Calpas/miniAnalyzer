Target  = do.exe

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
CXXFLAGS  = -std=c++0x -I$(CMSSW_RELEASE_BASE)/src/
TMVAFLAGS = -lTMVA -lTMVAGui
GCCFLAG   = -Wall -Wextra

all:$(Target)

do.exe: main.cc ../../func/myfunc.cc
	g++ -o $@ main.cc $(ROOTFLAGS) $(ROOTLIBS) $(CXXFLAGS) $(TMVAFLAGS) $(GCCFLAG)

