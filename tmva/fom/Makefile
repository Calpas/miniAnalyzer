Target  = do.exe

Objects = h2taus.o 

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

#CXXLIBS    = -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -lPhysicsToolsUtilities 
CXXFLAGS   = -std=c++0x -I$(CMSSW_RELEASE_BASE)/src/

all:$(Target)

do.exe: main.cc $(Objects) /lstore/cms/calpas/h2taus/CMSSW_8_0_20/src/miniAOD/miniAnalyzer/func/myfunc.cc 
	g++ -o $@ main.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS) $(CXXFLAGS) #$(CXXLIBS)

h2taus.o: h2taus.C h2taus.h 
	g++ -o $@ h2taus.C -c $(ROOTFLAGS) $(CXXFLAGS)
        
clean:
	rm -f do.exe h2taus.o

info:
	@echo SCRAM_ARCH = $(SCRAM_ARCH)
	@echo CMSSW_RELEASE_BASE = $(CMSSW_RELEASE_BASE)
