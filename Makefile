# branchrates Makefile

CPP  = g++
CC   = gcc

CXXINCLUDEDIRS =
CXXFLAGS = $(CXXINCLUDEDIRS) -Wall -g3 -ggdb -fno-inline
CXXLDFLAGS = -g3 -ggdb

OBJ  = BranchRateManager.o \
	   HistorySingleParameter.o \
	   Likelihood.o \
	   RatePMap.o \
	   main.o \
	   ML_multi.o \
	   ML_multi_DownhillSimplex.o \
	   ML_multi_Powell.o \
	   ML_multi_QuasiNewton.o \
	   ML_single.o \
	   ML_single_NewtonRaphson.o \
	   PhyloTree.o \
	   PhyloTreeNode.o \
	   RatePVector.o \
	   Simulate.o \
	   TaxonMatrix.o \
	   TraitMatrix.o


HEADER = BranchRateManager.h \
         GSL.h \
         HistoryMultiParameter.h \
         HistorySingleParameter.h \
         Likelihood.h \
         ML_multi.h \
         ML_multi_DownhillSimplex.h \
         ML_multi_Powell.h \
         ML_multi_QuasiNewton.h \
         ML_single.h \
         ML_single_NewtonRaphson.h \
         PhyloTree.h \
         PhyloTreeNode.h \
         RatePMap.h \
         RatePVector.h \
         Simulate.h \
         TaxonMatrix.h \
         TraitMatrix.h

BIN  = branchrates


all: $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(CXXLDFLAGS) $(OBJ) -o $@ $(LIBS)

$(OBJ): $(HEADER)

clean:
	rm -f $(OBJ) $(BIN)

