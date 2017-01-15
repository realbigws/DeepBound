

# make the subdir first
CXX = g++
FLAG = -O3
SHELL :=/bin/bash

all:
	$(CXX) distribute_ii_openmp.cpp $(FLAG) -fopenmp -o distribute_ii_openmp
	$(CXX) Fast_CAT.cpp $(FLAG) -fopenmp -o Fast_CAT
	$(CXX) Fast_CAT_v2.cpp $(FLAG) -fopenmp -o Fast_CAT_v2
	$(CXX) Generate_Feature.cpp $(FLAG) -fopenmp -o Generate_Feature
	$(CXX) Feature_Merge.cpp $(FLAG) -o Feature_Merge
	$(CXX) Feature_Merge_wind.cpp $(FLAG) -o Feature_Merge_wind
	$(CXX) Abundance_Estimate.cpp $(FLAG) -o Abundance_Estimate
	$(CXX) Decompose_Commented_File.cpp $(FLAG) -o Decompose_Commented_File
	mv distribute_ii_openmp Fast_CAT Fast_CAT_v2 Generate_Feature Feature_Merge Feature_Merge_wind Abundance_Estimate Decompose_Commented_File ../bin

# clean all the temp files
clean:
