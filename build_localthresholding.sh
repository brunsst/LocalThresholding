#!/bin/bash

############# Parameters ################
#########################################

#########################################
#########################################

mkdir -p "Build"

echo "building LocalThresholding"
echo "--------------------------------------------------"
echo "compiling auxiliary.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/auxiliary.cpp" -o $PWD"/Build/auxiliary.o"
echo "compiling hdcommunication.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/hdcommunication.cpp" -o $PWD"/Build/hdcommunication.o"
echo "compiling histogram.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/histogram.cpp" -o $PWD"/Build/histogram.o"
echo "compiling main.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/main.cpp" -o $PWD"/Build/main.o"
echo "linking LocalThresholding"
g++  -o $PWD/locthresh $PWD/Build/main.o  $PWD/Build/auxiliary.o $PWD/Build/hdcommunication.o $PWD/Build/histogram.o  -ltiff -lgomp
echo "--------------------------------------------------"

################################################################################################################################################################
################################################################################################################################################################
