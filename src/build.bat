@echo off

mkdir ..\build
pushd ..\build
cl /Zi ..\src\fmm_sim.cpp 
rem cl /O2i /INCREMENTAL:NO ..\src\fmm_sim.cpp
popd
