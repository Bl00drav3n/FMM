@echo off

mkdir ..\build
pushd ..\build
cl -Zi /INCREMENTAL:NO ..\src\fmm_sim.cpp 
popd
