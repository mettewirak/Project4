# Project4

Coursework in FYS4150 - Computational Physics at the University of Oslo.

We've here implemented the Metropolis algorithm (a Monte Carlo method) to model the phase transition of a ferromagnetic material. The Ising model is used to find the energy of each spin. 

All of our simulations are run in two dimensions, and we've used between 2 and 100 particles in each dimension. To be able to run the larger ones of these with a sufficiently high number of Monte Carlo cycles we've used MPI parallelisation. However, this is only implemented in the master-branch. We also used one computer without MPI, and thus have one version with all the MPI-parts of the code emitted. This can be found in the branch "mette". 



### How to run the code

If running with MPI on the school computers, this should be included in the .pro-file in QT:

TEMPLATE = app 

CONFIG += console c++11 

CONFIG -= app_bundle 

CONFIG -= qt 


SOURCES += main.cpp 


QMAKE_CXX = /usr/lib64/openmpi/bin/mpicxx 

QMAKE_CXX_RELEASE = $$QMAKE_CXX 

QMAKE_CXX_DEBUG = $$QMAKE_CXX 

QMAKE_LINK = $$QMAKE_CXX 

QMAKE_CC = /usr/lib64/openmpi/bin/mpicc 


QMAKE_CFLAGS += $$system(/usr/lib64/openmpi/bin/mpicc --showme:compile) 

QMAKE_LFLAGS += $$system(/usr/lib64/openmpi/bin/mpicxx --showme:link) 

QMAKE_CXXFLAGS += $$system(/usr/lib64/openmpi/bin/mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK 

QMAKE_CXXFLAGS_RELEASE += $$system(/usr/lib64/openmpi/bin/mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK -O3 

QMAKE_CXXFLAGS_DEBUG += $$system(/usr/lib64/openmpi/bin/mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK -O3 



One also has to change the run settings to the following: 

Executable: /usr/lib64/openmpi/bin/mpirun 

Command line arguments: -np 8 Project4 


Running the code on other computers with MPI, simply change the paths in both the .pro-file and the executable.
