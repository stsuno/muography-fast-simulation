# muography-fast-simulation

General description
===================
This software allows to simulate
  1. cosmic muon generation based on EXPACS [1],
  2. material construction in the front of muon detector and caluculate the enery-loss and scattering angle to total path length of the muon,
  3. modeling of various muon detector response and output hit-information in the file,
  4. reconstruction of the muon tracks based on the output file from the detector.

Code requires the ROOT package [2] provided by CERN [3], which is a free-software for academic use.
Within ROOT framework, it will be enable to see the graphical view of the simulation and to study the statistical analysis. 


How to run:
===========
There are two possible ways: 
  (1) compile executable by cmake
  (2) use root-prompt interpreter

(1) Compile executable by cmake

  This approach is relatively faster, however it will 
  be hard to see the images for the constructed geometries.
  It will be suitalble to run the large scale MC simulation.

  ~% cd work/
  ~% cmake ..
  ~% make

  The executables are located in the bin/ directory.
  Then, type

  ~% ./bin/SimulationWaterPipe()
  ~% ./bin/ReconstructionWaterPipe()

  The output files are created in the data/output/ directory.


(2) Use root-prompt interpreter

  This approach will be slower than the (1). But the constructed 
  geometries can be monitored once the code is run.
  It will be suitalble to develop the simulation.

  ~% cd work/
  ~% cp -R ../data/ ./

  Then, type

  ~% root

  In the root-prompt,

  root[0] gSystem->AddIncludePath("-I${PWD}/../")
  root[1] .L SimulationWaterPipe.cxx+

  Then the run simulation by typing,

  root[2] SimulationWaterPipe()

  In the same way,

  root[3] .L ReconstructionWaterPipe.cxx+
  root[3] ReconstructionWaterPipe()

  The output files are created in the data/output/ directory.



