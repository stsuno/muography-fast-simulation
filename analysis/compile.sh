g++ -c amain.cxx -I$ROOTSYS/include
g++ -c KomenAnalysis.cxx -I$ROOTSYS/include -I$PWD
g++ -c SakurajimaAnalysis.cxx -I$ROOTSYS/include -I$PWD
g++ -c DemoTest.cxx -I$ROOTSYS/include -I$PWD
g++ -c track_tree.C -I$ROOTSYS/include -I/$PWD
g++ -o makeHist amain.o KomenAnalysis.o SakurajimaAnalysis.o DemoTest.o track_tree.o `root-config --cflags --libs`

