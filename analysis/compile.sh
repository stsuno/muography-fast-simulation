g++ -c amain.cxx -I$ROOTSYS/include
g++ -c KomenAnalysis.cxx -I$ROOTSYS/include -I$PWD
g++ -c track_tree.C -I$ROOTSYS/include -I/$PWD
g++ -o makeHist amain.o KomenAnalysis.o track_tree.o -m64 -fPIC -pthread -I$ROOTSYS/include -L$ROOTSYS/lib -lCore -lRIO -lNet -lThread -lMathCore -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -pthread -lm -ldl -rdynamic

