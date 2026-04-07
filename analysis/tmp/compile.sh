g++ -c amain.cxx -I$ROOTSYS/include
g++ -c ZmumuAnalysis.cxx -I$ROOTSYS/include -I/home/tsuno/work/pixelOffline.01.19.2022/analysis/makeHist/
g++ -c makeTree.cxx -I$ROOTSYS/include -I/home/tsuno/work/pixelOffline.01.19.2022/analysis/makeHist/
g++ -c makeEventTree2.cxx -I$ROOTSYS/include -I/home/tsuno/work/pixelOffline.01.19.2022/analysis/makeHist/
g++ -c MyClassDict.cxx -I$ROOTSYS/include -I/home/tsuno/work/pixelOffline.01.19.2022/analysis/makeHist/
g++ -c eventtree.C -I$ROOTSYS/include -I/home/tsuno/work/pixelOffline.01.19.2022/analysis/makeHist/
g++ -o makeHist amain.o ZmumuAnalysis.o makeTree.o MyClassDict.o makeEventTree2.o eventtree.o -m64 -fPIC -pthread -I$ROOTSYS/include -L$ROOTSYS/lib -lCore -lRIO -lNet -lThread -lMathCore -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -pthread -lm -ldl -rdynamic

