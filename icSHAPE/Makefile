CXX        = g++

CXXFLAGS    = -O3 -std=c++0x -Wall

all: calcEnrich calcRT combineRTreplicates filterEnrich normalizeRTfile normedRT2bedGraph


calcEnrich: calcEnrich.cpp
	$(CXX) calcEnrich.cpp $(CXXFLAGS) -lPsBL -o calcEnrich 

calcRT: calcRT.cpp
	$(CXX) calcRT.cpp $(CXXFLAGS) -lPsBL -o calcRT 

combineRTreplicates: combineRTreplicates.cpp
	$(CXX) combineRTreplicates.cpp $(CXXFLAGS) -lPsBL -o combineRTreplicates 

filterEnrich: filterEnrich.cpp
	$(CXX) filterEnrich.cpp $(CXXFLAGS) -lPsBL -o filterEnrich

normalizeRTfile: normalizeRTfile.cpp
	$(CXX) normalizeRTfile.cpp $(CXXFLAGS) -lPsBL -o normalizeRTfile

normedRT2bedGraph: normedRT2bedGraph.cpp
	$(CXX) normedRT2bedGraph.cpp $(CXXFLAGS) -lPsBL -o normedRT2bedGraph


clean:
	rm calcEnrich || true
	rm calcRT || true
	rm combineRTreplicates || true
	rm filterEnrich || true
	rm normalizeRTfile || true
	rm normedRT2bedGraph || true


