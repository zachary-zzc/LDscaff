#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <stdlib.h>
#include <cstdlib>
#include <assert.h>
#include "tools.h"
#include "cal.h"

#include <lemon/matching.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/lgf_reader.h>
#include <lemon/math.h>

#define INF 65535
#define NSAMPLE 1000

using namespace std;
using namespace lemon;

typedef lemon::MaxWeightedMatching<lemon::SmartGraph, lemon::SmartGraph::EdgeMap<float> > MWM;

struct scafGraph {
    string *vexs;
    float **edge;
    int nvex, nedge;
    
    scafGraph(int n);
    ~scafGraph();
    void printGraph();
    string to_lgf();
};

// Core functions moved from graph.cc
scafGraph *createGraph(char** infList, int nsample, int nscaf, int nmarker, ofstream &logf);
int processMatching(char* sampleListFile, int nsample, int nscaf, int nmarker, const char* logFile);

#endif
