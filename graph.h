#ifndef _GRAPH_H
#define _GRAPH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <unordered_set>
#include <list>
#include <queue>
#include <algorithm>
#include <functional>
#include <utility>
#include <cmath>
#include <limits>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <string>
#include <cstdint>
#include <chrono>
#include <random>

using namespace std;

class Graph{
  public:
    int n;
    long long m;
    string file_name;

    int* outdegree;
    int* indegree;
    int* outoffset;
    int* inoffset;
    int* outedge;
    int* inedge;

    double* residue;
    double* reserve;
    vector<double> ground_truth;

    string getFilename();
    void setFilename(string name);

    Graph();
    ~Graph();
    void clear();
    void reset();
    void graphLoad(const string& readfilename, const string& attributefilename);
    void ResidueRemain();
    void ReserveTotal();
    void PPRatio();
    void Accuracy(double delta, double epsilon, double pf);

    inline int inDegree(int v)
    {
      return inoffset[v+1] - inoffset[v];
    }
    inline int outDegree(int v)
    {
      return outoffset[v+1] - outoffset[v];
    }
    void queryGenerate(int query_size, const string& queryfilename);
    void queryLoad(int query_size, vector<int>& query_nodes, const string& queryfilename);
    void MonteCarlo(int sourceNode, double alpha, double W);
    void PowerMethod(int sourceNode, double alpha, int max_iter);
    void ForwardPush(int sourceNode, double alpha, double rmax);
    void Fora(int sourceNode, double alpha, double W, double rmax);
    void CachePPR(int sourceNode, double alpha, double W, double rmax);
    void ResAcc(int sourceNode, int h, double alpha, double W, double rmax_hop, double rmax_f);
    void SpeedPPR(int sourceNode, double alpha, double W, double rmax);

    void RW(int sourceNode, double alpha, double W);
    void VCRW(int sourceNode, double alpha, double W);
};

#endif
