#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <set>
#include <functional>
#include <climits>
#include <ctime>
#include <stdlib.h>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <chrono>

#include "random.h"
#include "graph.h"
#include "util.h"
#include "helper.h"

using namespace std;

int main(int argc, char* argv[]){
  ios::sync_with_stdio(false);
  param = argParse(argc, argv);

  if (param.algorithm == "ToDirected") {
    ToDirected(param.graph_name, param.graph_name);
    return 0;
  }

  Graph g;
  clock_t start, end;
  if (!param.graph_name.empty()) {
    start=clock();
    g.graphLoad(param.graph_name, param.attribute_name);
    end=clock();
    cout << "Load Time: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
    cout << "--------------------------------------" << endl;
  } else {
    cout << "Error" << endl;
    return 0;
  }

  if (param.algorithm == "QueryGenerate") {
    if (param.query_name.empty()) {
      cout << "Input Query File Name." << endl;
    } else {
      g.queryGenerate(param.query_size, param.query_name);
    }
    return 0;
  }
  
  // Personalized PageRank Computation
  vector<int> query_nodes(param.query_size);
  g.queryLoad(param.query_size, query_nodes, param.query_name);
  double W = (2 * (2 * param.epsilon/3 + 2) * log(g.n) * g.n) / (param.epsilon * param.epsilon);
  double rmax = 1/sqrt(g.m * W);

  if (param.algorithm == "ForwardPush") {
    double total_time = 0;
    for (int sourceNode : query_nodes) {
      start = clock();
      g.ForwardPush(sourceNode, param.alpha, rmax); 
      end = clock();
      total_time += (double)(end-start)/CLOCKS_PER_SEC;
      cout << "      Computation Time: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
      cout << endl;
      g.reset();
    }
    cout << "Average Time: " << total_time / query_nodes.size() << endl;
    return 0;
  } 

  if (param.algorithm == "MonteCarlo") {
    double total_time = 0;
    for (int sourceNode : query_nodes) {
      start = clock();
      g.MonteCarlo(sourceNode, param.alpha, W); 
      end = clock();
      total_time += (double)(end-start)/CLOCKS_PER_SEC;
      cout << "      Computation Time: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
      cout << endl;
      g.reset();
    }
    cout << "Average Time: " << total_time / query_nodes.size() << endl;
    return 0;
  } 

  if (param.algorithm == "Fora") {
    double total_time = 0;
    for (int sourceNode : query_nodes) {
      start = clock();
      g.Fora(sourceNode, param.alpha, W, rmax); 
      end = clock();
      total_time += (double)(end-start)/CLOCKS_PER_SEC;
      cout << "      Computation Time: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
      cout << endl;
      g.reset();
    }
    cout << "Average Time: " << total_time / query_nodes.size() << endl;
    return 0;
  } 

  if (param.algorithm == "ResAcc") {
    double total_time = 0;
    for (int sourceNode : query_nodes) {
      start = clock();
      g.ResAcc(sourceNode, 2, param.alpha, W, 10e-14, rmax); 
      end = clock();
      total_time += (double)(end-start)/CLOCKS_PER_SEC;
      cout << "      Computation Time: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
      cout << endl;
      g.reset();
    }
    cout << "Average Time: " << total_time / query_nodes.size() << endl;
    return 0;
  }
  
  if (param.algorithm == "SpeedPPR") {
    rmax = 1/W;
    double total_time = 0;
    for (int sourceNode : query_nodes) {
      start = clock();
      g.SpeedPPR(sourceNode, param.alpha, W, rmax); 
      end = clock();
      total_time += (double)(end-start)/CLOCKS_PER_SEC;
      cout << "      Computation Time: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
      cout << endl;
      g.reset();
    }
    cout << "Average Time: " << total_time / query_nodes.size() << endl;
    return 0;
  }

  if (param.algorithm == "CachePPR") {
    double total_time = 0;
    for (int sourceNode : query_nodes) {
      start = clock();
      g.CachePPR(sourceNode, param.alpha, W, rmax); 
      end = clock();
      total_time += (double)(end-start)/CLOCKS_PER_SEC;
      cout << "      Computation Time: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
      cout << endl;
      g.reset();
    }
    cout << "Average Time: " << total_time / query_nodes.size() << endl;
    return 0;
  }

  return 0;
}
