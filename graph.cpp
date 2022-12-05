#include "graph.h"
#include "random.h"

Graph::Graph()
{
  n = m = 0;
}

Graph::~Graph()
{
  delete[] outdegree;
  delete[] indegree;
  delete[] outoffset;
  delete[] inoffset;
  delete[] outedge;
  delete[] inedge;
  delete[] residue;
  delete[] reserve;
}

void Graph::clear()
{
  n = 0;
  m = 0;
  file_name.clear();
  delete[] outdegree;
  delete[] indegree;
  delete[] outoffset;
  delete[] inoffset;
  delete[] outedge;
  delete[] inedge;
  delete[] residue;
  delete[] reserve;
}

void Graph::reset()
{
  for (int i = 0; i < n; i++) {
    residue[i] = reserve[i] = 0.0;
  }
}

void Graph::graphLoad(const string& readfilename, const string& attributefilename)
{
  ifstream ifs(readfilename, ios::in);
  if(!ifs) {
    cout << "Cannot open: " << readfilename << endl;
    exit(1);
  }

  ifstream ifs_attribute(attributefilename, ios::in);
  if (!ifs_attribute) {
    cout << "Cannot open: " << attributefilename << endl;
    exit(1);
  }
  ifs_attribute >> n >> m;

  outdegree = new int[n];
  indegree = new int[n];
  residue = new double[n];
  reserve = new double[n];
  outoffset = new int[n+1];
  inoffset = new int[n+1];

  vector<pair<int, int>> edge;
  edge.reserve(500000000);
  int from, to;
  while (ifs >> from >> to) {
    outdegree[from]++;
    indegree[to]++;
    edge.push_back(make_pair(from, to));
  }
  m = edge.size();

  outedge = new int[m];
  inedge = new int[m];

  outoffset[0] = 0;
  inoffset[0] = 0;

  for (int i = 0; i < n; i++) {
    outoffset[i+1] = outoffset[i] + outdegree[i];
    inoffset[i+1] = inoffset[i] + indegree[i];
  }

  vector<int> inpos(n);
  sort(edge.begin(), edge.end());
  for (long long i = 0; i < m; i++) {
    int from = edge[i].first;
    int to = edge[i].second;
    outedge[i] = to;
    inedge[inoffset[to] + inpos[to]] = from;
    inpos[to]++;
  }
}

void Graph::ResidueRemain()
{
  double remain = 0.0;
  for (int i = 0; i < n; i++) {
    remain += residue[i];
  }
  cout << "Residue Remain: " << remain << endl;
}

void Graph::ReserveTotal()
{
  double total = 0.0;
  for (int i = 0; i < n; i++) {
    total += reserve[i];
  }
  cout << "Reserve Total: " << total << endl;
}

void Graph::PPRatio()
{
  int overzero = 0;
  int n1 = 0;
  for (int i = 0; i < n; i++) {
    if (reserve[i] > 0.0) {
      overzero++;
    }
    if (reserve[i] > 1.0/n) {
      n1++;
    }
  }
  cout << "Reserve over zero ratio: " << (double)100 * overzero / n << endl;
  cout << "Reserve 1/n ratio: " << (double)100 * n1 / n << endl;
}

void Graph::Accuracy(double delta, double epsilon, double pf)
{
  int pf_count = 0;
  int delta_count = 0;
  for (int i = 0; i < n; i++) {
    if (ground_truth[i] > delta) {
      delta_count++;
      if (abs(ground_truth[i] - reserve[i]) <= epsilon * ground_truth[i]) {
        pf_count++;
      }
    }
  }
  if ((double)pf_count/delta_count >= 1.0-pf) {
    cout << "Accurate" << endl;
  } else {
    cout << "Inaccurate" << endl;
  }
}

void Graph::queryGenerate(int query_size, const string& queryfilename)
{
  ofstream ofs(queryfilename);
  if (!ofs) {
    cout << "Cannot open: " << queryfilename << endl;
    exit(0);
  }

  int query_count = 0;
  Random r = Random();
  while (query_count < query_size) {
    int query_node = r.generateInt() % n;
    if (outDegree(query_node) > 0) {
      ofs << query_node << endl;
      query_count++;
    }
  }
  ofs.close();
}

void Graph::queryLoad(int query_size, vector<int>& query_nodes, const string& queryfilename)
{
  ifstream ifs(queryfilename);
  if (!ifs) {
    cout << "Cannot open: " << queryfilename << endl;
    exit(0);
  }
  int line_count = 0;
  int query_node;
  while(line_count < query_size && ifs >> query_node) {
    query_nodes[line_count++] = query_node;
  }
  if (line_count != query_size - 1) {
    query_nodes.resize(line_count);
  }
  ifs.close();
}

void Graph::MonteCarlo(int sourceNode, double alpha, double W)
{
  Random r = Random();
  int total_rw = ceil(W);
  for (int i = 0; i < total_rw; i++) {
    int currentNode = sourceNode;
    while (r.generateReal() > alpha) {
      int len = outDegree(currentNode);
      if (len == 0) {
        currentNode = sourceNode;
      } else {
        int neighbor = r.generateInt() % len;
        currentNode = outedge[outoffset[currentNode] + neighbor];
      }
    }
    reserve[currentNode] += 1.0 / total_rw;
  }
}

void Graph::PowerMethod(int sourceNode, double alpha, int max_iter)
{
  for (int i = 0; i < n; i++) {
    residue[i] = 0.0;
  }
  residue[sourceNode] = 1.0;
  ground_truth.resize(n);
  fill(ground_truth.begin(), ground_truth.end(), 0);
  double finish_condition = (1/m < 10e-10) ? 1/m : 10e-10;
  double rsum = 1.0;
  for (int iter = 0; iter < max_iter; iter++) {
    for (int i = 0; i < n; i++) {
      double current_residue = residue[i];
      if (current_residue > 0.0) {
        double remain_residue = (1 - alpha) * current_residue;
        int out_deg = outDegree(i);
        ground_truth[i] += alpha * current_residue;
        residue[i] = 0.0;
        if (out_deg == 0) {
          residue[sourceNode] += remain_residue;
          continue;
        }
        double avg_push_residue = remain_residue / out_deg;
        for (int start = outoffset[i]; start < outoffset[i+1]; start++) {
          int next = outedge[start];
          residue[next] += avg_push_residue;
        }
      }
    }
  }
}

void Graph::ForwardPush(int sourceNode, double alpha, double rmax)
{
  static vector<bool> idx(n);
  fill(idx.begin(), idx.end(), false);
  idx[sourceNode] = true;

  residue[sourceNode] = 1.0;

  queue<int> active;
  active.push(sourceNode);
  while(active.size() > 0) {
    int tempNode = active.front();
    active.pop();
    idx[tempNode] = false;
    double current_residue = residue[tempNode];
    if(current_residue > 0.0){
      residue[tempNode] = 0.0;
      reserve[tempNode] += alpha * current_residue;
      double remain_residue = (1 - alpha) * current_residue;

      int out_deg = outDegree(tempNode);
      if (out_deg == 0) {
        residue[sourceNode] += remain_residue;
        if ((residue[sourceNode] / outDegree(sourceNode) >= rmax) && idx[sourceNode]!= true) {
          idx[sourceNode] = true;
          active.push(sourceNode);
        }
        continue;
      }
      double avg_push_residue = remain_residue / out_deg;
      for (int start = outoffset[tempNode]; start < outoffset[tempNode+1]; start++) {
        int next = outedge[start];
        residue[next] += avg_push_residue;
        if ((residue[next] / outDegree(next) >= rmax) && idx[next] != true) {
          idx[next] = true;
          active.push(next);
        }
      }
    }
  }
}

void Graph::Fora(int sourceNode, double alpha, double W, double rmax)
{
  static vector<bool> idx(n);
  std::fill(idx.begin(), idx.end(), false);
  idx[sourceNode] = true;

  residue[sourceNode] = 1.0;

  queue<int> active;
  active.push(sourceNode);
  while(active.size() > 0) {
    int tempNode = active.front();
    active.pop();
    idx[tempNode] = false;
    double current_residue = residue[tempNode];
    if(current_residue > 0.0){
      residue[tempNode] = 0.0;
      reserve[tempNode] += alpha * current_residue;
      double remain_residue = (1 - alpha) * current_residue;

      int out_deg = outDegree(tempNode);
      if (out_deg == 0) {
        residue[sourceNode] += remain_residue;
        if ((residue[sourceNode] / outDegree(sourceNode) >= rmax) && idx[sourceNode]!= true) {
          idx[sourceNode] = true;
          active.push(sourceNode);
        }
        continue;
      }
      double avg_push_residue = remain_residue / out_deg;
      for (int start = outoffset[tempNode]; start < outoffset[tempNode+1]; start++) {
        int next = outedge[start];
        residue[next] += avg_push_residue;
        if ((residue[next] / outDegree(next) >= rmax) && idx[next] != true) {
          idx[next] = true;
          active.push(next);
        }
      }
    }
  }
  RW(sourceNode, alpha, W);
}

void Graph::CachePPR(int sourceNode, double alpha, double W, double rmax)
{
  // BFS-based reordering
  int* neworder;
  int* c_queue;
  int* newoutoffset;
  int* newoutedges;

  vector<bool> found(n);
  c_queue = new int[n];
  neworder = new int[n];
  newoutoffset = new int[n+1];
  newoutedges = new int[m];

  found[sourceNode] = true;

  neworder[sourceNode]=0;
  newoutoffset[0]=0;

  vector<int> distance_index;
  distance_index.push_back(0);
  int left=0, right=1;
  int sentinel = 0;
  int current_hop = 0;
  c_queue[left] = sourceNode;
  while (right > left) {
    if (sentinel == left) {
      sentinel = right;
      distance_index.push_back(right);
    }
    int temp = c_queue[left];
    left++;
    int outdeg = outDegree(temp);
    newoutoffset[neworder[temp]+1]=newoutoffset[neworder[temp]] + outdeg;
    for (int start = outoffset[temp], i=0; start < outoffset[temp+1]; start++, i++) {
      int next = outedge[start];
      if (found[next] != true) {
        c_queue[right] = next;
        found[next] = true;
        neworder[next] = right++;
      }
      newoutedges[newoutoffset[neworder[temp]]+i]=neworder[next];
    }
  }
  distance_index.push_back(right);

  delete[] c_queue;
  delete[] neworder;

  int* swap_outoffset=outoffset;
  int* swap_outedge=outedge;
  int swap_n = n;
  int swap_sourceNode = sourceNode;
  outoffset=newoutoffset;
  outedge=newoutedges;
  n=right;
  sourceNode=0;

  // Hop-Extension Forward Push
  residue[sourceNode] = 1.0;

  int loop_end;
  int iter = 0;
  bool unit_flag;
  do {
    unit_flag = true;
    if (iter < distance_index.size()-1) {
      loop_end = distance_index[++iter];
    } else {
      loop_end = n;
    }
    for (int i = 0; i < loop_end; i++) {
      double current_residue = residue[i];
      if (current_residue > outDegree(i) * rmax) {
        unit_flag = false;
        double remain_residue = (1 - alpha) * current_residue;
        int out_deg = outDegree(i);
        if (out_deg == 0) {
          residue[i] = 0.0;
          reserve[i] += alpha * current_residue;
          residue[sourceNode] += remain_residue;
          continue;
        }
        residue[i] = 0;
        reserve[i] += alpha * current_residue;
        double avg_push_residue = remain_residue / out_deg;
        for (int start = outoffset[i]; start < outoffset[i+1]; start++) {
          int next = outedge[start];
          residue[next] += avg_push_residue;
        }
      }
    }
  } while (!unit_flag);

  // Vertex-Centric Random Walk
  VCRW(sourceNode, alpha, W);

  // Remedy
  outoffset = swap_outoffset;
  outedge = swap_outedge;
  n = swap_n;
  sourceNode = swap_sourceNode;
  delete[] newoutoffset;
  delete[] newoutedges;
}

void Graph::ResAcc(int sourceNode, int h, double alpha, double W, double rmax_hop, double rmax_f)
{
  unordered_set<int> khopset;  
  unordered_set<int> k1hopset;  
  khopset.insert(sourceNode);

  int* hops_from_source;
  hops_from_source = new int[n];
  for (int i = 0; i < n; i++) {
    hops_from_source[i] = numeric_limits<int>::max();
  }
  hops_from_source[sourceNode] = 0;

  static vector<bool> idx(n);
  std::fill(idx.begin(), idx.end(), false);

  queue<int> active;

  // sourceNode
  int out_deg = outDegree(sourceNode);
  reserve[sourceNode] = alpha;
  double remain_residue = 1 - alpha;
  if (out_deg == 0) {
    residue[sourceNode] += remain_residue;
  } else {
    double avg_push_residue = remain_residue / out_deg;
    for (int start = outoffset[sourceNode]; start < outoffset[sourceNode+1]; start++) {
      int next = outedge[start];
      hops_from_source[next] = 1;
      khopset.insert(next);
      residue[next] += avg_push_residue;
      if (h > 0 && next != sourceNode && outDegree(next) != 0 && residue[next] / outDegree(next) >= rmax_hop && idx[next] != true) {
        active.push(next);
        idx[next] = true;
      } else if (h == 0) {
        k1hopset.insert(next);
      }
    }
  }

  // kHopFWD Phase
  while(active.size() > 0) {
    int tempNode = active.front();
    active.pop();
    idx[tempNode] = false;
    double current_residue = residue[tempNode];
    if(current_residue > 0.0){
      residue[tempNode] = 0.0;
      reserve[tempNode] += alpha * current_residue;
      double remain_residue = (1 - alpha) * current_residue;

      int out_deg = outDegree(tempNode);
      if (out_deg == 0) {
        residue[sourceNode] += remain_residue;
        continue;
      }

      double avg_push_residue = remain_residue / out_deg;
      int temp_hops = hops_from_source[tempNode];
      for (int start = outoffset[tempNode]; start < outoffset[tempNode+1]; start++) {
        int next = outedge[start];
        residue[next] += avg_push_residue;
        if (next != sourceNode) {
          if (hops_from_source[next] > temp_hops + 1) {
            hops_from_source[next] = temp_hops + 1;
          }
          if (hops_from_source[next] <= h) {
            if (outDegree(next) != 0 && residue[next] / outDegree(next) >= rmax_hop) {
              if (idx[next] != true) {
                active.push(next);
                idx[next] = true;
              }
            }
            if (khopset.count(next) == 0) {
              khopset.insert(next);
              if (k1hopset.count(next) > 0) {
                k1hopset.erase(next);
              }
            }
          } else if (hops_from_source[next] == h + 1) {
            if (k1hopset.count(next) == 0) {
              k1hopset.insert(next);
            }
          }
        } else {
          continue;
        }
      }
    }
  }
  // loop at sourceNode
  for (int i = 0; i < n; i++) {
    if (outDegree(i) == 0 && residue[i] != 0.0) {
      reserve[i] += alpha * residue[i];
      residue[sourceNode] += (1 - alpha) * residue[i];
      residue[i] = 0.0;
    }
  }

  if (h > 0 && residue[sourceNode] != 0) {
    double min_rmax = 10e-14;
    unsigned long source_loop = (int) ceil(log(min_rmax * outDegree(sourceNode)) / log(residue[sourceNode]));
    double scale = (1 - pow(residue[sourceNode], source_loop - 1)) / (1 - residue[sourceNode]);
    for (auto iter = khopset.begin(); iter != khopset.end(); iter++) {
      int tempNode = *iter;
      reserve[tempNode] = reserve[tempNode] * scale;
      if (tempNode != sourceNode) {
        residue[tempNode] = residue[tempNode] * scale;
      } else {
        residue[tempNode] = pow(residue[tempNode], source_loop);
      }
    }
    for (auto iter = k1hopset.begin(); iter != k1hopset.end(); iter++) {
      int tempNode = *iter;
      residue[tempNode] = residue[tempNode] * scale;
    }
  }

  // OMFWD Phase
  fill(idx.begin(), idx.end(), false);
  queue<int>().swap(active);
  for (auto iter = k1hopset.begin(); iter != k1hopset.end(); iter++) {
    int k1node = *iter;
    active.push(k1node);
    idx[k1node] = true;
  }

  while(active.size() > 0) {
    int tempNode = active.front();
    active.pop();
    idx[tempNode] = false;
    double current_residue = residue[tempNode];
    if(current_residue > 0.0){
      residue[tempNode] = 0.0;
      reserve[tempNode] += alpha * current_residue;
      double remain_residue = (1 - alpha) * current_residue;

      int out_deg = outDegree(tempNode);
      if (out_deg == 0) {
        residue[sourceNode] += remain_residue;
        if ((residue[sourceNode] / outDegree(sourceNode) >= rmax_f) && idx[sourceNode]!= true) {
          idx[sourceNode] = true;
          active.push(sourceNode);
        }
        continue;
      }
      double avg_push_residue = remain_residue / out_deg;
      for (int start = outoffset[tempNode]; start < outoffset[tempNode+1]; start++) {
        int next = outedge[start];
        residue[next] += avg_push_residue;
        if ((residue[next] / outDegree(next) >= rmax_f) && idx[next] != true) {
          idx[next] = true;
          active.push(next);
        }
      }
    }
  }
  RW(sourceNode, alpha, W);
}

void Graph::SpeedPPR(int sourceNode, double alpha, double W, double rmax)
{
  static vector<bool> idx(n);
  std::fill(idx.begin(), idx.end(), false);
  idx[sourceNode] = true;

  residue[sourceNode] = 1.0;

  queue<int> active;
  active.push(sourceNode);

  double lambda = rmax * m;
  double scanThr = n/4;
  double rsum = 1.0;
  while(active.size() > 0 && active.size() < scanThr && rsum > lambda) {
    int tempNode = active.front();
    active.pop();
    idx[tempNode] = false;
    double current_residue = residue[tempNode];
    if(current_residue > 0.0){
      double remain_residue = (1 - alpha) * current_residue;

      int out_deg = outDegree(tempNode);
      residue[tempNode] = 0.0;
      reserve[tempNode] += alpha * current_residue;
      rsum -= alpha * current_residue;
      if (out_deg == 0) {
        residue[sourceNode] += remain_residue;
        if ((residue[sourceNode] / outDegree(sourceNode) >= rmax) && idx[sourceNode]!= true) {
          idx[sourceNode] = true;
          active.push(sourceNode);
        }
        continue;
      }
      double avg_push_residue = remain_residue / out_deg;
      for (int start = outoffset[tempNode]; start < outoffset[tempNode+1]; start++) {
        int next = outedge[start];
        residue[next] += avg_push_residue;
        if ((residue[next] / outDegree(next) >= rmax) && idx[next] != true) {
          idx[next] = true;
          active.push(next);
        }
      }
    }
  }

  if (rsum > lambda) {
    static int epoch_num = 8;

    for (int e = 1; e <= epoch_num; e++) {
      double rmax_epoch = pow(lambda, (double)e/epoch_num)/m;
      while(rsum > m * rmax_epoch) {
        for (int i = 0; i < n; i++) {
          double current_residue = residue[i];
          if (current_residue > 0.0) {
            double remain_residue = (1 - alpha) * current_residue;
            int out_deg = outDegree(i);
            if (out_deg == 0) {
              residue[i] = 0.0;
              reserve[i] += alpha * current_residue;
              rsum -= alpha * current_residue;
              residue[sourceNode] += remain_residue;
              continue;
            }
            if (current_residue / out_deg >= rmax_epoch) {
              residue[i] = 0;
              reserve[i] += alpha * current_residue;
              rsum -= alpha * current_residue;
              double avg_push_residue = remain_residue / out_deg;
              for (int start = outoffset[i]; start < outoffset[i+1]; start++) {
                int next = outedge[start];
                residue[next] += avg_push_residue;
              }
            }
          }
        }
      }
    }
  }
  queue<int>().swap(active);
  std::fill(idx.begin(), idx.end(), false);
  for (int i = 0; i < n; i++) {
    if (residue[i] >= rmax * outDegree(i)) {
      idx[i] = true;
      active.push(i);
    }
  }
  while(active.size() > 0) {
    int tempNode = active.front();
    active.pop();
    idx[tempNode] = false;
    double current_residue = residue[tempNode];
    if(current_residue > 0.0){
      double remain_residue = (1 - alpha) * current_residue;

      int out_deg = outDegree(tempNode);
      residue[tempNode] = 0.0;
      reserve[tempNode] += alpha * current_residue;
      rsum -= alpha * current_residue;
      if (out_deg == 0) {
        residue[sourceNode] += remain_residue;
        if ((residue[sourceNode] / outDegree(sourceNode) >= rmax) && idx[sourceNode]!= true) {
          idx[sourceNode] = true;
          active.push(sourceNode);
        }
        continue;
      }
      double avg_push_residue = remain_residue / out_deg;
      for (int start = outoffset[tempNode]; start < outoffset[tempNode+1]; start++) {
        int next = outedge[start];
        residue[next] += avg_push_residue;
        if ((residue[next] / outDegree(next) >= rmax) && idx[next] != true) {
          idx[next] = true;
          active.push(next);
        }
      }
    }
  }
  RW(sourceNode, alpha, W);
}

void Graph::RW(int sourceNode, double alpha, double W)
{
  Random r = Random();
  for (int i = 0; i < n; i++) {
    if (residue[i]  > 0.0) {
      int tempNode = i;
      double current_residue = residue[tempNode];
      int walknum = ceil(current_residue * W);
      for (int j = 0; j < walknum; j++) {
        int currentNode = tempNode;
        while (r.generateReal() > alpha) {
          int len = outDegree(currentNode);
          if (len == 0) {
            currentNode = sourceNode;
          } else {
            int neighbor = r.generateInt() % len;
            currentNode = outedge[outoffset[currentNode] + neighbor];
          }
        }
        reserve[currentNode] += current_residue / walknum;
      }
    }
  }
}

void Graph::VCRW(int sourceNode, double alpha, double W)
{
  Random r = Random();
  long global_total = 0;
  int* globalRW = new int[n];
  bool* localRW = new bool[n];
  for (int i = 0; i < n; i++) {
    localRW[i] = false;
    globalRW[i] = 0;
    if (residue[i] > 0.0) {
      localRW[i] = true;
      int Wi = floor(residue[i] * W);
      if (Wi >= 1) {
        global_total += Wi;
        globalRW[i] = Wi;
        residue[i] -= Wi/W;
      }
    }
  }

  // Independent Random Walk
  while (global_total > 0) {
    for (int i = 0; i < n; i++) {
      if (globalRW[i] > 0) {
        for (int j = 0; j < globalRW[i]; j++) {
          if (r.generateReal() > alpha) {
            int len = outDegree(i);
            if (len == 0) {
              globalRW[sourceNode]++;
            } else {
              int neighbor = r.generateInt() % len;
              int nextnode = outedge[outoffset[i] + neighbor];
              globalRW[nextnode]++;
            }
          } else {
            global_total--;
            reserve[i] += 1/W;
          }
        }
        globalRW[i] = 0;
      }
    }
  }

  // Dependent Random Walk
  for (int i = 0; i < n; i++) {
    if (localRW[i] == true) {
      int currentNode = i;
      double current_residue = residue[currentNode];
      while (r.generateReal() > alpha) {
        int len = outDegree(currentNode);
        if (len == 0) {
          currentNode = sourceNode;
        } else {
          int neighbor = r.generateInt() % len;
          currentNode = outedge[outoffset[currentNode] + neighbor];
        }
      }
      reserve[currentNode] += current_residue;
    }
  }
  delete[] globalRW;
  delete[] localRW;
}
