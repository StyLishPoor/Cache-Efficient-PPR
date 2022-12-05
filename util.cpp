#include "util.h"

void IdRandom (const string& readfilename, const string& attributefilename, const string& writefilename)
{
  ifstream ifs(readfilename, ios::in | ios::binary);
  if(!ifs) {
    cout << "Cannot open: " << readfilename << endl;
    exit(1);
  }

  ifstream ifs_attribute(attributefilename);
  if (!ifs_attribute) {
    cout << "Cannot open: " << attributefilename << endl;
    exit(1);
  }
  long long n, m;
  ifs_attribute >> n >> m;

  vector<int> order(n);
  iota(order.begin(), order.end(), 0);

  random_device seed_gen;
  mt19937 engine(seed_gen());
  shuffle(order.begin(), order.end(), engine);

  vector<pair<int, int>> edge(m);
  vector<pair<int, int>> random_edge(m);
  ifs.read((char *)&edge[0], edge.size() * sizeof(edge[0]));
  for (long long i = 0; i < m; i++) {
    int from = order[edge[i].first];
    int to = order[edge[i].second];
    random_edge[i] = make_pair(from, to);
  }

  ifs.close();
  ifs_attribute.close();

  ofstream ofs(writefilename, ios::out | ios::binary);
  if (!ofs) {
    cout << "Cannot open: " << writefilename << endl;
    exit(1);
  }

  ofs.write((char *)&random_edge[0], random_edge.size() * sizeof(random_edge[0]));
  ofs.close();
  cout << "Random Binary wrote" << endl;

}

void ToDirected(const string& readfilename, const string& writefilename)
{
  ifstream ifs(readfilename, ios::in);
  if(!ifs) {
    cout << "Cannot open: " << readfilename << endl;
    exit(1);
  }

  vector<pair<int, int>> edge;
  edge.reserve(500000000);
  int from, to;
  while (ifs >> from >> to) {
    edge.push_back(make_pair(from, to));
    edge.push_back(make_pair(to, from));
  }
  ifs.close();

  ofstream ofs(writefilename, ios::out);
  if (!ofs) {
    cout << "Cannot open: " << writefilename << endl;
  }

  for (long i = 0; i < edge.size(); i++) {
    ofs << edge[i].first << " " << edge[i].second << endl;;
  }
  ofs.close();
}
