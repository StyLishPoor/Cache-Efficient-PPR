#include <unordered_set>
#include "helper.h"
#include "random.h"

Param param;
Param argParse(int argnum, char** argname)
{
  vector<string> Algorithms {
    "MonteCarlo", "ForwardPush", "Fora", "ResAcc", "SpeedPPR", "CachePPR", "ToDirected", "QueryGenerate"
  };
  Param rtn;
  unordered_set<string> algos(Algorithms.begin(), Algorithms.end());
  int count;
  for (count = 1; count < argnum;) {
    string para = argname[count++];
    if (count == argnum) {
      cout << "Cannot Parse." << endl;
      usage_message();
      exit(0);
    }
    if (para[0] != '-') {
      cout << "Unknown Parameter." << endl;
      usage_message();
      exit(0);
    }
    string arg = argname[count++];
    if (para == "-algo") {
      if (algos.find(arg) == algos.end()) {
        cout << "Unknown Algorithm." << endl;
        usage_message();
        exit(0);
      }
      rtn.algorithm = arg;
      cout << "Algorithm: " << rtn.algorithm << endl;
    } else if (para == "-graph") {
      rtn.graph_name = arg;
      cout << "Input Graph: " << rtn.graph_name << endl;
    } else if (para == "-attribute") {
      rtn.attribute_name = arg;
      cout << "Attribute File: " << rtn.attribute_name << endl;
    } else if (para == "-query") {
      rtn.query_name = arg;
      cout << "Query File: " << rtn.query_name << endl;
    } else if (para == "-query_size") {
      regex size_check(R"([1-9]\d{0,})");
      if (regex_match(arg, size_check)) {
        rtn.query_size = stoi(arg);
        cout << "Query Size: " << rtn.query_size << endl;
      } else {
        cout << "Unknown Query Size." << endl;
        usage_message();
        exit(0);
      }
    } else if (para == "-alpha") {
      regex alpha_check(R"(0.\d{1,})");
      if (regex_match(arg, alpha_check)) {
        rtn.alpha = stof(arg);
        cout << "Alpha: " << rtn.alpha << endl;
      } else {
        cout << "Unknown Alpha." << endl;
        usage_message();
        exit(0);
      }
    } else if (para == "-epsilon") {
      regex epsilon_check(R"(0.\d{1,})");
      if (regex_match(arg, epsilon_check)) {
        rtn.epsilon = stof(arg);
        cout << "Epsilon: " << rtn.epsilon << endl;
      } else {
        cout << "Unknown Epsilon." << endl;
        usage_message();
        exit(0);
      }
    }
  }
  return rtn;
}

void usage_message() {
  cout << "Usage: ./SSPPR -algo <algorithm> -graph <graph-path> -attribute <attribute-path> -query <query-path> -quer_sizey <must be greater than 1> [-alpha <must be between 0 and 1>] [-epsilon <must be between 0 and 1>]" << endl;
}
