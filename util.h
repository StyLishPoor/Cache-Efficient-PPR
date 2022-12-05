#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <random>

using namespace std;

void IdRandom (const string& readfilename, const string& attributefilename, const string& writefilename);
void ToDirected(const string& readfilename, const string& writefilename);

#endif
