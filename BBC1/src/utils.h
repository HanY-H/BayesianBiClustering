#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <iostream>
#include <vector>
#include <string>

using namespace std;


template<class T>
class matrix : public vector<vector<T> > {
 public:
  matrix() {}
  matrix(int dim1, int dim2, T defVal = T())
    : vector<vector<T> >(dim1, vector<T>(dim2, defVal)) {
  }
};

namespace ctnutils {

template<class T>
void DispVector(const vector<T> & vec, string sep = "\t", bool line_break = true) {
  if (vec.size() == 0) {
    cout << "(empty vector)";
    return;
  }
  cout << "(";
  for (unsigned i = 0; i < vec.size() - 1; i++)
    cout << vec[i] << sep;
  cout << vec[vec.size() - 1] << ")";
  if (line_break) cout << endl;
}

}

namespace strutils {

inline string trim(string str, string ts = " \t\n") {
  str.erase(0, str.find_first_not_of(ts));       //prefixing spaces
  str.erase(str.find_last_not_of(ts) + 1);       //surfixing spaces

  return str;
}

vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string>  split(const string &s, char delim);

}

namespace numutils {

string EnoughSpaces(int MAX, int num);

}

#endif // UTILS_H_INCLUDED
