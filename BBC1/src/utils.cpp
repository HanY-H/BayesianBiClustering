#include "utils.h"

#include <sys/types.h>
#include <dirent.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>

using namespace std;


namespace numutils {

int CountDigits(int num) {
  if (num == 0)
    return 1;
  int num_digits = (num < 0) ? 1 : 0;
  while(num != 0) {
    num_digits ++;
    num /= 10;
  }
  return num_digits;
}

string EnoughSpaces(int MAX, int num) {
  string spaces = "";
  int num_spaces = MAX - CountDigits(num);
  for (int i = 0; i < num_spaces; i++)
    spaces = spaces + " ";
  return spaces;
}

bool is_number(const string& s) {
  string::const_iterator it = s.begin();
  while (it != s.end() && isdigit(*it)) ++it;
  return !s.empty() && it == s.end();
}
}
