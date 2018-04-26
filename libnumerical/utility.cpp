//
// Created by Todd Dupont
//

#include "utility.h"

void error(std::string s) {// write message and die
  std::cerr << "*** Error ***\n";
  std::cerr << s << std::endl;
  exit(1);
}

void warning(std::string s, std::ostream &out) {// write message but keep cranking
  out << "**warning: " << s << std::endl;
}

char *getCmdOption(char **begin, char **end, const std::string &option) {
  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) {
    return *itr;
  }
  return 0;
}

bool cmdOptionExists(char **begin, char **end, const std::string &option) {
  return std::find(begin, end, option) != end;
}