//
// Created by Todd Dupont
//

#include "utility.h"
#include "numvec.h"

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

void writeParams(std::string name, std::vector<std::string> params) {
  std::ofstream ofs(name);

  if (!ofs) {
    std::cout << "Error opening file for output" << std::endl;
    return;
  }
  for (int i = 0; i < params.size(); i++) {
    ofs << params[i];
    if (i != params.size() - 1)
      ofs << "\t";
  }
  ofs << std::endl;
  ofs.close();
}

void writeUpdateStep(std::string name, std::vector<double> a, double t) {
  std::ofstream ofs(name, std::fstream::app);

  if (!ofs) {
    std::cout << "Error opening file for output" << std::endl;
    return;
  }
  ofs << t << '\t';
  ofs << a;
  ofs.close();
}