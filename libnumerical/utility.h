//
// Created by Todd Dupont
//

#ifndef CS28515PROJ1_UTILITY_H
#define CS28515PROJ1_UTILITY_H

#include <cstdlib>
#include <iostream>
#include <string>

void error(std::string s);
void warning(std::string s, std::ostream &out = std::cerr);

char *getCmdOption(char **begin, char **end, const std::string &option);
bool cmdOptionExists(char **begin, char **end, const std::string &option);

#endif //CS28515PROJ1_UTILITY_H
