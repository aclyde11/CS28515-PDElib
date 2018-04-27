//
// Created by Todd Dupont
//

#ifndef CS28515PROJ1_UTILITY_H
#define CS28515PROJ1_UTILITY_H

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

void error(std::string s);
void warning(std::string s, std::ostream &out = std::cerr);

char *getCmdOption(char **begin, char **end, const std::string &option);
bool cmdOptionExists(char **begin, char **end, const std::string &option);

/*
 * Writes PDE parameters out to file for python plotting
 */
void writeParams(std::string name, std::vector<std::string> params);

/*
 * Writes vector U at t to data file
 */
void writeUpdateStep(std::string name, std::vector<double> a, double t);

#endif //CS28515PROJ1_UTILITY_H
