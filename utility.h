//
// Created by Austin Clyde on 4/4/18.
//

#ifndef CS28515PROJ1_UTILITY_H
#define CS28515PROJ1_UTILITY_H

#include <cstdlib>
#include <iostream>
#include <string>

void error( std::string s );
void warning( std::string s, std::ostream& out=std::cerr );
inline double sq( double t ) {return t*t;}

#endif //CS28515PROJ1_UTILITY_H
