/*
 * Panagiotis Skrimponis
 */

#pragma once
#include <vector>
#include <string>
#include <memory>
#include <complex>
#include <fstream>
#include <sstream> 
#include <iterator>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <sys/time.h>
#include <math.h>

enum OPTIONS_ENUM {
  DC=1, AC, TRAN, OPTIONS, PLOT, PRINT, SPD, ITER, SPARSE, METHOD, ITOL
};

#define str_toupper(str) {for (auto &it : str) it = toupper(it);}
#define print_vector(v, str) {std::cout << "\t\t-- "<< str <<" --\n"; for(unsigned int i = 1; i < vectorSize; i++) std::cout << v[i] << std::endl;}
#define read_line() {\
						getline(iFile, token);\
						std::replace (token.begin(), token.end(), '(' , ' ');\
						std::replace (token.begin(), token.end(), ')' , ' ');\
						std::replace (token.begin(), token.end(), ',' , ' ');\
						std::replace (token.begin(), token.end(), '=' , ' ');\
						last = token.find_last_not_of(' ');\
						token = token.substr(0, last+1);\
						str_toupper(token)\
					}
#include "cs_sparse_matrix.h"
#include "cs_complex_sparse_matrix.h"
#include "cs_matrix.h"