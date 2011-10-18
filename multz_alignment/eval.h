#ifndef _EVAL_H_
#define _EVAL_H_


#include <stdio.h>
#include <stdlib.h>
#include <string.h> // remember to add this!!
#include <time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <bitset>

using namespace std;

// evaluation function for an alignment
// return the max length of ultra conservative element
// return 0 if we cannot find any
// percent stores the percentage of conservative column
int eval(vector<string> &data_set, int &start, int &specie, int specie_thresh, int len_thresh);

// the linear algorithm for calculating the longest uc elem
int get_ucelem(int* arr, int size, int &start);




#endif