/** 
*	given the pre-processed data, output the needed data structure for evaluation
*/
#ifndef _PARSER_H_
#define _PARSER_H_


#include <stdio.h>
#include <stdlib.h>
#include <string.h> // remember to add this!!
#include <time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// the integrated solution for parsing
// always call this function
int parser( const int specie_thresh, const int len_thresh, const double percent_thresh );

// parseing the preparsed files
// read the files into main memory as bitmap
int parsing( const char* file_name, const int specie_thresh, const int len_thresh, const double percent_thresh);



#endif