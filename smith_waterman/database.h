#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>

using namespace std;

// define a database holding gestures
class SeqDB 
{
public:
  SeqDB(){}
  SeqDB(const char* file_name):file_name(file_name){}
  void load();
  void save();
  void append(string seq, string meta);
  void remove();  
  //  Gesture* findNext();  

  friend ostream& operator<< (ostream& out, SeqDB& db);

  // find the k-nearest neighbors
  //void findNearest(int k);  
  GestureType kNN(Gesture* query, int k);
  GestureType predict(Gesture *g);

private:
  vector<string> _seq;
  vector<string> _meta;
};


