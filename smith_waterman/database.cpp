#include "database.h"
#include <cmath>
#include <cstring>
#include <list>

// inline istream& operator>> (istream& in, Point& p)
// {
//   return in>>p.x>>p.y>>p.z>>p.pitch>>p.roll>>p.yaw;
// }


// inline ostream& operator<< (ostream& out, Point& p)
// {
//   return out<<p.x<<" "
// 	    <<p.y<<" "
// 	    <<p.z<<" "
// 	    <<p.pitch<<" "
// 	    <<p.roll<<" "
// 	    <<p.yaw;  
// }


// inline istream& operator>> (istream& in, Gesture& g)
// {
//   return in;
// }


// inline ostream& operator<< (ostream& out, Gesture& g)
// {
//   switch (g.type_gesture)
//     {
//     case PUSH:
//       out<<"@ PUSH"<<endl;
//       break;
//     case WAVE:
//       out<<"@ WAVE"<<endl;
//       break;
//     default:
//       out<<"@ UNKNOWN"<<endl;
//       break;
//     }
  
//   for (vector<Point*>::iterator iter = g.data.begin(); iter != g.data.end(); ++iter)
//     {
//       out<<**iter<<endl;
//     }
//   return out<<endl;
// }




inline ostream& operator<< (ostream& out, GestureDB& db)
{
  Gesture *ptr_g;
  int count = 0;
  
  // we also want to add meta data to the record
  while ( (ptr_g = db.findNext()) && (++count) )
    out<<"# gesture #"<<count<<endl<<*ptr_g<<endl;
  return out;
}


// // visualize the gesture
// void Gesture::draw(){}


// // append a new point at the end of the gesture sequence
// void Gesture::append(Point* point)
// {
//   data.push_back(point);
// }
  

/////////////////////////////////////////////////
/// Distance
/////////////////////////////////////////////////
// inline float dist_l1(Gesture* g1, Gesture* g2)
// {
//   float res = 0.0f;
//   for (int i=0; i<min(g1->numPoint(), g2->numPoint()); ++i)
//     {
//       res += abs((*g1)[i] - (*g2)[i]);
//     }

//   return res;  
// }


// inline float dist_l2(Gesture* g1, Gesture* g2)
// {
//   float res = 0.0f;
//   for (int i=0; i<min(g1->numPoint(), g2->numPoint()); ++i)
//     {
//       res += pow((*g1)[i] - (*g2)[i], 2);
//     }
  
//   return sqrt(res);  
// }
  

// inline float score_l1(Point* p1, Point* p2)
// {
//   return ( abs(p1->x - p2->x) + 
// 	   abs(p1->y - p2->y) +
// 	   abs(p1->z - p2->z) +
// 	   abs(p1->pitch - p2->pitch) + 
// 	   abs(p1->roll - p2->roll) + 
// 	   abs(p1->yaw - p2->yaw) );
// }


// inline float score_l2(Point* p1, Point* p2)
// {
//   return 1/(1 + exp(sqrt( pow(p1->x - p2->x, 2) + 
// 			  pow(p1->y - p2->y, 2) +
// 			  pow(p1->z - p2->z, 2) +
// 			  pow(p1->pitch - p2->pitch, 2) + 
// 			  pow(p1->roll - p2->roll, 2) + 
// 			  pow(p1->yaw - p2->yaw, 2) )));
//}

// inline float score(Point* p1, Point* p2)
// {
//   return (p1 == NULL || p2 == NULL)? -1 : score_l2(p1, p2);
// }

float smith_waterman(Gesture* g1, Gesture* g2)
{
  int len1 = g1->numPoint();
  int len2 = g2->numPoint();
  float *score_buf = new float[len1];
  memset((void*)score_buf, 0, len1*sizeof(float)); 
  //for (int i=0; i<len1; score_buf[i]=4e7, ++i);

  float prev_l = 0;
  float best = 0;
  for ( int j=1; j<len2; ++j )
    for ( int i=1; i<len1; ++i )
      {	
	float tmp = max( max( score( (*g1)[i], (*g2)[j-1] ) + score_buf[i],
			      score( (*g1)[i-1], (*g2)[j] ) + prev_l ),
			 score( (*g1)[i-1], (*g1)[j-1] ) + score_buf[i-1] );
	score_buf[i-1] = prev_l;
	prev_l = tmp;
	best = max(best, tmp);
      }
  return best;
}



////////////////////////////////////////////////////////////////////////
// laod and store the database on disk
////////////////////////////////////////////////////////////////////////

inline string chomp(string& str)
{                                                                                
  int begin_idx = str.find_first_not_of(" \t");                                 
  int length = str.find_last_not_of(" \t") - begin_idx + 1;                     
  str = str.substr(begin_idx, length);                                          
  return str;
}

void SeqDB::load()
{
  ifstream fin;
  fin.open(file_name);
  
  string line;
  bool is_reading = 0;
  //stringstream ss ( stringstream::in | stringstream::out );
  // Gesture* ptr_g = new Gesture();
  // Point* ptr_p = new Point();	 

  string sequence;
  while (getline(fin, line))
    {
      // an empty line marks the end of a sequnce record
      if ( line.size() == 0 || chomp(line).size() == 0 )
	{
	  // if there are more empty lines in the beginning
	  // we just don't care	  
	  if ( is_reading )
	    {
	      // add the data into our database 
	      //data.push_back(ptr_g);
	      _seq.push_back(sequence);
	      sequence.clear();
	      is_reading = false;
	    }
	 	  
	  continue;
	}

      // # for comment and meta data for human reader
      if ( line[0] == '#' )
	continue; 

      // this marks the entrance of a new input
      if ( line[0] == '>' )
	{
	  // start a new sequence
	  is_reading = true;
	  _meta.push_back(line);
	  sequence.clear();
	}
      else if ( is_reading )
	sequence += input;

      // if we are in the middle of a record
      //ss << line;
      //ss >> *ptr_p;
      //ss.clear(); // this is vital...!!!
      //cout<<*ptr_p<<endl;
     
      // // append
      // ptr_g->append(ptr_p);	  
      // ptr_p = new Point();	  
    }  

  // there might not be an extra line in the end
  if ( is_reading )
      _seq.push_back(sequence);      

  //cout<<*this<<endl;

  fin.close();
}
  

void SeqDB::save()
{
  ofstream fout;
  // all the previously loaded examples will also be written back
  fout.open(file_name);
  fout<<*this<<endl;
  fout.close();
}


// void SeqDB::append(Gesture* g)
// {
//   data.push_back(g);
// }
 

// void SeqDB::remove()
// {
  
// }


// // find the next element in the sequence
// Gesture* SeqDB::findNext()
// {
//   static size_t curr_idx = 0;
//   if ( curr_idx >= data.size() )
//     {
//       curr_idx = 0;
//       return NULL;
//     }
//   return data[curr_idx++];
// }


inline float distance(Gesture* g1, Gesture* g2)
{
  //return dist_l2(g1, g2);
  return smith_waterman(g1,g2);
}


// find the k-nearest neighbors
GestureType SeqDB::kNN(Gesture* query, int k)
{
  //map<float, Gesture*> knn_map;
  list<Gesture*> knn;
  list<float> weight;
  
  printf("\t\t database size: %d\n", data.size());
  Gesture* curr_ges = data[0]; //findNext();
  knn.push_front(curr_ges);
  weight.push_front( distance(query, curr_ges) );     

  //while ( curr_ges = findNext() )
  for ( int idx=1; idx<data.size(); ++idx )
    {
      curr_ges = data[idx];

      // insert the gesture to the appropriate position
      float curr_dist = distance(query, curr_ges);
      list<float>::iterator weight_iter = weight.begin();
      list<Gesture*>::iterator knn_iter = knn.begin();	  
      for (; knn_iter != knn.end(); ++knn_iter, ++weight_iter )	  
	if ( *weight_iter < curr_dist )
	  break;	  

      if ( weight_iter == weight.begin() )
	{
	  weight.push_front(curr_dist);
	  knn.push_front(curr_ges);
	}
      else
	{
	  weight.insert(--weight_iter, curr_dist);
	  knn.insert(--knn_iter, curr_ges);
	}

      if ( knn.size() > k )
	{
	  weight.pop_back();
	  knn.pop_back();
	}
    }
  
  int count[NUM_GESTURETYPE] = {0};
  list<float>::iterator weight_iter = weight.begin(); 
  list<Gesture*>::iterator knn_iter = knn.begin();
  for (; knn_iter != knn.end(); ++knn_iter, ++weight_iter )	  
    {
      printf("\t\t kNN: gesture type %d, smith_waterman %f\n", 
	     (int)(*knn_iter)->getType(), *weight_iter);
      
      ++count[(int)(*knn_iter)->getType()];
    }

  GestureType majority = (GestureType) 0;
  int max_count = 0;
  for ( int i=0; i<NUM_GESTURETYPE; ++i )
    {
      if ( max_count < count[i] )
	{
	  max_count = count[i];
	  majority = (GestureType) i;
	}
    }
  return majority;
}


// // // test the result
// int main()
// {
//   Point* ptr_p = new Point();//(1,2,3,4,54,6);
//   Gesture* ptr_g = new Gesture();
//   SeqDB* ptr_db = new SeqDB("tmp_gesture.db");
//   ptr_db->load();  

//   // cin>>*ptr_p;
//   // cout<<*ptr_p<<endl;
//   // delete ptr_p;
  
//   for ( int k=0; k<2; ++k )
//     {
//       for ( int i=0; i<3; ++i)
// 	{
// 	  ptr_p = new Point(1,2,3,4,54,6);
// 	  ptr_g->append(ptr_p);
// 	}
//       ptr_db->append(ptr_g);
//       ptr_g = new Gesture();
//     }  

//   //cout<<*ptr_db;
//   ptr_db->save();
//   //cout<<*ptr_g;
// }


