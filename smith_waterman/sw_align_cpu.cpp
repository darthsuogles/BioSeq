#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <map>
#include <ctime>
#include <sys/time.h>
using namespace std;

// debug mode: set it true if you want to see the calculation steps
#define debug_mode false   

//************************ extra credit modes *********************
// pruning level: 
// stopping condition for each pairwise local alignment
//			0 for no pruning
//          1 for safe pruning
//          2 for pruning using heuristics
#define pruning_level 0

//*******************global variable declarations**************************

int g_mismatch_score = -4;

// BLOSUM62 is used for evaluation
const int size_blosum  = 24;
int BLOSUM62[size_blosum*size_blosum] = {
  4,  -1,  -2,  -2,   0,  -1,  -1,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   0,  -3,  -2,   0,  -2,  -1,   0,  -4,   
  -1,   5,   0,  -2,  -3,   1,   0,  -2,   0,  -3,  -2,   2,  -1,  -3,  -2,  -1,  -1,  -3,  -2,  -3,  -1,   0,  -1,  -4,   
  -2,   0,   6,   1,  -3,   0,   0,   0,   1,  -3,  -3,   0,  -2,  -3,  -2,   1,   0,  -4,  -2,  -3,   3,   0,  -1,  -4,   
  -2,  -2,   1,   6,  -3,   0,   2,  -1,  -1,  -3,  -4,  -1,  -3,  -3,  -1,   0,  -1,  -4,  -3,  -3,   4,   1,  -1,  -4,   
  0,  -3,  -3,  -3,   9,  -3,  -4,  -3,  -3,  -1,  -1,  -3,  -1,  -2,  -3,  -1,  -1,  -2,  -2,  -1,  -3,  -3,  -2,  -4,   
  -1,   1,   0,   0,  -3,   5,   2,  -2,   0,  -3,  -2,   1,   0,  -3,  -1,   0,  -1,  -2,  -1,  -2,   0,   3,  -1,  -4,   
  -1,   0,   0,   2,  -4,   2,   5,  -2,   0,  -3,  -3,   1,  -2,  -3,  -1,   0,  -1,  -3,  -2,  -2,   1,   4,  -1,  -4,   
  0,  -2,   0,  -1,  -3,  -2,  -2,   6,  -2,  -4,  -4,  -2,  -3,  -3,  -2,   0,  -2,  -2,  -3,  -3,  -1,  -2,  -1,  -4,   
  -2,   0,   1,  -1,  -3,   0,   0,  -2,   8,  -3,  -3,  -1,  -2,  -1,  -2,  -1,  -2,  -2,   2,  -3,   0,   0,  -1,  -4,   
  -1,  -3,  -3,  -3,  -1,  -3,  -3,  -4,  -3,   4,   2,  -3,   1,   0,  -3,  -2,  -1,  -3,  -1,   3,  -3,  -3,  -1,  -4,   
  -1,  -2,  -3,  -4,  -1,  -2,  -3,  -4,  -3,   2,   4,  -2,   2,   0,  -3,  -2,  -1,  -2,  -1,   1,  -4,  -3,  -1,  -4,   
  -1,   2,   0,  -1,  -3,   1,   1,  -2,  -1,  -3,  -2,   5,  -1,  -3,  -1,   0,  -1,  -3,  -2,  -2,   0,   1,  -1,  -4,   
  -1,  -1,  -2,  -3,  -1,   0,  -2,  -3,  -2,   1,   2,  -1,   5,   0,  -2,  -1,  -1,  -1,  -1,   1,  -3,  -1,  -1,  -4,   
  -2,  -3,  -3,  -3,  -2,  -3,  -3,  -3,  -1,   0,   0,  -3,   0,   6,  -4,  -2,  -2,   1,   3,  -1,  -3,  -3,  -1,  -4,   
  -1,  -2,  -2,  -1,  -3,  -1,  -1,  -2,  -2,  -3,  -3,  -1,  -2,  -4,   7,  -1,  -1,  -4,  -3,  -2,  -2,  -1,  -2,  -4,   
  1,  -1,   1,   0,  -1,   0,   0,   0,  -1,  -2,  -2,   0,  -1,  -2,  -1,   4,   1,  -3,  -2,  -2,   0,   0,   0,  -4,   
  0,  -1,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   5,  -2,  -2,   0,  -1,  -1,   0,  -4,   
  -3,  -3,  -4,  -4,  -2,  -2,  -3,  -2,  -2,  -3,  -2,  -3,  -1,   1,  -4,  -3,  -2,  11,   2,  -3,  -4,  -3,  -2,  -4,   
  -2,  -2,  -2,  -3,  -2,  -1,  -2,  -3,   2,  -1,  -1,  -2,  -1,   3,  -3,  -2,  -2,   2,   7,  -1,  -3,  -2,  -1,  -4,   
  0,  -3,  -3,  -3,  -1,  -2,  -2,  -3,  -3,   3,   1,  -2,   1,  -1,  -2,  -2,   0,  -3,  -1,   4,  -3,  -2,  -1,  -4,   
  -2,  -1,   3,   4,  -3,   0,   1,  -1,   0,  -3,  -4,   0,  -3,  -3,  -2,   0,  -1,  -4,  -3,  -3,   4,   1,  -1,  -4,   
  -1,   0,   0,   1,  -3,   3,   4,  -2,   0,  -3,  -3,   1,  -1,  -3,  -1,   0,  -1,  -3,  -2,  -2,   1,   4,  -1,  -4,   
  0,  -1,  -1,  -1,  -2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -2,   0,   0,  -2,  -1,  -1,  -1,  -1,  -1,  -4,   
  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,   1
};


// character to index mapping
map<char, int> map_index_blosum;

// best numerical score
int best_score = 0;

// best protein sequence for matching
string best_result;

// best protein's metadata
string best_meta;

//************************function declarations****************************

// local alignment function for string
//int alignment( const char* query, const char* base, int len_query, int len_base); 
int sw_align_kernel(string S, string T);

// evaluation function for alignment
int score(char ch1, char ch2);

// the intergrated solution
void enquery(char* file_name, char* query);

// output the result alignment details
void alignment_print(const char* S, const char* T);

double getElapsedTime(timeval& t1, timeval& t2)
{
  double elapsedTime;
  elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
  return elapsedTime;
}

int main()
{
  time_t start_time ,end_time;

  // load the index map into main memory
  map_index_blosum['A']=0;
  map_index_blosum['R']=1;
  map_index_blosum['N']=2;
  map_index_blosum['D']=3;
  map_index_blosum['C']=4;
  map_index_blosum['Q']=5;
  map_index_blosum['E']=6;
  map_index_blosum['G']=7;
  map_index_blosum['H']=8;
  map_index_blosum['I']=9;
  map_index_blosum['L']=10;
  map_index_blosum['K']=11;
  map_index_blosum['M']=12;
  map_index_blosum['F']=13;
  map_index_blosum['P']=14;
  map_index_blosum['S']=15;
  map_index_blosum['T']=16;
  map_index_blosum['W']=17;
  map_index_blosum['Y']=18;
  map_index_blosum['V']=19;
  map_index_blosum['B']=20;
  map_index_blosum['Z']=21;
  map_index_blosum['X']=22;
  map_index_blosum['*']=23;
    
  char* query = "MTHQSHAYHMVKPSPWPLTGALSALLMTSGLAMWFHFHSMTLLMLGLLTNTLTMYQWWRDVTRESTYQGHHTPPVQKGLRYGMILFITSEVFFFAGFFWAFYHSSLAPTPQLGGHWPPTGITPLNPLEVPLLNTSVLLASGVSITWAHHSLMENNRNQMIQALLITILLGLYFTLLQASEYFESPFTISDGIYGSTFFVATGFHGLHVIIGSTFLTICFIRQLMFHFTSKHHFGFEAAAWYWHFVDVVWLFLYVSIYWWG";

  // community prokaryote
  timeval t1, t2;

  //start_time = time(NULL);
  gettimeofday(&t1, NULL);
  enquery("NC_004347.faa", query);
  enquery("NC_004349.faa", query);
  //end_time = time(NULL);
  gettimeofday(&t2, NULL);
  
  printf("time elapsed for community prokaryote %.4lf s\n", getElapsedTime(t1, t2)/1000.0);
  printf("best score: %d\nbest result: %s\n\n", best_score, best_result.c_str());
  printf("result as follows: \n %s\n", best_meta.c_str());
  alignment_print(query, best_result.c_str());
  
  best_score = 0;
  best_result.clear();
  best_meta.clear();
  printf("\n***************************************************************************************\n");

  // my own prokaryote
  gettimeofday(&t1, NULL);
  enquery("NC_010622.faa", query);
  enquery("NC_010623.faa", query);
  enquery("NC_010625.faa", query);
  enquery("NC_010627.faa", query);	
  gettimeofday(&t2, NULL);
  
  printf("time elapsed for my own prokaryote %.4lf s\n", getElapsedTime(t1, t2)/1000.0);
  printf("best score: %d\nbest result: %s\n\n", best_score, best_result.c_str());
  printf("result as follows: \n %s\n", best_meta.c_str());
  alignment_print(query, best_result.c_str());

  return 0;
}


void enquery(char* file_name, char* query)
{
  ifstream fin;
  fin.open(file_name);

  string input;
  string sequence;
  string seq_meta;

  int num_protein = 0;
  while ( getline(fin, input) )
    {
      if (input.size() == 0)
	break;

      // the starting line of a protein coding
      if ( input[0] == '>' )
	{
	  // it is possible that an protein coding is just finished
	  if ( sequence.size() > 0 )
	    {	
	      int score = sw_align_kernel(query, sequence);

	      if (score > best_score)
		{
		  best_score = score;
		  best_result = sequence;					
		  best_meta = seq_meta;
		}
	      
	      //printf( "file: %s, protein #%d, current score: %d, best score: %d\n", file_name, ++num_protein, score, best_score );
	      sequence.clear();
	    }

	  seq_meta = input; /// Store the meta data for next sequence
	}

      else
	sequence += input;
    }	
  
  /// Compare the last sequence and update the data
  int score = sw_align_kernel(query, sequence);  
  if (score > best_score)
    {
      best_score = score;
      best_result = sequence;					
      best_meta = seq_meta;
    }
}

// actual local alignment algorithm
//int alignment(const char* S, const char* T, int len_query, int len_base) 
int sw_align_kernel(string S, string T)
{
  int m = S.size();
  int n = T.size();

  int mismatch = g_mismatch_score;
  int prev = 0; /// We need one more element to store the previous data
  int* score_table = new int[n+1];
  memset(score_table, 0, sizeof(int)*(n+1));  
  
  int max_val = -1e9;
  char* sp = &S[0];
  for ( int i=0; i<m; ++i )
    {
      prev = 0; 
      char* tp = &T[0];
      for ( int j=1; j<=n; ++j )
	{
	  int vl = prev + mismatch;
	  int vd = score_table[j-1] + score(*sp, *tp++);
	  int vu = score_table[j] + mismatch;
	  
	  score_table[j-1] = prev;
	  prev = max(max(vl, 0), max(vd, vu));	
	  max_val = max(max_val, prev);
	}
      ++sp;
    }
  
  delete [] score_table;
  return max_val;
}


int score(char ch1, char ch2)
{
  if (ch1 == '-' || ch2 == '-')
    return -4;

  int x, y;

  if ( map_index_blosum.count(ch1) == 0 )
    x = 22;
  else
    x = map_index_blosum[ch1];
  if ( map_index_blosum.count(ch2) == 0 )
    y = 22;
  else
    y = map_index_blosum[ch2];

  return BLOSUM62[ y*size_blosum + x];
}


void alignment_print(const char* S, const char* T)
{
  int size_S = strlen(S);
  int size_T = strlen(T);

  int* ptr = (int*)malloc( sizeof(int) * 2 * (size_S+1) );
  int* ptr0 = ptr;
  int* ptr1 = ptr + size_S + 1;

  memset( ptr, 0, sizeof(int) * 2 * (size_S+1) );

  // reconstruction table
  char* reconstr = (char*)malloc( sizeof(char) * (size_S + 1) * (size_T + 1) );

  *reconstr = 'e';
  for ( int j = 0; j <= size_T; ++j )
    {
      for ( int i = 0; i <= size_S; ++i )
	{
	  if ( i == 0 && j != 0 )
	    *( reconstr + j * ( size_S + 1 ) ) = 'u';
	  else if (i!=0 && j == 0 )
	    *( reconstr +i ) = 'l';
	  else
	    continue;
	}
    }

  // do the dp and update the reconstruction table
  int max_score = -1e9;
  for ( int j =1; j <= size_T; ++j )
    {
      for ( int i = 1; i <= size_S; ++i )
	{
	  // left, up and diagnol values
	  int vl = *( ptr1 + i - 1 ) + g_mismatch_score;
	  int vu = *( ptr0 + i ) + g_mismatch_score;
	  int vd = *( ptr0 + i - 1 ) + score( *( S + i - 1 ), *( T + j - 1 ) );

	  /// Keep a record of the maximum score seen so far
	  max_score = max( max_score, *(ptr1+i) = max( max(vd, 0), max(vl, vu) ) );
	  
	  if ( vl ==  *(ptr1+i) )
	    *( reconstr + j * (size_S + 1) + i ) = 'l' ;
	  if ( vu ==  *(ptr1+i) )
	    {
	      *( reconstr + j * (size_S + 1) + i ) = 'u';
	    }
	  if ( vd ==  *(ptr1+i) )
	    {
	      *( reconstr + j * (size_S + 1) + i )  = 'd';
	    }
	  if ( 0 == *(ptr1 + i) )
	    {
	      *( reconstr + j * (size_S + 1) + i ) = 'e';
	    }
	}
      int *tmp = ptr0;
      ptr0 = ptr1;
      ptr1 = tmp;
    }
  
  // reconstruction
  char* iter = reconstr + (size_T + 1)*(size_S + 1) - 1;
  int idx_T = size_T - 1;
  int idx_S = size_S - 1;
  int end_T = 0, end_S = 0;
  string al_S, al_T, al_mid;
  bool last_match = false;
	
  // the alignment should start from the last match
  while ( *iter != 'e' )
    {
      if ( *iter == 'd' )
	{
	  last_match = true;
	  al_S = al_S.insert(0, 1, S[idx_S]);
	  al_T = al_T.insert(0, 1, T[idx_T]);

	  if ( S[idx_S] == T[idx_T] ) 
	    al_mid = al_mid.insert(0, 1, T[idx_T]);
	  else if ( score(S[idx_S], T[idx_T]) > 0 )
	    al_mid = al_mid.insert(0, 1,  '+');
	  else
	    al_mid = al_mid.insert(0, 1,  ' ');		

	  if (end_T == 0)
	    end_T = idx_T ;
	  if (end_S == 0)
	    end_S = idx_S ;

	  --idx_T;
	  --idx_S;
	  iter -= size_S + 2 ;
	}

      else if ( *iter == 'l' )
	{
	  if ( last_match )
	    {
	      al_S = al_S.insert(0, 1, S[idx_S]);
	      al_T = al_T.insert(0, 1, '-');
	      al_mid = al_mid.insert(0, 1, ' ');
	    }
	  --idx_S;
	  iter -= 1;
	}

      else if ( *iter == 'u' )	
	{
	  if ( last_match )
	    {
	      al_T = al_T.insert(0, 1,T[idx_T]);
	      al_S = al_S.insert(0, 1,'-');
	      al_mid = al_mid.insert(0, 1, ' ');
	    }
	  --idx_T;
	  iter -= size_S + 1 ;
	}
      
      else
	printf("Error: the reconstruction table contains unidentified character -- %c\n", *iter);
    }

  string prefix_mid = "    ";
  int num_space = max(idx_S + 1, idx_T+1);
  if (num_space != 0)
    {
      while ( num_space != 0 )
	{
	  num_space /= 10;
	  prefix_mid += " ";
	}
    }
  else
    prefix_mid += " ";

  al_mid += prefix_mid;

  printf("score:%d\n Q %d %s %d\n     %s\n S %d %s %d\n", max_score, idx_S + 2, al_S.c_str(), end_S+1, al_mid.c_str(), idx_T+2,al_T.c_str(), end_T+1);
}
