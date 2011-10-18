#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <map>
#include <ctime>
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
int alignment( const char* query, const char* base, int len_query, int len_base); 

// evaluation function for alignment
int score(char ch1, char ch2);

// the intergrated solution
void enquery(char* file_name, char* query);

// output the result alignment details
void alignment_print(const char* S, const char* T);

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
  start_time = time(NULL);
  enquery("NC_004347.faa", query);
  enquery("NC_004349.faa", query);
  end_time = time(NULL);
  
  printf("time elapsed for community prokaryote %f\n", difftime(end_time, start_time));
  printf("best score: %d\nbest result: %s\n\n", best_score, best_result.c_str());
  printf("result as follows: \n %s\n", best_meta.c_str());
  alignment_print(query, best_result.c_str());
  
  best_score = 0;
  best_result.clear();
  best_meta.clear();
  printf("\n***************************************************************************************\n");

  // my own prokaryote
  start_time = time(NULL);
  enquery("NC_010622.faa", query);
  enquery("NC_010623.faa", query);
  enquery("NC_010625.faa", query);
  enquery("NC_010627.faa", query);	
  end_time = time(NULL);

  printf("time elapsed for my own prokaryote %f\n", difftime(end_time, start_time));
  printf("best score: %d\nbest result: %s\n\n", best_score, best_result.c_str());
  printf("result as follows: \n %s\n", best_meta.c_str());
  alignment_print(query, best_result.c_str());


  //alignment_print("abcxdex", "xxxcde");

  return 0;
}



void enquery(char* file_name, char* query)
{
  ifstream fin;
  fin.open(file_name);

  string input;
  string sequence;

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
	      int score = alignment(query, sequence.c_str(), strlen(query), sequence.size());

	      if (score > best_score)
		{
		  best_score = score;
		  best_result = sequence;					
		  best_meta = input;
		}


	      //printf( "file: %s, protein #%d, current score: %d, best score: %d\n", file_name, ++num_protein, score, best_score );

	      sequence.clear();
	    }
	}

      else
	sequence += input;
    }
	
}



// actual local alignment algorithm
int alignment( const char* query, const char* base, int len_query, int len_base) 
{
  //int* ptr0 = (int*)malloc(sizeof(int)*(len_query+1));
  //int* ptr1 = (int*)malloc(sizeof(int)*(len_query+1));

  int *ptr = (int*)malloc(sizeof(int)*2*(len_query+1));
  int *ptr0 = ptr;
  int *ptr1 = ptr0 + len_query + 1;

  //memset(ptr0, 0, sizeof(int)*(len_query+1));
  //memset(ptr1, 0, sizeof(int)*(len_query+1));
  memset(ptr, 0, sizeof(int)*2*(len_query+1));

  for ( int j = 1; j<=len_base; ++j )
    {
      for ( int i=1; i<=len_query; ++i )
	{
	  // left, up and diagnol values
	  int vl = *( ptr1 + i - 1 ) + score( *( query + i - 1 ) , '-');
	  int vu = *( ptr0 + i ) + score( '-', *( base + j - 1 ) );
	  int vd = *( ptr0 + i - 1 ) + score( *( query + i - 1 ), *( base + j - 1 ) );

	  *(ptr1+i) = max( max(max(vl, vu), max(vu, vd)), max(max(vd, 0), max(0, vl)));		


#if pruning_level == 1
	  if ( *(ptr1+i) +11*max(len_query - i, len_base - j)  < best_score)
	    return -1;
#endif
	  
#if pruning_level == 2
	  if ( *(ptr1+i) +7*max(len_query - i, len_base - j)  < best_score)
	    return -1;
#endif
	}

      // after one cycle, switch ptr0 and ptr1
      int* tmp = ptr1;
      ptr1 = ptr0;
      ptr0 = tmp;
    }

  int ret = 0;
  /*if ( len_base % 2 )
    ret = *( ptr1 + len_query );
    else
    ret = *( ptr0 + len_query );*/
  ret = *( ptr + len_query + ( len_base%2 ) * ( len_query + 1 ) );

  free(ptr);
  //free(ptr1);
  return ret;
}


int score(char ch1, char ch2)
{
  if (ch1 == '-' || ch2 == '-')
    return -4;

  int x = map_index_blosum[ch1];
  int y = map_index_blosum[ch2];

  if ( map_index_blosum.count(ch1) == 0 )
    x = 22;
  if ( map_index_blosum.count(ch2) == 0 )
    y = 22;

  return BLOSUM62[ y*size_blosum + x];

  //if ( ch1 == ch2 )
  //	return 2;
  //else
  //	return -1;
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
      //*( reconstr + j * ( size_S + 1 ) ) = 'u';
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
  for ( int j =1; j <= size_T; ++j )
    {
      for ( int i = 1; i <= size_S; ++i )
	{
	  // left, up and diagnol values
	  int vl = *( ptr1 + i - 1 ) + score( *( S + i - 1 ) , '-');
	  int vu = *( ptr0 + i ) + score( '-', *( T + j - 1 ) );
	  int vd = *( ptr0 + i - 1 ) + score( *( S + i - 1 ), *( T + j - 1 ) );

	  *(ptr1+i) = max( max(max(vl, vu), max(vu, vd)), max(max(vd, 0), max(0, vl)));		
			
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

  int local_score = *( ptr + size_S + ( size_T%2 ) * ( size_S + 1 ) );

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

  printf("score:%d\n Q %d %s %d\n     %s\n S %d %s %d\n", local_score, idx_S + 2, al_S.c_str(), end_S+1, al_mid.c_str(), idx_T+2,al_T.c_str(), end_T+1);
}
