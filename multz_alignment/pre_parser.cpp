// preparser for multz alignment

#include <stdio.h>
#include <stdlib.h>
#include <string.h> // remember to add this!!
#include <time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>

using namespace std;

//************************global variables***************************
const int specie_thresh = 40;
const int len_thresh = 100;
string data_path = "/projects/instr/09wi/nobackup/cse427/multiz44way/";
//string data_path = "";

//*********************function declarations*************************
// the intergrated solution for pre-parsing
// always call this function
int pre_parser();

// pre-processing of genome alignment files
int pre_parsing(const char* file_name, const char* out_file );


// change this function for settings
int pre_parser(int range_start, int range_end)
{
  // valid block counter
  int num_block = 0;

  // for the 24 human genomes
  for (int i = range_start; i <= range_end; ++i)
    {
      time_t start=0, end=0;
      string input_name, output_name;
      string name_num;
      if ( i < 10 )
        {
          name_num = char(i+int('0'));
        }
      else if ( i < 23 )
        {
          name_num.push_back(char(i/10 + int('0')));
          name_num.push_back(char(i%10 + int('0')));
        }
      else if ( i == 23 )
        {
          name_num = "X";
        }
      else
        {
          name_num = "Y";
        }

      //*****************input/output names*******************
      input_name = data_path + "chr" + name_num + ".maf";
      /*output_name = "chr" + name_num + "_preproc.txt";*/
      output_name = "chr" + name_num + "_preproc.txt";

      start = time(NULL);
      num_block += pre_parsing(input_name.c_str(), output_name.c_str());
      end = time(NULL);

      printf("time spent for %s is %f sec\n",
             input_name.c_str(),difftime(end, start));
    }
  return num_block;
}

// pre-processing of genome alignment files
int pre_parsing(const char* file_name, const char* out_file)
{

  ifstream fin;
  fin.open(file_name);
  ofstream fout;
  fout.open(out_file);

  int num_specie = 0;
  int align_len = 0;
  int align_block_num = 0;

  bool is_human = true;
  bool valid_read = true;

  vector<string> block_buffer;
  string input;

  while ( getline(fin, input) )
    {

      // if we met an empty line, do the data filtering
      if ( !input.size() )
        {
          // reject the cases that are not valid
          if ( num_specie < specie_thresh || !valid_read )
            {
              // do nothing
            }

          // valid cases
          else
            {
              ++align_block_num;
	      
	      cout<<"#block "<<align_block_num
		  <<", length "<<align_len
		  <<", #species "<<num_specie<<endl;
	      
	      // print out the original alignment block
	      for (int i=0; i<block_buffer.size() || (fout<<endl && 0); ++i)
		{
		  fout<<block_buffer[i]<<endl;
		}
		
	    }
	  
          // reset the counters and flags
          num_specie = 0;
          align_len = 0;
          if ( block_buffer.size() )
	    block_buffer.clear();
          valid_read = true;
          is_human = true;
        }

      // if the line starts with an s, read in the data
      else if ( valid_read && input[0] == 's' )
        {
          if ( is_human )
            {
	      is_human = false;

              char is_s;
              int plus_start;
              int length = 0;
              char sign;
              int minus_start;
              char* align_specie = (char*)malloc(32*sizeof(char));
	      char* line_buffer = (char*)malloc(input.size()*sizeof(char));

              sscanf(input.c_str(), "%c %s %d %d %c %d %s",
                     &is_s,
                     align_specie,
                     &plus_start,
                     &length,
                     &sign,
                     &minus_start,
                     line_buffer) ;
              if ( (align_len = strlen(line_buffer)) < len_thresh )
                valid_read = false;

	      free(align_specie);
	      free(line_buffer);
            }
	  
	  if ( valid_read )
	    {
	      block_buffer.push_back(input);
	      ++num_specie;
	    }
	  
        }
    }

  
  if (block_buffer.size()) // clear up the rest
    {
      if ( num_specie < specie_thresh || !valid_read )
	{
	  // do nothing
	}
      
      // valid cases
      else
	{
	  ++align_block_num;
	  
	  cout<<"#block "<<align_block_num
	      <<", length "<<align_len
	      <<", #species "<<num_specie<<endl;
	  
	  // print out the original alignment block
	  for (int i=0; i<block_buffer.size() || (fout<<endl && 0); ++i)
	    {
	      fout<<block_buffer[i]<<endl;
	    }
	}
    }
  
  return align_block_num;
  
}


int to_int(char* num)
{
  int size = strlen(num);
  int ret = 0;
  for (int i=0; i<size; ++i)
    {
      ret = ret*10 + (int)(num[i]-'0');
    }
  return ret;
}


int main(int argc, char** argv)
{
  if ( argc!=3 )
    printf("usage: pre_parser <start chromsome> <end chromsome>\n(from 1 to 24)\n");
  else
    {
      int a = to_int(argv[1]);
      int b = to_int(argv[2]);
      printf("# valid blocks: %d\n", pre_parser(a, b));
    }
  return 0;
}
