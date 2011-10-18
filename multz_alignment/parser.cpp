#include "parser.h"
#include "eval.h"

using namespace std;
enum {A, T, C, G};

ofstream deout("debug.txt");
ofstream fout("alignment.txt");


int validate( vector<string> data_set, int eval_start, int eval_len, int eval_num_specie, string human_coordinate )
{
	int eval_end = eval_start + eval_len - 1;
	// used as the ultra-conserve column marker
	char mark_char[1000] = {'\0'};
	int pos_count = 0;
	int remove_count = 0;

	for ( int j=eval_start; j<=eval_end; ++j )
	{
		bool is_space = true;
		int nuclitide_count[4] = {0};
		for ( int i=0; i<eval_num_specie; ++i )
		{
			char ch = data_set[i][j];
			switch(ch)
			{
			case 'A':
				++nuclitide_count[A];
				is_space = false;
				break;
			case 'T':
				++nuclitide_count[T];
				is_space = false;
				break;
			case 'C':
				++nuclitide_count[C];
				is_space = false;
				break;
			case 'G':
				++nuclitide_count[G];
				is_space = false;
				break;
			case '-':
				is_space = true;
			default:
				break;
			}
		}

		int max_count = max(	 max(	max(nuclitide_count[A], nuclitide_count[T]), 
			max(nuclitide_count[T], nuclitide_count[C]) ), 
			nuclitide_count[G] );

		if ( max_count == eval_num_specie )
		{
			mark_char[j] = '*';
			++pos_count;
		}
		else
			mark_char[j] = ' ';
	}

	double eval_percent = pos_count / (double) eval_len;

	fout<<human_coordinate<<endl
		<<eval_len<<" columns"<<endl
		<<eval_num_specie<<" species"<<endl
		<<eval_percent*100<<"% columns identical"<<endl;

	// print the alignment body
	for ( int i=0; i<data_set.size(); ++i )
	{
		for (int j=eval_start; j<=eval_end; ++j)
		{	
			if ( mark_char[j] != 'n' )
				fout<<data_set[i][j];
		}
		fout<<endl;
	}

	// print the star markers
	for (int k =eval_start; k<= eval_end; ++k)
	{
		if ( mark_char[k] != 'n' )
			fout<<mark_char[k];
	}
	fout<<endl<<endl;

	return 1;
}



// get human coordinate
string get_coord(string onum, int length, char sign)
{
	int size = onum.size();

	int num = 0;
	for (int i=size-5; i<size; ++i)
	{
		num = num*10 + (int)(onum[i]-'0');
	}

	if (sign == '+')
		num += length-1;
	else
		num -= length+1;

	string tmp1 = onum;
	for (int i=size-1; i>=size-5; --i)
	{
		tmp1[i] = (char)(num%10+'0');
		num /= 10;
	}

	return tmp1;
}

int parser( const int specie_thresh, const int len_thresh, const double percent_thresh )
{
	string data_path = "data/";

	int total_block_chosen = 0;

	for (int i=1; i<=22; ++i)
	{
		string num;
		if ( i < 10 )
			num.push_back( char(int('0') + i) );
		else
		{
			num.push_back( char(i/10 + int('0')) );
			num.push_back( char(i%10 + int('0')) );
		}
		total_block_chosen += parsing((data_path+"chr"+num+"_preproc.txt").c_str(), specie_thresh, len_thresh, percent_thresh);
		//total_block_chosen += parsing((data_path+"chr"+num+"_preproc44.txt").c_str(), specie_thresh, len_thresh, percent_thresh);
	}
	total_block_chosen += parsing((data_path+"chrX_preproc.txt").c_str(), specie_thresh, len_thresh, percent_thresh);
	//parsing((data_path+"chrX_preproc44.txt").c_str(), specie_thresh, len_thresh, percent_thresh);
	total_block_chosen += parsing((data_path+"chrY_preproc.txt").c_str(), specie_thresh, len_thresh, percent_thresh);
	//parsing((data_path+"chrY_preproc44.txt").c_str(), specie_thresh, len_thresh, percent_thresh);

	// debuggin
	//parsing((data_path+"debug_preproc.txt").c_str(), specie_thresh, len_thresh, percent_thresh);

	return total_block_chosen;
}


int parsing( const char* file_name, const int specie_thresh, const int len_thresh, const double percent_thresh)
{
	ifstream fin;
	fin.open(file_name);

	// figures to be collected
	string human_coordinate;

	// for the overall stat
	int total_num_block = 0;
	int chosen_num_block = 0;

	// validations
	bool read_on = true;
	bool is_human = true;
	vector<int> valid_count;

	// the stored data_set
	vector<string> data_set;

	// I/O stuff
	string input;

	// for record holder
	int best_align_len = 0;

	while ( getline( fin, input) )
	{
		if ( input.size() == 0 ) // after finishing a block reading
		{
			++total_num_block;

			// filter out those less-than-specie_thresh-identical columns
			int valid_len = 0;
			for ( int i=0; i<valid_count.size(); ++i )
			{
				if ( valid_count[i] >= specie_thresh )
					++valid_len;
			}

			if ( valid_len < len_thresh )
				read_on = false;

			// the reading block is finished, call other routines
			if ( read_on )
			{
				// call the evaluation function
				// and print the result
				int eval_start = 0;
				int eval_end = 0;
				int eval_len = 0;
				double eval_percent = 0;
				int eval_num_specie = 0;

				// run the evaluation function, get the evaluated data
				// this action will modify the data set, remove those lines that are not good
				vector<string> original_data_set = data_set;				
				eval_len = eval(data_set, eval_start, eval_num_specie, specie_thresh, len_thresh);
				eval_end = eval_start + eval_len - 1;

				static int here = 0;
				if ( eval_len )
				{
					cout<<eval_len<<" "<<eval_num_specie<<" "<<++here<<endl;
					deout<<eval_len<<" "<<eval_num_specie<<" "<<here<<endl;

					// amend human coordinate
					int delm1 = human_coordinate.find_first_of(':');
					int delm2 = human_coordinate.find_first_of('-');
					string hc_prefix = human_coordinate.substr(0, delm1);
					string hc_start = human_coordinate.substr(delm1 + 1, delm2-delm1-1);
					string hc_end = human_coordinate.substr(delm2+1, human_coordinate.size()-delm2-1);

					int start_offset = 0;
					for ( int i=0; i<eval_start; ++i )
					{
						if ( original_data_set[0][i] != '-' )
							++start_offset;
					}
					int end_offset = 0;
					for ( int i=eval_end; i<original_data_set[0].size(); ++i )
					{
						if ( original_data_set[0][i] != '-' )
							++end_offset;
					}

					
					human_coordinate = hc_prefix + ":" + get_coord(hc_start, start_offset, '+') +"-"+ get_coord(hc_end, end_offset, '-');

					if ( eval_len > best_align_len ) // just for record
						best_align_len = eval_len;

					validate(data_set, eval_start, eval_len, eval_num_specie, human_coordinate);
					++chosen_num_block;
				}
			}

			// reset the parameters
			read_on = true;
			is_human = true;
			human_coordinate.clear();
			data_set.clear();
			valid_count.clear();
			valid_len = 0;
		}

		// the point where we are reading data blocks, change it as preprocessed value changes
		else if ( input[0] == 's' ) // all valid entry start with 's'
		{	
			// first line of the alignment block
			if ( is_human ) // we got the human DNA sequence
			{
				is_human = false;

				// get the human coordinate 
				string local_buffer = input.substr(0, input.find_last_of(' '));
				int tmp = input.find_first_of(' ', 2);
				local_buffer = local_buffer.substr(tmp, local_buffer.size()-tmp);
				tmp = local_buffer.find_first_not_of(' ');
				local_buffer = local_buffer.substr(tmp, local_buffer.size()-tmp);

				char* onum_plus = (char*)malloc( 32 * sizeof(char) );
				char* onum_minus = (char*)malloc( 32 * sizeof(char) );
				char sign;
				int seq_len;

				sscanf( local_buffer.c_str(), "%s %d %c %s",
					onum_plus, &seq_len, &sign, onum_minus );
				string prefix = file_name;
				tmp = prefix.find_first_of('/')+1;
				prefix = prefix.substr(tmp, prefix.size() - tmp); 
	
				if ( sign == '+' )
				{
					string onum_str = onum_plus;
					human_coordinate = prefix.substr(0, prefix.find_first_of('_'))+":" + onum_str + "-" + get_coord(onum_plus, seq_len, '+');
				}
				else if ( sign == '-' )
				{
					string onum_str = onum_minus;
					human_coordinate = prefix.substr(0, prefix.find_first_of('_'))+":" + get_coord(onum_minus, seq_len, '-') + "-" + onum_str ;
				}

				free(onum_plus);
				free(onum_minus);
			}


			// chomp the prefix, leave the aligned sequence
			int tmp = input.find_last_of(' ')+1;
			input = input.substr(tmp, input.size() - tmp);

			if ( ! valid_count.size() ) // initialize the valid length and valid array
			{
				for ( int i=0; i<input.size();++i )
				{
					valid_count.push_back(0);
				}
			}

			// make the input consists only upper case letters
			// also check whether the block is worth further analysis
			for (int i=0; i<input.size(); ++i)
			{
				char ch = input[i];

				if ( ch != '-' )
					++valid_count[i];

				input[i] = toupper(input[i]);
			}

			// update the data set
			data_set.push_back(input);

		}

	}

	printf("\n%s: the number of alignment chosen %d out of %d, best length found %d\n\n", file_name, chosen_num_block, total_num_block, best_align_len);
	
	return chosen_num_block;
}


