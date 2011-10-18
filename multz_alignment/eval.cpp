#include "eval.h"

using namespace std;

const int BIT_SIZE = 44;

extern ofstream deout;

int get_ucelem(int arr[], int size, int &start)
{
	// calculate prefix sum
	int *prefix_sum = (int*)malloc(size*sizeof(int));
	for (int i=0; i<size; ++i)
		prefix_sum[i] = arr[i];
	for (int i=1; i<size; ++i)
	{
		prefix_sum[i] += prefix_sum[i-1];
	}

	// calculate prefix minimum list and suffix maximum list
	int *prefix_min = (int*)malloc(size*sizeof(int));
	int *suffix_max = (int*)malloc(size*sizeof(int));
	prefix_min[0] = prefix_sum[0];
	suffix_max[size-1] = prefix_sum[size-1];

	int min = prefix_min[0];
	for (int i=1; i<size; ++i)
	{
		if ( prefix_sum[i]<min )
			min = prefix_sum[i];
		prefix_min[i] = prefix_sum[i];
	}
	int max = suffix_max[size-1];
	for (int i=size-2; i>=0; --i)
	{
		if ( prefix_sum[i] > max )
			max = prefix_sum[i];
		suffix_max[i] = max;
	}

	// calculate the best distance
	int best_dist = 0;
	int best_i = 0;
	int best_j = 0;
	for ( int i =0; i<size; ++i )
		for (int j = size-1; j>0; --j)
		{
			if ( prefix_min[i] <= suffix_max[j] )
			{
				int dist = j - i + 1 ;
				if ( dist > best_dist )
				{
					best_dist = dist;
					best_i = i;
					best_j = j;
				}
			}

		}

		start = best_i;
	
		free(prefix_sum);
		free(prefix_min);
		free(suffix_max);

		return best_dist;
}

// delete all-gap columns
void trim( vector<string> &data_set )
{
	int size = data_set[0].size();
	int num_row = data_set.size();
	
	vector<string> noval_set;

	vector<int> gap_count;
	for ( int i=0; i<size; ++i )
		gap_count.push_back(0);

	for ( int i=0; i<num_row; ++i )
	{
		for ( int j=0; j<size; ++j )
		{
			if ( data_set[i][j] == '-' )
				++gap_count[j];
		}
	}
	
	for ( int i=0; i<num_row; ++i )
	{
		string tmp;
		noval_set.push_back(tmp);
		for ( int j=0; j<size; ++j )
		{
			if ( gap_count[j] < num_row  )
				noval_set[i].push_back( data_set[i][j] );
		}
	}

	data_set = noval_set;
}

// take the data set, calculate the valid subset that satisfy the specie and length constriants
int eval(vector<string> &data_set, int &start, int &specie, int specie_thresh, int len_thresh)
{
	int num_specie = data_set.size();

	// if the number of species is less than the threshold
	// simply drop this answer
	if ( num_specie < specie_thresh )
		return 0;

	// trim the data set
	trim(data_set);
	int size = (data_set[0]).size();

	// data related
	int *dominator = (int*)malloc(5*size*sizeof(int));
	char *dom_letter = (char*)malloc(size*sizeof(char));
	memset(dominator, 0, 5*size*sizeof(int));

	// find the dominating letter
	for (vector<string>::iterator iter=data_set.begin();
		iter!=data_set.end(); ++iter)
	{
		for (int i=0; i<size; ++i)
		{
			char ch = (*iter)[i];
			switch (ch)
			{
			case 'A':
				++dominator[i];
				break;
			case 'T':
				++dominator[i+size];
				break;
			case 'C':
				++dominator[i+size*2];
				break;
			case 'G':
				++dominator[i+size*3];
				break;
			case '-':
				++dominator[i+size*4];
				break;
			default:
				break;
			}
		}
	}

	for (int i=0; i<size; ++i)
	{
		int a_count = dominator[i];
		int t_count = dominator[i + size];
		int c_count = dominator[i + size*2];
		int g_count = dominator[i + size*3];
		int m_count = dominator[i + size*4];

		char ch;
		if (a_count>=specie_thresh)
		{
			ch = 'A';
		}
		else if (t_count>=specie_thresh)
		{
			ch = 'T';
		}
		else if (c_count>=specie_thresh)
		{
			ch = 'C';
		}
		else if (g_count>=specie_thresh)
		{
			ch = 'G';
		}
		else if (m_count>=specie_thresh)
		{
			ch = '-';
		}
		else
		{
			ch = 'n';
		}
		dom_letter[i] = ch;
	}
	free(dominator);


	//2. convert the data set to bitmap
	//	every column has a bitmap, some space for improvement
	vector< bitset<BIT_SIZE> >  bitmap;
	for (int i=0; i<size; ++i)
	{
		bitset<BIT_SIZE> tmp;
		bitmap.push_back(tmp);
	}
	for (int j = 0; j<num_specie; ++j)
	{	
		for (int i=0; i<size; ++i)
		{
			char ch = dom_letter[i];
			if ( ch!='-' && ch!='n' && data_set[j][i] == ch )
				bitmap[i][j] = 1;
			else
				bitmap[i][j] = 0;
		}
	}


	// run the linear algorithm that maximizes the length
	int *arr = (int*)malloc(size*sizeof(int));
	int pos_count = 0;
	for (int i=0; i<size; ++i)
	{
		// if the column has dominate nucleotide
		if ( bitmap[i].count() >= specie_thresh )
		{	
			arr[i] = 1;
			++pos_count;
		}
		else
			arr[i] = -4;
	}
	free(dom_letter);	

	// get the first level answer from the linear algorithm
	int ucelem_start = 0;
	int ucelem_len = get_ucelem(arr, size, ucelem_start);

	// if the maximun possible length is smaller than our expectation
	// we simply don't consider it anymore
	if ( ucelem_len < len_thresh )
	{
		start = 0;
		specie = 0;
		return 0;
	}

	// use the bitmap to see the common species
	// only those within this range can be selected as candidates
	typedef vector< bitset<BIT_SIZE> >::iterator bit_iter;
	bitset<BIT_SIZE> ANDy;
	bitset<BIT_SIZE> ORy(0);
	ANDy.set();
	int valid_col_count = 0;
	for ( int i=ucelem_start; i<ucelem_start+ucelem_len; ++i )
	{
		if ( bitmap[i].count() < specie_thresh )
			continue;
		++valid_col_count;
		ANDy &= (bitmap[i]);
		ORy |= (bitmap[i]);
	}
	bitset<BIT_SIZE> XORy = ANDy ^ ORy;
	int and_count = ANDy.count();
	int or_count = ORy.count();
	int xor_count = XORy.count();

	// if current species 
	if ( and_count >= specie_thresh )
	{
		specie = and_count;
		start = ucelem_start;
		return ucelem_len;
	}

	else // remove some species in the alignment and find if there is a good region
	{
		// select columns with all the species presented here
		int* selector = (int*)malloc(size*sizeof(int));
		for ( int i=0; i<size; ++i )
		{
			// impose a stronger constraint, consider all the species
			if (  bitmap[i].count() == num_specie )
				selector[i] = 1;
			else
				selector[i] = -4;
		}
		int local_ucelem_start = 0;
		int local_ucelem_len = get_ucelem(selector, size, local_ucelem_start); 
		free(selector);
		//cout<<local_ucelem_len<<endl;

		// if we got a good enough alignment with these species
		if ( local_ucelem_len >= len_thresh )
		{
			start = local_ucelem_start;
			specie = num_specie;
			return local_ucelem_len;
		}

		// if the alignment is not good enough
		// remove some species from the current data set and try again
		else 
		{
			// greedy
			if ( or_count - specie_thresh )
			{
				// priority queue
				map<float, int> specie_plist;

				// find the specie with smallest occurrance
				for ( int i =0; i<BIT_SIZE; ++i )
				{
					if ( XORy[i] )
					{
						float local_count = 0;
						for ( int j =0; j<bitmap.size(); ++j )
						{
							local_count += bitmap[j][i];
						}
						while ( specie_plist.count(local_count) )
							local_count += 0.07;							
						specie_plist[local_count] = i;
					}
				}

				int best_len = -1;
				int best_start = 0;
				int best_specie = 0;
				vector<string> best_set;
				for ( map<float, int>::iterator iter = specie_plist.begin(); (iter!=specie_plist.end()) || (ucelem_len = 0); ++iter )
				{
					// remove that specie from our data set		
					vector<string> local_set = data_set;
					int local_specie = 0;
					int local_start = 0;
					local_set.erase( local_set.begin() + iter->second );

					// run the algorithm in the new data set
					int local_len = eval(local_set, local_start, local_specie, specie_thresh , len_thresh);
					
					if  ( local_len >= len_thresh )
					{
						data_set = local_set;
						start = local_start;
						specie = local_specie;
						return local_len;
					}
				}
				
				start = 0;
				specie = 0;
				return 0;
			}

			// if there is no more species to remove, we failed
			else
			{
				start = 0;
				specie = 0;
				return 0;
			}
		}

	}
}
