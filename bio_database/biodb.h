#include <db.h>
#include <vector>
#include <string>

using namespace std;

class SeqDB
{
  /* // both dbs use the same unique query tag for each protein as the key */
  /* DB* sequence_db; */
  /* DB* meta_db; */
  
  /* char *db_home_dir; */
  /* vector<string> seq_db_names; // locate the database in raw format */

public:
  // initialize the database by loading the configure files
  // if the database doest not exist, we will create a new one
  int init(const char* config_file = NULL);
  void shutdown();

  // import new data into the database
  //int import(const char* config_file);
  int open();
  void close();

  /// parse the FASTA type data in .faa format
  void protein_parser( char* file_name );

  /// Read the Multiple Sequence Alignment data structure
  //void msa_parser( char* file_name, vector<string> &data_set, int num_species );

private:

  vector<string> file_list;
  DB_ENV* _env;
  DB* _sequence_db;  
  size_t _num_byte;

  ////////////////////////////////////////////////////////////////////////
  // Helper Functions
  ////////////////////////////////////////////////////////////////////////  
  int open(DB** dbp, const char* data_file);
 
};

