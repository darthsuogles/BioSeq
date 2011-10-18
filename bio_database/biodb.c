#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "biodb.h"

using namespace std;

/**
 *\brief Environment initialization
 */
/// initialize the database by loading the configure files
/// If database doest not exist, we will create a new one
int SeqDB::init(const char* config_file)  
{
  static int is_open = false;
  if ( is_open )
    return 0;

  if ( config_file != NULL )
    {
      // TODO: read the config file
    }
  
  u_int32_t env_flags =     
    DB_CREATE | /* create the st_g_enviornment variable if it does not exist */
    DB_INIT_MPOOL; /* initialize the in-memory cache */

  u_int32_t db_flags = DB_CREATE;
  int ret;  

  // create an st_g_environment object and initialize it for error reporting
  ret = db_env_create(&_env, 0);
  if ( ret != 0 )
    {
      fprintf(stderr, "Error: can't create environment object\n%s", db_strerror(ret));
      return ret;
    }
    
  // open an st_g_environment 
  // once it is opened, an database could be opened inside it
  ret = _env->open(_env,  /* st_g_env handle */
		   "prokeryotic_proteins", /* st_g_env home directory */
		   env_flags, /* st_g_env open flags */
		   0); /* file mode (0 for default) */
  if ( ret != 0 )
    {
      fprintf(stderr, "Error: cannot open environment\n%s", strerror(ret));
      return ret;
    }  

  // with st_g_environment, support multi-threading and multi-database
  ret = db_create(&_sequence_db, _env, 0);
  if ( ret != 0 )
    {
      //fprintf(err_file_ptr, "%s, %s", program_name, db_strerror(ret));
      fprintf(stderr, "Error: can't create database\n%s", db_strerror(ret));
      return ret;
    }

  /* ret = db_create(&_meta_db, _env, 0); */
  /* if ( ret != 0 ) */
  /*   { */
  /*     //fprintf(err_file_ptr, "%s, %s", program_name, db_strerror(ret)); */
  /*     fprintf(stderr, "Error: can't create database\n%s", db_strerror(ret)); */
  /*     return ret; */
  /*   } */

  // create database if the one we specify does not exist
  ret = _sequence_db->open(_sequence_db, /* DB structure pointer */
			   NULL, /* Transaction pointer */
			   "sequence.db", /* Data file */
			   NULL, /* Optional database logical name */
			   DB_BTREE, /* Database access method */
			   db_flags, /* Open flags */
			   0); /* File mode */ 
  if ( ret != 0 )
    {
      //fprintf(err_file_ptr, "%s, %s", program_name, db_stderror(ret));
      fprintf(stderr, "Error: can't open sequence database\n%s", db_strerror(ret));
      return ret;
    }


  /* // create database if the one we specify does not exist */
  /* ret = _meta_db->open(_meta_db, /\* DB structure pointer *\/ */
  /* 		       NULL, /\* Transaction pointer *\/ */
  /* 		       "meta.db", /\* Data file *\/ */
  /* 		       NULL, /\* Optional database logical name *\/ */
  /* 		       DB_BTREE, /\* Database access method *\/ */
  /* 		       db_flags, /\* Open flags *\/ */
  /* 		       0); /\* File mode *\/  */
  /* if ( ret != 0 ) */
  /*   { */
  /*     //fprintf(err_file_ptr, "%s, %s", program_name, db_stderror(ret)); */
  /*     fprintf(stderr, "Error: can't open meta database\n%s", db_strerror(ret)); */
  /*     return ret; */
  /*   }   */

  is_open = true;
  return 0;
}


/**
 *\brief Close the database
 */
void SeqDB::shutdown()
{

  if ( _sequence_db != NULL )
    _sequence_db->close(_sequence_db, 0);

  /* if ( _meta_db != NULL ) */
  /*   _meta_db->close(_meta_db, 0); */

  if ( _env != NULL )   
    _env->close(_env, 0); 
  _env = NULL;
}


void SeqDB::protein_parser( char* file_name )
{
  ifstream fin;
  fin.open(file_name);

  //vector<string> data_set;
  //vector<string> header;

  string protein;
  string metastr;
  bool first_time = true;

  string input;
  while( getline(fin, input) )
    {
      if ( !input.size() )
	continue;
		
      if ( input[0] == '>' )
	{
	  /* //header.push_back(input); */
	  /* if ( protein.size() ) */
	  /*   data_set.push_back(protein); */
	  /* protein.clear(); */
	  metastr = input;
	  if ( first_time )
	    {
	      first_time = false;
	      continue;
	    }
	  
	  // put the sequence into the database
	  DBT meta, seq;;
	  memset(&meta, 0, sizeof(DBT));       
	  memset(&seq, 0, sizeof(DBT));
	  meta.data = (void*) metastr.c_str();
	  meta.size = metastr.size() * sizeof(char);
	  seq.data = (void*) protein.c_str();
	  seq.size = protein.size() * sizeof(char);

	  _sequence_db->put(_sequence_db, NULL, &meta, &seq, DB_NOOVERWRITE);	  
	  _num_byte += seq.size;
	  
	  metastr.clear();
	  protein.clear();
	}

      else
	{
	  protein += input;
	}
    }
	
  // get the last one in
  if ( protein.size() )
    {
      // put the sequence into the database
      DBT meta, seq;;
      memset(&meta, 0, sizeof(DBT));       
      memset(&seq, 0, sizeof(DBT));
      meta.data = (void*)metastr.c_str();
      meta.size = metastr.size() * sizeof(char);
      seq.data = (void*)protein.c_str();
      seq.size = protein.size() * sizeof(char);

      _sequence_db->put(_sequence_db, NULL, &meta, &seq, DB_NOOVERWRITE);	  
      _num_byte += seq.size;
    }
  
}


/// remove the leading and trailing zeros of an input string
string chomp(string& str)
{
  if ( str.size() == 0 )
    return str;

  int begin_idx = str.find_first_not_of(" \t\r\n");
  if ( begin_idx == -1 )
    return "";
  int end_idx = str.find_last_not_of(" \t\r\n");
  return str.substr(begin_idx, end_idx-begin_idx+1);    
}


/* /\** */
/*  *\brief Parse multiple multiple alignment structure */
/*  *\/ */
/* void SeqDB::msa_parser( char* file_name, vector<string> &data_set, int num_species ) */
/* { */
/*   ifstream fin; */
/*   fin.open(file_name); */

/*   // current line of msa */
/*   int curr_line = 0; */

/*   // initialize data_set */
/*   data_set.clear(); */
/*   for (int i=0; i<num_species; ++i) */
/*     { */
/*       string tmp; */
/*       data_set.push_back(tmp); */
/*     } */

/*   string input; */
/*   while( getline(fin, input) ) */
/*     { */
/*       chomp(input); */
/*       if ( !input.size() || input[0] == '#' ) */
/* 	{ */
/* 	  // set the next line to zero */
/* 	  if ( curr_line != num_species && curr_line != 0 ) */
/* 	    cerr<<"line count error, "<<curr_line<<endl; */
/* 	  curr_line = 0; */
/* 	  continue; */
/* 	} */

/*       else */
/* 	{ */
/* 	  // ignore the first element in the line */
/* 	  // coz it is in capital letter and we don't like it */
/* 	  bool is_update = false; */
/* 	  for (int i=1; i<input.size(); ++i) */
/* 	    { */
/* 	      char ch = input[i]; */

/* 	      // ignore lower cases and underscore */
/* 	      if ( ch == ' '|| */
/* 		   ch == '_' || */
/* 		   ch == '.' || */
/* 		   ch == ':' || */
/* 		   isdigit(ch) || */
/* 		   islower(ch) ) */
/* 		continue; */
	      
/* 	      is_update = true; */
/* 	      data_set[curr_line].push_back(ch); */
/* 	    } */
/* 	  if ( !is_update ) */
/* 	    ++curr_line; */
/* 	} */

/*     } */
/* } */


int main()
{
  // open the database named mein_database.db
  // if the database does not exist, create it
  /* init(); */
  
  /* DB *db; */
  /* open(&db, "mein_database.db");   */
  /* close(db); */

  /* shutdown(); */

  SeqDB db;
  db.init();
  db.protein_parser("../data/NC_004347.faa");
  db.protein_parser("../data/NC_004349.faa");
  db.protein_parser("../data/NC_010622.faa");
  db.protein_parser("../data/NC_010623.faa");
  db.protein_parser("../data/NC_010625.faa");
  db.protein_parser("../data/NC_010627.faa");  

  db.shutdown();
  return 0;
}
