/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */

#ifndef DATA_UTILITY_H
#define	DATA_UTILITY_H


#include <unordered_map>
#include "data_struct.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <random>
using namespace std;


#ifndef WIN32
#include<sys/resource.h>
#endif

#define TEXT_COL_MAX			1075
#define MAX_FILENAME_LENG		256			//maximum filename length.
#define	MAX_DESCRIPTION_LENG	256

#define	KEY_RANGE_SIZE			1
#define K_NN_PARAMETER			10


#define	CONFIG_FILE				"IRTree_config.txt"
#define	COLOCATION_CONFIG_FILE	"config.txt"
#define	COLOCATION_STAT_FILE	"stat.txt"
#define	COLOCATION_RESULT_FILE	"result.txt"
#define	COLOCATION_RESULT_FILE2	"result_all.txt"
#define COLOCATION_PATTERN_FILE "detail.txt"

#define GEN_IRTREE_CONFIG		"gen_irtree_config.txt"

#define SYN_CONFIG_FILE			"syn_config.txt"



//The structure for storing the frequency of the keywords.
//
typedef struct key_freq
{
	float*		freq;
	KEY_TYPE	key;
}	key_freq_t;


//The structure for storing the objects.
typedef struct 
{
	int			dim;
	int			obj_n;
	obj_t*		obj_v;

	int			key_n;

}	data_t;


//The structure for storing the configuration info.
typedef	struct
{
    //	int		M;
    //	int		m;
    int		dim;
    //int	height_max;
    int		obj_n;
    int		key_n;
    
    int		split_opt;
    
    //char	data_file[ MAX_FILENAME_LENG];
    char	loc_file[ MAX_FILENAME_LENG];
    char	doc_file[ MAX_FILENAME_LENG];
    
    char	tree_file[ MAX_FILENAME_LENG];
    
    //IR-tree augmentation.
    
}	IRTree_config_t;


//The structure for storing the configuration information for the co-location pattern mining problem.
typedef struct colocation_config
{
	//Dataset specific.
	int		dim;
	int		obj_n; //|O|
	int		key_n; // total number of keywords appear in O

	char	loc_file[ MAX_FILENAME_LENG];
	char	doc_file[ MAX_FILENAME_LENG];
    char	tree_file[ MAX_FILENAME_LENG];

	//Co-location mining.
	int		alg_opt; //Method 1-5
    int     cost;
    
    int     query_n;
    
	float    min_sup;
    float    min_conf;
    float    dist_thr;

    int		tree_tag;
    
}	colocation_config_t;

//The structure of the config file for generating synthetic data.
typedef struct syn_config
{
	
	char	loc_file[ MAX_FILENAME_LENG];
	char	doc_file[ MAX_FILENAME_LENG];

    //parameters
    int Ncoloc;    //number of colocation
    
    int lambda1;    //size of each colocation pattern
    int lambda2;   //size of each table instance of a colocation pattern
    
    int seed1;      //random seed for lambda1
    int seed2;      //random seed for lambda2
    
    double D1, D2;
    double d;
    
    double r_noise_fea;   //number of *additional* features for noise / number of original features
    double r_global;      //number of *additional* instances from noise features
    double r_local;       //number of *additional* instances from original features
    
    double m_overlap;
    double m_clump;
	
	int m_clump_type;	//0: fixed, 1: uniform range [1,m_clump]
	
}	syn_config_t;

//
typedef struct key_freq_pair
{
	int		key;
	int		freq;
}	key_freq_pair_t;


extern colocation_stat_t stat_v;

//extern IRTree_t IRTree_v;

#ifndef WIN32

void GetCurTime( rusage* curTime);

void GetTime( struct rusage* timeStart, struct rusage* timeEnd, float* userTime, float* sysTime);

#endif

IRTree_config_t* read_config_irtree( );

colocation_config_t* read_config_colocation( );

void add_keyword_entry( k_node_t* &k_node_v, KEY_TYPE key);

void copy_k_list( k_node_t* k_head1, k_node_t* k_head2);

void print_k_list( k_node_t* k_head, FILE* o_fp);

void print_loc( loc_t* loc_v, FILE* o_fp);

void print_colocation_stat( colocation_config_t* cfg, int cnt);

void alloc_obj( obj_t* obj_v, int dim);

data_t*	read_data_irtree( IRTree_config_t* cfg);

data_t* alloc_data( int num);

data_t*	read_data_colocation( colocation_config_t* cfg);

void release_k_list( k_node_t* k_node_v);

void release_data( data_t* data_v);



//Query & dataset generation utilities.
int rand_i( int min, int max);

float rand_f( float min_v, float max_v);

bool is_old( int* rand_v, int cur, int var);


//Generate synthetic datasets.
float gaussian_f( float mean, float s_var);

range* collect_data_range( data_t* data_v);

void print_data( data_t* data_v, syn_config_t* s_cfg);

void gen_syn_data( syn_config_t* s_cfg);

void batch_gen_syn_data( );

void batch_gen_syn_data2( );

void gen_syn_data2( syn_config_t* s_cfg);

void gen_syn_noise(int num_noise_obj, FILE* loc_fp, FILE* doc_fp, int& objectID, int feaStart, int feaEnd, double rectX, double rectY, double d);

void gen_syn_obj(FILE* loc_fp, FILE* doc_fp, int& objectID, int fea, double rectX, double rectY, double d);

void gen_scalability_data(data_t* data_v);
void gen_scalability_data2(data_t* data_v);
#endif




