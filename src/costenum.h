
/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */

#ifndef COSTENUM_H
#define	COSTENUM_H

#include "data_struct.h"
#include "irtree.h"
#include "b_heap.h"
#include "frac.h"


//#define max(a, b) ((a) >= (b) ? (a) : (b))
//#define min(a, b) ((a) <= (b) ? (a) : (b))


extern	int	cost_tag;
extern  float dist_thr;


//#define B_KEY_TYPE	double

//#define	COSTENUM	1


/*
//The structure of an entry in the IF structure.
typedef struct IF_entry
{
	obj_set_t*	p_list;
	KEY_TYPE	key;
}	IF_entry_t;

//Inverted File (IF) structure.
typedef struct IF
{
	int			key_n;
	bst_t*		bst_v;
	//IF_entry_t*	entry_v;
}	IF_t;
*/

extern IRTree_t IRTree_v;


loc_t* alloc_loc( int dim);

loc_t* get_obj_loc( obj_t* obj_v);

loc_t* copy_loc( loc_t* loc_v);

void release_loc( loc_t* loc_v);

psi_t* alloc_psi( );

void add_psi_entry( psi_t* psi_v, KEY_TYPE key);

void release_psi( psi_t* psi_v);

query_t* alloc_query( );

void print_query( query_t* q, FILE* o_fp);

query_t* read_query( FILE* i_fp);

void release_query( query_t* q);

disk_t* alloc_disk( int dim);

void set_disk( disk_t* disk_v, loc_t* loc_v, B_KEY_TYPE radius);

disk_t* const_disk( loc_t* loc_v, B_KEY_TYPE radius);

void release_disk( disk_t* disk_v);

obj_set_t* alloc_obj_set( );

obj_set_t* copy_obj_set( obj_set_t* obj_set_v);

bool has_same_content_obj_set( obj_set_t* obj_set_v1, obj_set_t* obj_set_v2);

void remove_identical_obj( obj_set_t* obj_set_v);

void print_obj_set( obj_set_t* obj_set_v, FILE* o_fp);

void release_obj_set( obj_set_t* obj_set_v);

//IF_entry_t* alloc_IF_entry( );

//void release_IF_entry( IF_entry_t* IF_entry_v);

//IF_t* alloc_IF( );

//void release_IF( IF_t* IF_v);

B_KEY_TYPE calc_dist_loc( loc_t* loc_v1, loc_t* loc_v2);

B_KEY_TYPE calc_dist_obj( obj_t* obj_v1, obj_t* obj_v2);

bool has_key_obj( obj_t* obj_v, KEY_TYPE key);

bool is_relevant_obj( obj_t* obj_v, query_t* q);

bool is_covered_obj_set( obj_set_t* obj_set_v, query_t* q);

BIT_TYPE has_key_node( node_t* node_v, KEY_TYPE key);

BIT_TYPE is_relevant_node( node_t* node_v, query_t* q);

B_KEY_TYPE calc_minDist( range* MBR, loc_t* loc_v);

B_KEY_TYPE calc_maxDist( range* MBR, loc_t* loc_v);

B_KEY_TYPE calc_minDist_node( node_t* node_v, loc_t* loc_v);

//obj_t* NN_key( loc_t* loc_v, KEY_TYPE key);

bool is_overlap( range* MBR, disk_t* disk_v); 

obj_t* const_NN_key( loc_t* loc_v, KEY_TYPE key, disk_t* disk_v);

obj_set_t* const_k_NN_key( loc_t* loc_v, KEY_TYPE key, disk_t* disk_v, int k);

bool is_enclosed( range* MBR, disk_t* disk_v);

void add_obj_set_entry( obj_t* obj_v, obj_set_t* obj_set_v);

void remove_obj_set_entry( obj_set_t* obj_set_v);

void retrieve_sub_tree( node_t* node_v, obj_set_t* &obj_set_v, query_t* q);

void range_query_sub( node_t* node_v, disk_t* disk_v, obj_set_t* &obj_set_v, query_t* q);

obj_set_t* range_query( disk_t* disk_v, query_t* q);

void refine_region( obj_set_t* obj_set_v, disk_t* disk_v);

psi_t* get_psi_obj( obj_t* obj_v);

psi_t* const_psi( k_node_t* k_head);

k_node_t* key_exclusion( k_node_t* k_head1, k_node_t* k_head2);

k_node_t* key_exclusion( k_node_t* k_head1, KEY_TYPE key);

psi_t* psi_exclusion( psi_t* psi_v1, psi_t* psi_v2);

void psi_exclusion( psi_t* psi_v, FEA_TYPE fea);

bool check_feasibility( obj_set_t* O_t, psi_t* psi);

void obj_filter_range( obj_set_t* &O_t, disk_t* disk_v);

bst_t* const_IF( obj_set_t* obj_set_v, psi_t* psi_v);

obj_set_t* const_obj_set( tri_t* triplet_v);

bool check_dist_constraint( obj_set_t* obj_set_v, obj_t* obj_v, obj_t* o, B_KEY_TYPE d);

bst_node_list_t* update_IF_obj( bst_t* IF_v, obj_t* obj_v);

void release_bst_node_list( bst_node_list_t* bst_node_list_v);

void restore_IF_bst_node_list( bst_t* IF_v, bst_node_list_t* bst_node_list_v);

void combine_obj_set( obj_set_t* obj_set_v1, obj_set_t* obj_set_v2);

obj_set_t* const_feasible_set_sub( bst_t* IF_v, obj_set_t* S_0, obj_t* o, B_KEY_TYPE d);

obj_set_t* const_feasible_set( obj_set_t* O_t, psi_t* psi, obj_t* o, B_KEY_TYPE d);

obj_set_t* AchievabilityCheck( tri_t triplet_v, query_t* q);

B_KEY_TYPE comp_farthest( obj_set_t* obj_set_v, query_t* q);

B_KEY_TYPE comp_diameter( obj_set_t* obj_set_v);

B_KEY_TYPE comp_cost( int cost_tag, obj_set_t* obj_set_v, query_t* q);

obj_set_t* comp_bounds( query_t* q, B_KEY_TYPE &LB, B_KEY_TYPE &UB, disk_t* disk_v);

b_heap_t* heap_sort_obj_set( obj_set_t* obj_set_v, query_t* q);

void obj_exclusion_disk( obj_set_t* obj_set_v, disk_t* disk_v);

void obj_inclusion_disk( obj_set_t* obj_set_v, disk_t* disk_v);

void process_obj_pair( obj_t* o_1, obj_t* o_2, B_KEY_TYPE dist, bst_t* obj_pair);

bst_t* ini_obj_pair( obj_set_t* region, obj_t* o, B_KEY_TYPE d_l, B_KEY_TYPE d_u, query_t* q);

void update_obj_pair( bst_t* obj_pair, obj_set_t* region, obj_t* o, B_KEY_TYPE d_l, B_KEY_TYPE d_u, query_t* q);

obj_set_t* ConstructFeasibleSet_Exact( obj_t* o, query_t* q, bst_t* obj_pair);

obj_set_t* ConstructFeasibleSet_Appro( obj_t* o, query_t* q);

obj_set_t* CostEnum( query_t* q, int s_tag, int prune_tag, obj_set_t* obj_set_v, disk_t* disk_v);

obj_set_t* CostEnum_Exact( query_t* q, int prune_tag);

obj_set_t* CostEnum_Appro( query_t* q);

//void test_CostEnum( int s_tag);



//Query related.
/*
query_t* gen_query( int key_n, range* MBR, data_t* data_v, int low, int high);

void collect_key_freq( data_t* data_v, int* key_freq, int size);

int get_max_key( data_t* data_v);

int compare (const void * a, const void * b);

int compare_key_freq (const void * a, const void * b);

query_t* gen_query2( int key_n, range* MBR, data_t* data_v, int low, int high);

//query_t* gen_query3( int key_n, range* MBR, data_t* data_v, int low, int high);

query_t** gen_query_set( int query_n, int key_n, range* MBR, data_t* data_v, int low, int high);

int	get_key_num( int* freq_v, int size);

query_t** gen_query_set2( int query_n, int key_n, range* MBR, data_t* data_v, int low, int high);

//For gen_query_set3.

loc_t* gen_loc( obj_t* obj_v1, obj_t* obj_v2, KEY_TYPE key1, KEY_TYPE key2);

obj_set_t* retrieve_obj_key( data_t* data_v, KEY_TYPE key);

loc_t* gen_loc( KEY_TYPE key1, KEY_TYPE key2, data_t* data_v);

void gen_loc( query_t* q, data_t* data_v, range* MBR);

query_t** gen_query_set3( int query_n, int key_n, range* MBR, data_t* data_v, int low, int high);
*/
#endif
