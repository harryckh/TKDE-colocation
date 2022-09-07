
/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */

/*
 *	The implementation of the binary search tree follows pseudo-code 
 *	in the textbook "Introduction to Algorithms".
 */

#ifndef BST_H
#define	BST_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <limits.h>

#include "bit_operator.h"
#include "data_struct.h"
//#include "data_struct2.h"

//The structure of the node in bst.
typedef struct bst_node
{
	struct bst_node*	p;
	struct bst_node*	left;
	struct bst_node*	right;

	//IF specific:
	KEY_TYPE			key;

	BIT_TYPE			p_list;		//in the IR-tree.
	obj_set_t*			p_list_obj;	//a normal one.
	double 				totalWeight;

	//ConstructFeasibleSet-Exact specific:
	//KEY2_TYPE			key2;
	obj_t*				obj_v1;
	obj_t*				obj_v2;
    
	//Cao-Exact specific:
	//cns_t*	cns_v;
	void*				cns_v;

	//data_utility specific:
	KEY_TYPE			key_id;
	
    ///SKECa+
    int freq;
    
}	bst_node_t;

//The structure of a bst.
typedef struct bst
{
	bst_node_t* root;
	
	//Problem specific information.
	int			node_n;
	KEY_TYPE	min;
	KEY_TYPE	max;

}	bst_t;

//The structure for storing a list of bst_nodes.
typedef struct bst_node_list
{
	bst_node_t*				bst_node_v;
	struct bst_node_list*	next;
}	bst_node_list_t;


bst_t* bst_ini( );

void bst_release_sub( bst_node_t* x);

void bst_release( bst_t* T);

void bst_insert( bst_t* T, bst_node_t* z);

void bst_transplant( bst_t* T, bst_node_t* u, bst_node_t* v);

bst_node_t* bst_get_min( bst_node_t* x);

bst_node_t* bst_get_max( bst_node_t* x);

void bst_delete( bst_t* T, bst_node_t* x);

void bst_update( bst_t* T, bst_node_t* x);

void in_order_walk( bst_node_t* x);

bst_node_t* bst_successor( bst_node_t* x);

bst_node_t* bst_predecessor( bst_node_t* x);

void print_bst( bst_t* T);

void test_bst( );

//Updated from 27 August, 2012 and on.
bool in_order_sub( bst_node_t* &x, int &tag);

void in_order_walk_non_recur( bst_node_t* x);

bool get_next_in_order( bst_node_t* &x, int &tag);

bst_node_t* bst_search( bst_t* bst_v, KEY_TYPE key);

void bst_trim( bst_t* bst_v, KEY_TYPE u_bound);

//Augmentating the bst to be a b_heap.
//bst_noded_t* bst_b_h_get_top( bst_t* bst_v);

int bst_search_range( bst_t* bst_v, KEY_TYPE low, KEY_TYPE high, bst_node_t* &low_n, bst_node_t* &high_n);

int bst_search_percentile_range( bst_t* bst_v, int low, int high, bst_node_t* &low_n, bst_node_t* &high_n);


//For data_utility.h only.
bst_node_t* bst_locate( bst_t* bst_v, KEY_TYPE key_id);


#endif
