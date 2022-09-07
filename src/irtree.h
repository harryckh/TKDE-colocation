/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */

/*
 *	IR-tree.h defines the structures for implementing the IR-tree.
 *
 *	1. The "delete" operation of the IR-tree is not mentioned in the Cong et al.'s paper,
 *	since this operation is not used in their algorithms.
 *	In our implementation, this operation has been implemented.
 *
 *	2. Possible ways for implementing the "delete" operatoin:
 *	When an obj/sub-tree is removed, maintain the set of keywords 
 *	that contained in this obj/sub-tree but not contained by other entries of the same node.
 */

#ifndef IRTREE_H
#define IRTREE_H

#include "bit_operator.h"
#include "bst.h"
#include "data_utility.h"
#include "data_struct.h"

//#define IRTree_OR_R_STAR_TREE	1

#define M				32
#define m1				20

#define MAX_LINE_LENG	1027

//#define DATA_FILE		"IRTree_DataFile.txt"		//
#define SAVE_FILE		"IRTree.txt"				//savefile of rtree.


//The node structure of the node of the rtree.
typedef struct node         
{
	int		num;            //The number of items in the node.
	int		loc;            //The location in its parent's node.

	int		level;			//from 'leaf(0)' to 'root', increased by 1 level by level.
	node*	parent;			//pointer to the parent node.
	void**	child;			//pointer array(because the inner-node and leaf-node has different chidren) 
	range** MBRs;			//closing MBRange within 'x' dimension.

	//Extension for IR-tree.
	//BIT_TYPE *keys;
		//bit array is adopted for representing the posting list.
		//adjacent list is another data structure option.
	//k_node_t** key;

	bst_t*	bst_v;
	
}	node_t;

//The info structure of the rtree.
typedef struct	IRTree
{
	int		dim;
	
	int		split_opt;
	
	node_t*	root;			//root node.
	int		height;			//height.
	int		obj_n;			//the number of the objects indexed by the tree.

	int		inner_n;		//the number of inner-nodes.
	int		leaf_n;			//the number of leaf-nodes.

	//bool*	is_first;		//used for the r*-tree.
	//int	height_max;		//used for the r*-tree.

	//IF info.
	//int		key_n;
	
}	IRTree_t;

//
struct objPointer	//typedef 
{
	obj_t*		p;
	objPointer* next;
};

//
struct nodePointer	//typedef 
{
	node_t*			p;
	nodePointer*	next;
};


extern IRTree_t IRTree_v;
	

/*---------------------Calculation Related APIs----------------------*/

//double abs( double var);

void CreateNode( node_t* &pnode);

void ReleaseNode( node_t* pnode, int tag);

int IsOverlapped( range* MBR1, range* MBR2);

bool IsSame( range* MBR1, range* MBR2);

double GetArea( node_t* pnode);	

double GetArea_2( range* MBR1, range* MBR2);

double GetIncArea( node_t* pnode, range* MBR);

void NeedLeastArea( node_t* pnode, range* obj, int& loc);

void calc_MBR( node_t* pnode, range* &MBR);

range* get_MBR_node( node_t* pnode, int dim);

bool UpdateMRB( node_t* pnode, range* &MBR);
 
/*----------------------Insertion related APIs----------------------*/

void ChooseLeaf( node_t* root, range* MBR, node_t* &leaf_node);

void PickSeeds( node_t* pnode, range* MBR, int& g1, int& g2);

void LinearPickSeeds( node_t* pnode, range* MBR, int& g1, int& g2);

int PickNext( node_t* pnode, node_t* node_1, node_t* node_2, node_t* &node_c, range* MBR, bool* done, int rTag);

int LinearPickNext( node_t* pnode, node_t* node_1, node_t* node_2, node_t* &node_c, bool* done, int rTag);

void AssignEntry( node_t* pnode, node_t* node_c, void* obj, range* MBR);

void SplitNode( node_t* &pnode, node_t* &node_2, void* obj, int opt);

//void QuadraticSplit( node_t* &pnode, node_t* &node_2, void* obj);

//void LinearSplit( node_t* &pnode, node_t* &node_2, void* obj);

void AdjustTree( node_t* pnode);

node_t* CreateNewRoot( node_t* node_1, node_t* node_2);

node_t* InsertSub( node_t* pnode, void* obj);

void Insert( obj_t* obj);

/*---------------------Deletion Related APIs------------------------*/

node_t* Search( node_t* root, range* MBR, int& loc);

void ChooseInner( node_t* root, range* MBR, node_t* &inner_node);

void CondenseSub( node_t* pnode, objPointer* h1, nodePointer* h2 );

void CondenseTree( node_t* pnode);

int Delete( obj_t* obj); 

/*---------------------Other Operation APIs------------------------*/

void ini_tree( );

bool print_and_check_tree( int o_tag, const char* tree_file);

void read_node( FILE* i_fp, node_t* node_v);

obj_t* const_obj( node_t* node_v, int loc, int& id_cnt);

void read_tree( char* ir_tree_f);

void build_IRTree( data_t* data_v);

void free_IRTree_sub( );

void free_IRTree( );

//void bp_free_IRTree( );

bool loc_to_key( int key, int loc, node_t* node_v);

bool loc_oriented_check( node_t* node_v);

bool key_to_loc( int key, int loc, node_t* node_v);

bool key_oriented_check( node_t* node_v);

bool check_IF_sub( node_t* node_v);

void check_IF( );

void print_IF( node_t* node_v, FILE* fp, int tag);

void print_MBR( range* MBR, int dim);

void test_IRTree( int o_tag);

/*---------------------IR-tree augmentation APIs------------------------*/

k_node_t* collect_keywords_list( k_node_t* k_head);

k_node_t* collect_keywords_bst( bst_t* bst_v);

k_node_t* collect_keywords_fea(FEA_TYPE fea);

//void release_k_list( k_node_t* k_head);

int add_IF_entry( bst_t* bst_v, KEY_TYPE key, int loc);

void adjust_IF( node_t* c_node, void* obj, int loc);

void const_IF( node_t* node_v);

void delete_loc_IF( bst_t* bst_v, int loc);

void adjust_parent_IF( node_t* p_node, int loc);

data_t*	IRTree_read_data( IRTree_config_t* cfg);

void IRTree_free_data( data_t* data_v);

//Generate the IR-trees of the datasets.
void gen_irtree( );


#endif
