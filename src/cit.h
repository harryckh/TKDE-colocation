
#ifndef CIT_H
#define	CIT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <limits.h>

#include <vector>
#include <unordered_set>

#include "bit_operator.h"
#include "data_struct.h"
//#include "data_struct2.h"

//The structure of the node in cit.
typedef struct cit_node
{
	struct cit_node*	parent;
	struct cit_node*	child;

	int 				child_num;

	obj_t*				obj_v;

	struct cit_node*	next;

	bool				flag; //true if its in a row instance

}	cit_node_t;

//The structure of a cit.
typedef struct cit
{
	cit_node_t* root;
	
	//Problem specific information.
	int			node_n;

}	cit_t;


cit_t* cit_ini( );

void cit_release_sub( cit_node_t* x);

void cit_release( cit_t* T);

cit_node_t* alloc_cit_node() ;

cit_node_t* cit_insert(cit_t *T, cit_node_t *x, obj_t *obj_v) ;

void print_tree(cit_node_t *x);


void print_ancestor_sub(cit_node_t *cit_node_v);

void print_ancestor(cit_node_t* cit_node_v);

bool check_tree_sub(cit_node_t *cit_node_v, fsi_t *fsi_v, int depth) ;

bool check_tree(cit_t *T, fsi_t *fsi_v);

#endif






