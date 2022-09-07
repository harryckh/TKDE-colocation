
/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */

/*
*	1. The mechanism of "cmp" function and the "void*" type involved there.
*	2. The choice min-heap or max-heap can be determined by changing the function of "cmp".
*	3. The heap array's subscript starts from *ONE*.
*	4. The concept of "heap array" and "object array", which work collaborarily.
*	5. The size of the heap is static.
*		--Consider make it dynamic.
*
*	5. [Updated on 17 Feb. 2012]
*		(a). Variable MIN_HEAP_OPT is variable for specifying min-heap or max-heap.
*		(b). Object type should maintain an element 'key' for the heap structure.
*/


#ifndef B_HEAP_H
#define	B_HEAP_H

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

#include "irtree.h"
//#include "cao_alg.h"
#include "data_struct2.h"

#define		FLOAT				float
//#define		B_KEY_TYPE			double

#define		B_CONFIG_FILE		"b_config.txt"
#define		MAX_FILENAME_LENG	256

#define		MIN_HEAP_OPT			1
#define		ENLARGE_SCALE			1000

#define		INI_HEAP_SIZE			10000

/*
#define		B_INPUT_FILE	"b_input.txt"
#define		B_OUTPUT_FILE	"b_output.txt"
*/

//The structure of the object that indexed by the binary heap.
typedef	struct h_obj
{
	int			loc;				

	//object-specific variables.
	B_KEY_TYPE	key;	//The key of the heap.
	
	node_t*		node_v;	//For NN_key query.
	obj_t*		obj_v;	//For heap_sort.
	
	cns_t*		cns_v;	//For Cao_Exact.

	void*		element;//For Cao_Appro2.
	int			e_tag;	//e_tag = 1: element is a node_t*;
						//e_tag = 2: element is a obj_t*.

	//void*		cns_v;	//For Cao_Exact;
	
    
    
    bool        type;    //for SKEC+ //true = in, false = out
    B_KEY_TYPE  cir_x;
    B_KEY_TYPE  cir_y;
    
} h_obj_t;


//The structure of the node of a heap.
typedef	struct b_heap
{
	int			size;			//the length of the two array.
	int			rear;			//rear points the index of the 'last' element.
	//int			num;			//the number of elements in the queue.
	h_obj_t*	obj_arr;		//the object array.
	
	int*		h_arr;			//the heap array.

} b_heap_t;



//-------------------------------------

//int cmp( const void* obj_v, int fir, int sec);

int cmp_min( int* array, int n1, int n2, h_obj_t* obj_v);

int cmp_max( int* array, int n1, int n2, h_obj_t* obj_v);

//-------------------------------------

b_heap_t* alloc_b_heap( int size);

void release_b_heap( b_heap_t* b_h);

void b_t_heapify( int* array, int cur, h_obj_t* obj_v);

void t_b_heapify( int* array, int cur, int rear, h_obj_t* obj_v);

void b_h_insert( b_heap_t* &b_h, int n);

int b_h_get_top( b_heap_t* b_h);

void b_h_update_key( b_heap_t* b_h, int k);

//-------------------------------------

void print_b_heap( FILE* o_fp, b_heap_t* b_h);

void test_b_heap( );

//-------------------------------------
//Updated form 31 August and on.
void enlarge_b_heap( b_heap_t* &b_h);

bool b_h_is_empty( b_heap_t* b_h);

#endif
