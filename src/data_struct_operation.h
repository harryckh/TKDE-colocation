//
//  data_struct_operation.h
//  colocation
//
//  Created by Harry on 26/10/2018.
//  Copyright © 2018 Harry. All rights reserved.
//

#ifndef data_struct_operation_h
#define data_struct_operation_h

#include "data_utility.h"
#include "costenum.h"
#include "irtree.h"

/* ------------------------- create/add/release  ------------------------- */

fsi_t* alloc_fsi( int key_n);

fsi_set_t* alloc_fsi_set();

//fsi_set_list_t* alloc_fsi_set_list();

void release_fsi( fsi_t* fsi_v);

void release_fsi_set( fsi_set_t* fsi_set_v);

//void release_fsi_set_list( fsi_set_list_t* fsi_set_list_v);

void add_fsi_set_entry( fsi_set_t* fsi_set_v, fsi_t* fsi_v);

void remove_fsi_set_entry(fsi_set_t *fsi_set_v) ;

fsi_t* add_fsi(fsi_t* fsi_v, FEA_TYPE fea);

fsi_t* copy_fsi(fsi_t* fsi_v);

bst_t* const_IF(obj_set_t* obj_set_v);

void release_IF_sub( bst_node_t* x);

void release_IF( bst_t* T);

/* ------------------------- print / checking / convert  ------------------------- */

void print_fsi(fsi_t* fsi_v, FILE* o_fp);

void print_fsi_set(fsi_set_t *fsi_set_v, FILE *o_fp);

void print_fsi_set(fsi_set_t** result, int numOfFea, FILE* o_fp);

void print_maximal_fsi_set(fsi_set_t** result, int numOfFea, FILE* o_fp);

bool is_maximal(fsi_t* fsi_cur, fsi_set_t* fsi_set_v);

bool is_subset(fsi_t* fsi_v1, fsi_t* fsi_v2);

bool has_fea(fsi_t* fsi_v, FEA_TYPE fea);

bool check_obj_set_equal_interest(obj_set_t* v1, obj_set_t* v2);

bool check_feasibility(obj_set_t* O_t, fsi_t* fsi_v);

psi_t* fsi_to_psi(fsi_t* fsi_v, FEA_TYPE feaToExclude);

/* ------------------------- IR-tree  ------------------------- */

void retrieve_sub_tree( node_t* node_v, obj_set_t* &obj_set_v);

void range_query_sub( node_t* node_v, disk_t* disk_v, obj_set_t* &obj_set_v);

obj_set_t* range_query( disk_t* disk_v);

#endif /* data_struct_operation_h */
