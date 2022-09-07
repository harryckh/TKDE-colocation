/*
 *	Author: Harry
 *	Email: khchanak@cse.ust.hk
 */

#ifndef yoo_alg_h
#define yoo_alg_h

#include <stdio.h>
#include "data_utility.h"
#include "irtree.h"
#include "costenum.h"
#include <vector>

fsi_set_t** joinless_mining(data_t* data_v, int numOfObj, int numOfFea);

unordered_map<FEA_TYPE, vector<obj_set_t*>> gen_star_neighborhoods(data_t* data_v);

vector<obj_set_t*> gen_one_star_neighborhood(obj_t* obj_v, obj_set_t* obj_set_v);

B_KEY_TYPE comp_PI(fsi_t* fsi_v, obj_set_t* O);

B_KEY_TYPE comp_sup(fsi_t* fsi_v, obj_set_t* O);

void check_star_instance(fsi_t* fsi_v,  unordered_map<FEA_TYPE, vector<obj_set_t*>> SN);

vector<obj_set_t*> filter_star_instance(obj_t* obj_v, bst_t* inverted_list);

void filter_star_instance_sub( bst_t* IF_v, obj_set_t* S_0, vector<obj_set_t*>* obj_set_list_v );

void filter_clique_instance(fsi_t* fsi_cur, fsi_set_t* fsi_set_v, vector<obj_set_t*>* obj_set_list_v);

bool check_obj_set_equal(obj_set_t* v1, obj_set_t* v2, FEA_TYPE fea);

bool obj_exist(vector<obj_set_t*> obj_set_list_v, obj_t* obj_v);

#endif /* yoo_alg_hpp */
