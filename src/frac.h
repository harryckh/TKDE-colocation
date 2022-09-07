/*
 *	Author: Harry
 *	Email: khchanak@cse.ust.hk
 */

#ifndef frac_h
#define frac_h

#include <stdio.h>
#include "data_utility.h"
#include "irtree.h"
#include "costenum.h"
#include "guo_mck.h"
#include "data_struct_operation.h"
#include "BronKerbosch.h"

//#include <bits/stdc++.h>
#include <unordered_set>

#include <vector>
#include <set>
#include <list>
#include <algorithm>

extern bool debug_mode;
extern IRTree_t IRTree_v;
extern bst_t* IF_v;

extern int cost_tag;
extern double fea_highest_freq;
extern float dist_thr;
extern float min_sup;

extern vector<list<int>> adjacency_list;


/*------------------------ Apriori ------------------------ */

fsi_set_t** apriori(int alg_opt, int numOfObj, int numOfFea);

fsi_set_t* const_L1_apriori();

FEA_TYPE join_check( fsi_t* fsi_v1, fsi_t* fsi_v2);

bool all_subesets_exist(fsi_set_t* fsi_set_v, fsi_t* fsi_v);

/*------------------------ Support Computation ------------------------ */

B_KEY_TYPE comp_support(int alg_opt, fsi_t* fsi_v, int numOfObj, unordered_set<int>* RI1, unordered_set<int>* RI2);

B_KEY_TYPE comp_support_partitioning(int alg_opt, fsi_t* fsi_v, int numOfObj);

B_KEY_TYPE comp_support_construction(int alg_opt, fsi_t* fsi_v, int numOfObj);

B_KEY_TYPE comp_support_enum(int alg_opt, fsi_t* fsi_v, int numOfObj);

void enum_sub(bst_t* IF_v, obj_set_t* S_0, obj_t* o, int& cnt);

B_KEY_TYPE comp_support_fraction(int alg_opt,  fsi_t* fsi_v, int numOfObj, unordered_set<int>* RI1, unordered_set<int>* RI2);

/*------------------------ RI ------------------------ */

int filters(int alg_opt, fsi_t *fsi_v, obj_t *obj_v, int i,
		unordered_set<int> *RI, unordered_set<int> *RIno,
		unordered_set<int> *RI1, unordered_set<int> *RI2) ;
bool check_row_instance(int alg_opt,   fsi_t* fsi_v, obj_t* obj_v, unordered_set<int>* RI, unordered_set<int>* RIno);

bool combinatorial( fsi_t* fsi_v, obj_t* obj_v, obj_set_t* S2, unordered_set<int>* RI, unordered_set<int>* RIno);

obj_set_t* combinatorial_sub( bst_t* IF_v, obj_set_t* S_0, obj_t* o, B_KEY_TYPE d);

bool bst_check_plist_obj_n(bst_node_t* bst_node_v);

bool bst_check_plist(bst_t* bst_v, fsi_t* fsi_v, FEA_TYPE fea);

bool dia( fsi_t* fsi_v, obj_t* obj_v, obj_set_t* S2, unordered_set<int>* RI, unordered_set<int>* RIno);

bool mck(fsi_t* fsi_v, obj_t* obj_v, obj_set_t* S2, unordered_set<int>* RI, unordered_set<int>* RIno);

// Method 4

bool filter_and_verify(int alg_opt2,  fsi_t* fsi_v, obj_t* obj_v, unordered_set<int>* RI, unordered_set<int>* RIno);

bool check_Nof_feasibility(fsi_t* fsi_v, obj_t* obj_v);

bool check_Nof2_feasibility(fsi_t* fsi_v, obj_t* obj_v);

int pruning_NNset( query_t* q, obj_t* obj_v);

/*------------------------ FractionScore ------------------------ */

B_KEY_TYPE min_frac_receive(fsi_t* fsi_v, obj_t* obj_v);

void precomputation(data_t* data_v, KEY_TYPE dist_thr);

bool fea_freq_check(FEA_TYPE fea) ;

bst_t* const_IF( obj_set_t* obj_set_v, psi_t* psi_v, unordered_set<int>* RIno);

/*------------------------ Maximal ------------------------ */


fsi_set_t** maximal(int alg_opt, int numOfObj, int numOfFea);

fsi_set_t* const_L_k_apriori(int alg_opt, int numOfObj, fsi_set_t *L, int i,
		B_KEY_TYPE &sup_avg) ;

vector<list<int>> genGraph(int numOfFea, fsi_set_t *fsi_set_v);

fsi_set_t* candidate_maximal_patterns(set<vector<int>> cliques, int& maxSize);

fsi_set_t** distribute_by_size(fsi_set_t* fsi_set_v, int maxSize);

void find_maximal(int minSize, int maxSize, int numOfObj, fsi_set_t **cmp,
		fsi_set_t **result);

fsi_t* contain_fsi(fsi_set_t *fsi_set_v, fsi_t *fsi_v);
bool is_subset_fsi(fsi_t *fsi_v2, fsi_t *fsi_v);
bool checked_subset_fsi(fsi_set_t *fsi_set_checked, fsi_t *fsi_v);
bool checked_subset_fsi(fsi_set_t **fsi_set_v, int maxSize, fsi_t *fsi_v) ;
#endif /* frac_h */
