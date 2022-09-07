/*
 * yao_maximal.h
 *
 *  Created on: 19 Apr 2022
 *      Author: kai-ho
 */

#ifndef YAO_MAXIMAL_ALG_H_
#define YAO_MAXIMAL_ALG_H_

#include "frac.h"

#include "cit.h"

//pair structure.
typedef struct pair {
	obj_t *o_1;
	obj_t *o_2;
} pair_t;

struct PairHash {
public:
	size_t operator()(const pair_t &str) const {
//        int size = str.o_1->id + str.o_2->id;
		return std::hash<int>()(str.o_1->id);
	}
};
// Custom comparator that compares the string objects by length
struct PairEqual {
public:
	bool operator()(const pair_t &str1, const pair_t &str2) const {
		return (str1.o_1->id == str2.o_1->id && str1.o_2->id == str2.o_2->id);
	}
};

fsi_set_t** maximal_yao(data_t *data_v, int numOfObj, int numOfFea);

unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>> const_ins_table2(
		data_t *data_v);

void find_maximal_tree(int minSize, int maxSize, int numOfObj, fsi_set_t **cmp,
		fsi_set_t **result,
		unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>>* insTable2);

B_KEY_TYPE calc_sup_tree(fsi_t *fsi_cur,
		unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>> *insTable2);

cit_t* const_cit(fsi_t *fsi_v,
		unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>>* insTable2);

string create_key_for_table(FEA_TYPE fea1, FEA_TYPE fea2);

cit_node_t* is_child(cit_node_t *cit_node_v, obj_t *obj_v);

void set_ancestor_flag(cit_node_t *cit_node_v, unordered_set<int> *RI);

void cit_search_sub(cit_node_t *cit_node_v, int i, vector<cit_node_t*> &list,
		obj_t *obj_v, int depth);

vector<cit_node_t*> cit_search(cit_t *T, int i, obj_t *obj_v);

bool check_clique(cit_t *T, cit_node_t *cit_node_v, obj_t *obj_v,
		unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>>* insTable2);

void add_sup_tree_sub(cit_node_t *cit_node_v, int depth, fsi_t *fsi_v,
		double *sups);

B_KEY_TYPE add_sup_tree(cit_t *T, fsi_t *fsi_v);

void set_RI_sub(cit_node_t* cit_node_v,	unordered_set<int> *RI, int num, int depth);


void set_RI(cit_t* T, 	unordered_set<int> *RI,int num);



bool is_instance_sub(cit_node_t *cit_node_v, int num, int depth);

bool is_instance(cit_t *T, int num, int targetDepth, obj_t *obj_v) ;

void test_insTable2(
		unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>> insTable2,
		data_t *data_v);

#endif /* YAO_MAXIMAL_ALG_H_ */
