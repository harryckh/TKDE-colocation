/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */

#ifndef DATA_STRUCT_H
#define DATA_STRUCT_H

#include <unordered_map>
//#include "bst.h"
#include <vector>
#include <unordered_set>

//#define	KEY_TYPE	float
#define	KEY_TYPE	double
#define B_KEY_TYPE double
//#define KEY2_TYPE	double

#define FEA_TYPE    unsigned int

//----------------

//The structure for storing a set of features.
typedef struct fsi {
	int fea_n;
	FEA_TYPE *feaset; //array

	B_KEY_TYPE sup;

	struct fsi *next;   //use when fsi_t in fsi_set_t

	std::vector<struct obj_set*> obj_set_list_v; //used in joinless

	std::unordered_set<int> *RI = NULL;

	std::unordered_set<int> *parentRI = NULL; //for maximal

} fsi_t;

//The structure for storing a set of feature sets.
typedef struct fsi_set {
	int fsi_n;
	fsi_t *head;

	struct fsi_set *next;
} fsi_set_t;

/*
 //The structure for storing a set of sets of feature sets.
 typedef struct fsi_set_set
 {
 int             fsi_set_n;
 fsi_set_t*      head;

 }	fsi_set_list_t;
 */

//----------------

//From costenum.h
//Inverted File related structures.
typedef struct k_node {
	KEY_TYPE key;
	struct k_node *next;

	int freq;

} k_node_t;

//The structure for storing a set of keywords.
typedef struct psi {
	int key_n;
	k_node_t *k_head;

} psi_t;

//From irtree.h
//The interval structure of one dimension.
typedef struct {
	float min;
	float max;

} range;

//The structure of the object.
typedef struct obj {
	int id;
	range *MBR;

	FEA_TYPE fea;
	double weight;

	//  std::unordered_map<FEA_TYPE, int>* N_o_f; //storing feature count pair

	struct bst *inverted_list;          //storing the inverted list in Disk(o,d)

	//  std::unordered_map<FEA_TYPE, int>* N_o_f2; //storing feature count pair for D(o,d/2)

	//   std::unordered_map<FEA_TYPE, float>* frac_f; //storing fraction received

//    int* N_o_f;
	short int *N_o_f2 = NULL;
	float *frac_f = NULL; //fraction received from each label

} obj_t;

//The node structure for storing an object pointer.
typedef struct obj_node {
	obj_t *obj_v;
	struct obj_node *next;

	//for const_k_NN_key only.
	B_KEY_TYPE dist;

	//for SKECa+
	B_KEY_TYPE maxInvalidRange;
	bool type; //true=in, false = out
	B_KEY_TYPE angle;
	struct obj_set *range = NULL;

} obj_node_t;

//The structure for storing a group of objects.
typedef struct obj_set {
	int obj_n;
	obj_node_t *head;

} obj_set_t;

//From costenum.h.

//The structure for storing a location.
typedef struct loc {
	int dim;
	float *coord;
} loc_t;

//Structure for storing a disk.
typedef struct disk {
	loc_t *center;
	B_KEY_TYPE radius;
} disk_t;

//Structure for storing a query.
typedef struct query {
	//loc.
	loc_t *loc_v;

	//doc.
	//int			key_n;
	//k_node_t*	doc_v;

	psi_t *psi_v;
} query_t;


//Triplet structure.
typedef struct tri {
	obj_t *o;
	obj_t *o_1;
	obj_t *o_2;
} tri_t;

//From data_utility.h
//The data structure for storing the statistics.
typedef struct colocation_stat {
	//cost.
	B_KEY_TYPE aver_cost;
	float aver_size;

	//time.
	float pre_time;
	float q_time;
	float irtree_build_time;

	//memory.
	double memory_v;
	double memory_max;

	double tree_memory_v;
	double tree_memory_max;

	//float		data_memory_v;
	//float		data_memory_max;

	//algorithm related.
	float n_1_sum;
	float achi_sum;

	float O_t_size_aver;
	float O_t_size_sum;

	float psi_n_aver;
	float psi_n_sum;

	float O_simp_size;

	//Cao-Exact specific.
	float node_set_n;
	float feasible_set_n;

	//Cao-Appro2
	float n_k;

	//Additional.
	float ratio_min;
	float ratio_max;
	float ratio_aver;
	float ratio_dev;

	//Method 4
	double S0_sum;
	double S1_sum;
	double S2_sum;
	double S3_sum;
	double S4_sum;
	double S5_sum;
	double S6_sum;

	double S31_sum;
	double S41_sum;


	double S0_time;
	double S1_time;
	double S2_time;
	double S3_time;
	double S4_time;
	double S5_time;
	double S6_time;

} colocation_stat_t;

#endif
