/*
 *	Author: Harry
 *	Email: khchanak@cse.ust.hk
 */

#include "yoo_alg.h"
#include "irtree.h"
#include "frac.h"
#include "data_utility.h"
#include <iostream>
#include <fstream>
#include <random>
#include "SmallestEnclosingCircle.h"
#include <cstddef>
#include <cmath>
#include "stringgenc.h"
#include "BronKerbosch.h"
#include "yao_maximal_alg.h"

using namespace std;
using std::size_t;
using std::vector;

IRTree_t IRTree_v;
bst_t *IF_v;

colocation_stat_t stat_v;

bool debug_mode = false;

//----------------------------------
//user parameter
float dist_thr;
float min_sup;
int cost_tag; //1= participation based; 2=fraction based
double fea_highest_freq;
double fea_highest_obj_n;
double alpha; //for maximal pattern mining

//------------------------------------
void colocation_patterns_support(int cost);

void colocation();

void build_IF(data_t *data_v);

double get_highest_freq(bst_node_t *bst_node_v);
double get_highest_obj_n(bst_node_t *bst_node_v);
double print_freq(bst_node_t *bst_node_v);

fsi_set_t* read_patterns();

fsi_set_t** find_maximal_fsi_set(fsi_set_t **prev, int numOfFea);
int main(int argc, char *argv[]) {

	//     batch_gen_syn_data2();
	colocation();

	return 0;
}

void colocation() {
	int i;
	colocation_config_t *cfg;
	FILE *r_fp, *r_fp2;
	data_t *data_v;
	fsi_set_t **result;
	fsi_set_t **tmp;

	memset(&stat_v, 0, sizeof(colocation_stat_t));

	//Read the config.
	printf("Reading configuration ...\n");
	cfg = read_config_colocation();

	cost_tag = cfg->cost;
	dist_thr = cfg->dist_thr;
	min_sup = cfg->min_sup;
	//  min_conf = cfg->min_conf;

	//Read the data.
//	printf("Reading data ...\n");
	data_v = read_data_colocation(cfg);


#ifndef WIN32
	float sys_t, usr_t, usr_t_sum = 0;
	struct rusage IR_tree_sta, IR_tree_end;

	GetCurTime(&IR_tree_sta);
#endif

	//Option 1: Build the tree from scratch.
	//Build the IR-tree.
	if (cfg->tree_tag == 0) {
//		printf("Building IR-tree ...\n");
		build_IRTree(data_v);

		print_and_check_tree(1, cfg->tree_file);
		//check_IF( );
	} else {
		//Option 2: Read an existing tree.
//		printf("Reading IR-Tree ...\n");
		read_tree(cfg->tree_file);
	}

#ifndef WIN32
	GetCurTime(&IR_tree_end);
	GetTime(&IR_tree_sta, &IR_tree_end, &stat_v.irtree_build_time, &sys_t);
#endif

	//Build Inverted file
	build_IF(data_v);

	//---

	if (cost_tag == 1) {
		printf("participation based: min_sup:%lf\n", min_sup);
	} else {
		fea_highest_freq = get_highest_freq(IF_v->root);
		fea_highest_obj_n = get_highest_obj_n(IF_v->root);

//	print_freq(IF_v->root);
		printf("Number of fea:%d\tHighest freq:%lf\t highest obj_n:%lf\n",
				IF_v->node_n, fea_highest_freq, fea_highest_obj_n);
		printf("min_sup:%lf\t", min_sup);
		printf("alpha:%.3lf\n", alpha);
	}
	//---
	//    //Get the whole range.
	//    MBR = get_MBR_node( IRTree_v.root, IRTree_v.dim);
	printf("precomputation ...\n");

#ifndef WIN32

	float pre_t = 0;
	struct rusage pre_sta, pre_end;
	GetCurTime(&pre_sta);
#endif

	//Pre-computation of N_o_f
	precomputation(data_v, cfg->dist_thr);

#ifndef WIN32
	GetCurTime(&pre_end);
	GetTime(&pre_sta, &pre_end, &pre_t, &sys_t);
	printf("pre time:%0.5lf\n", pre_t);
	stat_v.pre_time = pre_t;

#endif
	printf("Running apriori (Method %d):\n", cfg->alg_opt);

	for (i = 0; i < cfg->query_n; i++) //i==1
			{

#ifndef WIN32
		struct rusage query_sta, query_end;
		GetCurTime(&query_sta);

#endif


		if (cfg->alg_opt == 5)
			result = joinless_mining(data_v, cfg->obj_n, cfg->key_n);
		else if (cfg->alg_opt == 46)
			result = maximal(cfg->alg_opt, cfg->obj_n, cfg->key_n);
		else if (cfg->alg_opt == 47) {
			tmp = apriori(41, cfg->obj_n, cfg->key_n);

#ifndef WIN32
			struct rusage query_sta1, query_end1;
			GetCurTime(&query_sta1);
#endif
			result = find_maximal_fsi_set(tmp, cfg->key_n);

#ifndef WIN32
			GetCurTime(&query_end1);
			GetTime(&query_sta1, &query_end1, &usr_t, &sys_t);

			printf("find maximal Time: %f\n\n", usr_t);
			//pre-computation time should also be included in the query time

#endif
		} else if (cfg->alg_opt == 48) {
			result = maximal_yao(data_v, cfg->obj_n, cfg->key_n);
		} else
			result = apriori(cfg->alg_opt, cfg->obj_n, cfg->key_n);

#ifndef WIN32
		GetCurTime(&query_end);

		GetTime(&query_sta, &query_end, &usr_t, &sys_t);
		usr_t_sum += usr_t;

		printf("qTime: %f\n", usr_t);
		printf("Total time: %f\n\n", usr_t + pre_t);
		//pre-computation time should also be included in the query time
		stat_v.q_time = usr_t_sum / (i + 1);

#endif

		//print size for each L_i
		int cnt = 0;
		for (int i = 0; i < cfg->key_n; i++) {
			if (result[i] != NULL && result[i]->fsi_n > 0) {
				printf("L_%d fsi_n:%d\n", i + 1, result[i]->fsi_n);
				cnt += result[i]->fsi_n;
			}
		}
		printf("total fsi_n:%d\n", cnt);

		//Print the query results.
		if (i == 0) {
			if ((r_fp = fopen(COLOCATION_RESULT_FILE, "w")) == NULL) {
				fprintf(stderr, "Cannot open the colocation_result file.\n");
				exit(0);
			}

			//Print the query result.
			print_fsi_set(result, cfg->key_n, r_fp);
			fclose(r_fp);

			//also print all patterns
			if (cfg->alg_opt == 47) {
				if ((r_fp2 = fopen(COLOCATION_RESULT_FILE2, "w")) == NULL) {
					fprintf(stderr,
							"Cannot open the colocation_result2 file.\n");
					exit(0);
				}

				print_fsi_set(tmp, cfg->key_n, r_fp2);

				fclose(r_fp2);
			}

		}
		//release the result memory
		for (int k = 0; k < cfg->key_n; k++)
			if (result[k] != NULL)
				release_fsi_set(result[k]);

		print_colocation_stat(cfg, i + 1);

		//reset some statistcs for the next query
		stat_v.S0_sum = 0.0;
		stat_v.S1_sum = 0.0;
		stat_v.S2_sum = 0.0;
		stat_v.S3_sum = 0.0;
		stat_v.S5_sum = 0.0;
		stat_v.S6_sum = 0.0;

		stat_v.S0_time = 0.0;
		stat_v.S1_time = 0.0;
		stat_v.S2_time = 0.0;
		stat_v.S3_time = 0.0;
		stat_v.S4_time = 0.0;
		stat_v.S5_time = 0.0;
		stat_v.S6_time = 0.0;

	}

	free(cfg);
}

/*
 * Construct the inverted file @IF_v based on the data
 */
void build_IF(data_t *data_v) {

	bst_node_t *bst_node_v;

	IF_v = bst_ini();

	//Insert all the objects to construct the IF
	for (int i = 0; i < data_v->obj_n; i++) {
		//        if(i%100==0)
		{

			bst_node_v = bst_search(IF_v, data_v->obj_v[i].fea);

			if (bst_node_v != NULL) {
				add_obj_set_entry(&data_v->obj_v[i], bst_node_v->p_list_obj);
				bst_node_v->totalWeight += data_v->obj_v[i].weight;
			} else //bst_node_v = NULL
			{
				//Insert a new keyword in the binary tree.
				bst_node_v = (bst_node_t*) malloc(sizeof(bst_node_t));
				memset(bst_node_v, 0, sizeof(bst_node_t));

				/*s*/
				stat_v.memory_v += sizeof(bst_node_t);
				if (stat_v.memory_v > stat_v.memory_max)
					stat_v.memory_max = stat_v.memory_v;
				/*s*/

				//Update the posting list.
				bst_node_v->key = data_v->obj_v[i].fea;
				bst_node_v->p_list_obj = alloc_obj_set();
				bst_node_v->totalWeight = data_v->obj_v[i].weight;

				add_obj_set_entry(&data_v->obj_v[i], bst_node_v->p_list_obj);
				bst_insert(IF_v, bst_node_v);
			}
		}
	}
}
double print_freq(bst_node_t *bst_node_v) {

	if (bst_node_v == NULL)
		return 0;
	printf("%.3lf\t%.3lf\n", bst_node_v->key, bst_node_v->totalWeight);

	print_freq(bst_node_v->left);
	print_freq(bst_node_v->right);

	return 0;
}
double get_highest_freq(bst_node_t *bst_node_v) {

	if (bst_node_v == NULL)
		return 0;

	return fmax(bst_node_v->totalWeight,
			fmax(get_highest_freq(bst_node_v->left),
					get_highest_freq(bst_node_v->right)));
}

double get_highest_obj_n(bst_node_t *bst_node_v) {

	if (bst_node_v == NULL)
		return 0;

	return fmax(bst_node_v->p_list_obj->obj_n,
			fmax(get_highest_obj_n(bst_node_v->left),
					get_highest_obj_n(bst_node_v->right)));
}

fsi_set_t* read_patterns() {
	fsi_set_t *result;
	fsi_t *fsi_v, *fsi_cur;
	ifstream readConfig;

	string STRING;
	vector<string> STRING2, arg;

	result = alloc_fsi_set();

	readConfig.open("pattern.txt");

	//fsi_cur always at the end of the list
	fsi_cur = result->head;
	int i = 0;
	while (!readConfig.eof()) {
		getline(readConfig, STRING);
		//		cout << "STRING: "<< STRING << endl;
		if (STRING.size() == 0)
			break;

		//putting config into vars.
		STRING2 = split(STRING, "  ");

		int size = (int) STRING2.size() - 1;
		fsi_v = alloc_fsi(size);

		for (int j = 0; j < STRING2.size() - 1; j++) {
			//			cout << "STRING2  " << j << ": "<< STRING2[j] << endl;
			fsi_v->feaset[j] = atoi(STRING2[j].c_str());
		}

		//		print_fsi(fsi_v, stdout);

		fsi_cur->next = fsi_v;
		fsi_cur = fsi_cur->next;
		i++;
	}

	return result;
}

/*
 * Additional method for computing support of given fsi_sets
 * It reads "pattern.txt" line by line
 * Each line corrsponds to one fsi_set, the 0 to n-1 numbers are the features
 * Ouput result in console
 * Note that d and min-sup need to be set in config.txt
 //	1: participation
 //	2: fraction-score
 //	3: partitioning
 //	4: construction
 //	5: enumeration
 // 6: gen_scalability_data2 (put in here since it need to read data first)
 */
void colocation_patterns_support(int cost) {
	colocation_config_t *cfg;
	data_t *data_v;

	memset(&stat_v, 0, sizeof(colocation_stat_t));

	//Read the cofig.
	printf("Reading configuration ...\n");
	cfg = read_config_colocation();

	//cost_tag = cfg->cost;
	cost_tag = cost;
	dist_thr = cfg->dist_thr;
	min_sup = cfg->min_sup;
	//min_conf = cfg->min_conf;

	//Read the data.
	printf("Reading data ...\n");
	data_v = read_data_colocation(cfg);

	//Option 1: Build the tree from scratch.
	//Build the IR-tree.
	if (cfg->tree_tag == 0) {
		printf("Building IR-tree ...\n");
		build_IRTree(data_v);

		print_and_check_tree(1, cfg->tree_file);
		//check_IF( );
	} else {
		//Option 2: Read an existing tree.
		printf("Reading IR-Tree ...\n");
		read_tree(cfg->tree_file);
	}

	//Build Inverted file
	build_IF(data_v);
	fea_highest_freq = get_highest_freq(IF_v->root);
	printf("Highest total weight:%lf\n", fea_highest_freq);

	//---
	if (cost_tag == 1)
		printf("participation based:\n");
	else if (cost_tag == 2) {
		printf("Fraction based:\n");
		//Pre-computation of N_o_f
		precomputation(data_v, cfg->dist_thr);

	} else if (cost_tag == 3)
		printf("partitioning based: \n");
	else if (cost_tag == 4)
		printf("construction based: \n");
	else if (cost_tag == 5)
		printf("enumeration based: \n");

	//---
	if (cost == 6)
		gen_scalability_data2(data_v);
	//---

	fsi_set_t *patterns = read_patterns();

	for (int i = 3; i <= 3; i++) {
		cost_tag = i;
		printf("----------- cost_tag:%d -------\n", cost_tag);
		if(i==2)continue;
		fsi_t *fsi_cur = patterns->head->next;
		while (fsi_cur != NULL) {
//			print_fsi(fsi_cur, stdout);

			comp_support(cfg->alg_opt, fsi_cur, cfg->obj_n, NULL, NULL);
			print_fsi(fsi_cur, stdout);
			fsi_cur = fsi_cur->next;
		}
	}
	release_fsi_set(patterns);
	free(cfg);
}

fsi_set_t** find_maximal_fsi_set(fsi_set_t **prev, int numOfFea) {
	fsi_set_t **result;

	// Initialize the structure for storing L_1, L_2, ..., L_{|F|}
	result = (fsi_set_t**) malloc(numOfFea * sizeof(fsi_set_t*));
	memset(result, 0, numOfFea * sizeof(fsi_set_t*));

	/*s*/
	stat_v.memory_v += numOfFea * sizeof(fsi_set_t*);
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	for (int i = 0; i < numOfFea; i++) {
		result[i] = alloc_fsi_set();
	}

	for (int i = 1; i < numOfFea; i++) {
		if (prev[i] != NULL) {
			fsi_t *fsi_v = prev[i]->head->next;

//			if (fsi_v != NULL)
//				fprintf(o_fp, "%d\n", result[i]->fsi_n);

			while (fsi_v != NULL) {
				if (i == numOfFea - 1 || is_maximal(fsi_v, prev[i + 1])) {
//					print_fsi(fsi_v, o_fp);
					fsi_t *fsi_v2 = copy_fsi(fsi_v);
					fsi_v2->sup = fsi_v->sup;
					add_fsi_set_entry(result[fsi_v2->fea_n - 1], fsi_v2);
				}
				fsi_v = fsi_v->next;
			}
//			fprintf(o_fp, "\n");
		}
	}
	return result;
}
