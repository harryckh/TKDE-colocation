/*
 * yao_maximal.cpp
 *
 *  Created on: 19 Apr 2022
 *      Author: kai-ho
 */

#include "yao_maximal_alg.h"

/*
 *  [Yao Maximal]
 */

fsi_set_t** maximal_yao(data_t *data_v, int numOfObj, int numOfFea) {
	fsi_set_t **maximalResult; // storing overall result
	fsi_set_t *resultL1; // storing L1 result
	fsi_set_t *resultL2; // storing L2 result

	unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>> insTable2;

// Initialize the structure for storing L_1, L_2, ..., L_{|F|}
	maximalResult = (fsi_set_t**) malloc(numOfFea * sizeof(fsi_set_t*));
	memset(maximalResult, 0, numOfFea * sizeof(fsi_set_t*));

	/*s*/
	stat_v.memory_v += (numOfFea) * sizeof(fsi_set_t*);
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	for (int i = 0; i < numOfFea; i++) {
		maximalResult[i] = alloc_fsi_set();
	}

#ifndef WIN32
	struct rusage query_sta, query_end;
	float sys_t, usr_t, usr_t_sum = 0;
	GetCurTime(&query_sta);
#endif

	// L_1
	resultL1 = const_L1_apriori();

#ifndef WIN32
	GetCurTime(&query_end);

	GetTime(&query_sta, &query_end, &usr_t, &sys_t);
	printf("L_1:%d\n", resultL1->fsi_n);
	GetCurTime(&query_sta);
#endif
	fsi_set_t *fsi_set_v = NULL;
	fsi_set_t **cmp = NULL;
	int maxSize = 0;

	//Step 1
	maxSize = 0;
	insTable2 = const_ins_table2(data_v);

#ifndef WIN32
	GetCurTime(&query_end);

	GetTime(&query_sta, &query_end, &usr_t, &sys_t);
	printf("P1 \ttime:%0.5lf\t insTable2.size():%d\n", usr_t, insTable2.size());

	GetCurTime(&query_sta);
#endif

//	test_insTable2(insTable2, data_v);
//	exit(-1);
	///===========================================
	//Step 2
	double sup_avg = 0;
	resultL2 = const_L_k_apriori(46, numOfObj, resultL1, 0, sup_avg);

	///===========================================
	//Step 3
	fsi_set_v = NULL;

//			printf("generate graph from size %d!\n", k + 1);
	//gen graph from L_k
	adjacency_list.clear();
	adjacency_list = genGraph(numOfFea, resultL2);

	//Step 2.2 find candidate maximal cliques
	set<vector<int>> maximalClqiues = run_BronKerbosch();

	fsi_set_v = candidate_maximal_patterns(maximalClqiues, maxSize);
	cmp = distribute_by_size(fsi_set_v, maxSize);
	release_fsi_set(fsi_set_v);

	printf("clique maxSize:%d\n", maxSize);

//			for (int i = 0; i < maxSize; i++) {
//				printf("CMP_%d:%d\n", i + 1, cmp[i]->fsi_n);
//			}
#ifndef WIN32
	GetCurTime(&query_end);

	GetTime(&query_sta, &query_end, &usr_t, &sys_t);
	printf("P2 \ttime:%0.5lf\n", usr_t);

	GetCurTime(&query_sta);
#endif
	///===========================================



//	printf("S%d\n", 3);
//	print_fsi_set(fsi_set_v,stdout);

//==================================
//Step 4 compute sup for each candidate
	find_maximal_tree(2, maxSize, numOfObj, cmp, maximalResult, &insTable2);
//============================================



//	printf("---\n");
//	printf("S%d\n", 4);
//	cnt = 0;
//	fsi_t *fsi_v0 = fsi_set_cur2->head->next;
//	while (fsi_v0 != NULL) {
//		print_fsi(fsi_v0, stdout);
//
//		fsi_v0 = fsi_v0->next;
//		cnt++;
//	}

#ifndef WIN32
	GetCurTime(&query_end);

	GetTime(&query_sta, &query_end, &usr_t, &sys_t);
	printf("P3 \ttime:%0.5lf\t\n", usr_t);
#endif

	return maximalResult;
}

//generate a set of feature f_i star neighborhood
unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>> const_ins_table2(
		data_t *data_v) {
	unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>> insTable2;
	vector<obj_set_t*> temp;
	disk_t *disk_v;
	loc_t *loc_v;
	obj_set_t *obj_set_v;
	obj_t *obj_v;

//	printf("const_ins_table2\n");

	//for each object
	for (int i = 0; i < data_v->obj_n; i++) {

		//------
		FEA_TYPE fea1 = data_v->obj_v[i].fea;
		bst_node_t *bst_node_v = bst_search(IF_v, fea1);
		if (bst_node_v != NULL) {
//			double n = bst_node_v->p_list_obj->obj_n;
			double n = bst_node_v->totalWeight;
			if (n < fea_highest_freq * min_sup) {
				continue;
			}
		}
		//----

//		printf("i:%d\n",i);

		obj_v = &data_v->obj_v[i];

		//find all objects in range D(o,d)
		loc_v = get_obj_loc(obj_v);

		disk_v = alloc_disk(IRTree_v.dim);
		set_disk(disk_v, loc_v, dist_thr);

		obj_set_v = range_query(disk_v);

		//  printf("i:%d\tobj_v:%d\n",i,obj_v->id);

//        temp = new vector<obj_set_t*>();

		/*s*/
//		stat_v.memory_v += sizeof(vector<obj_set_t*> );
//		if (stat_v.memory_v > stat_v.memory_max)
//			stat_v.memory_max = stat_v.memory_v;
		/*s*/

		//=======================
		//for each obj in obj_set_v
		obj_node_t *obj_node_v = obj_set_v->head->next;
		while (obj_node_v != NULL) {
			///added: skip if fea must not frequent
			FEA_TYPE fea2 = obj_node_v->obj_v->fea;
			if (bst_search(IF_v, fea2)->totalWeight
					>= fea_highest_freq * min_sup && fea1 > fea2) {
				//Initialize the S_0.
//				obj_set_t *S_0 = alloc_obj_set();
//				add_obj_set_entry(obj_node_v->obj_v, S_0);
//				add_obj_set_entry(obj_v, S_0);

				pair_t pair_v;
				pair_v.o_1 = obj_node_v->obj_v;
				pair_v.o_2 = obj_v;

				/*s*/
				stat_v.memory_v += sizeof(pair_t);
				/*s*/

				string key = create_key_for_table(fea1, fea2);
				auto got = insTable2.find(key);
				if (got == insTable2.end()) //if vector this key is not exist yet
						{
//					vector<obj_set_t*> temp;
//					temp.push_back(S_0);
					unordered_set<pair_t, PairHash, PairEqual> temp;
					temp.insert(pair_v);
					insTable2.insert( { key, temp });

				} else {
//					got->second.push_back(S_0);
					got->second.insert(pair_v);
				}
			}
			obj_node_v = obj_node_v->next;
		}
		//=======================

		release_disk(disk_v);
		release_obj_set(obj_set_v);
		release_loc(loc_v);
	} //end for each object
	/*s*/
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return insTable2;
}

//order-aware key creation
string create_key_for_table(FEA_TYPE fea1, FEA_TYPE fea2) {

	if (fea1 < fea2)
		return to_string(fea1) + "," + to_string(fea2);
	else
		return to_string(fea2) + "," + to_string(fea1);

}

//fsi_set_t* find_size2_pattern(
//		unordered_map<string,unordered_set<pair_t, PairHash, PairEqual>>> insTable2,
//		fsi_set_t *resultL1) {
//
//	fsi_set_t *resultL2 = alloc_fsi_set();
//	fsi_t *fsi_v1, *fsi_v2;
//	FEA_TYPE fea1, fea2;
//
//	fsi_v1 = resultL1->head->next;
//	while (fsi_v1->next != NULL) {
//		fea1 = fsi_v1->feaset[0];
//
//		fsi_v2 = fsi_v2->next;
//		while (fsi_v2 != NULL) {
//			fea2 = fsi_v2->feaset[0];
//
//			//for each possible combination, find instance at table
//			//and obtain its support value
//
//			string key = create_key_for_table(fea1, fea2);
//
//			auto got = insTable2.find(key);
//			if (got != insTable2.end()) {
//				for (auto i : got->second) {
//
//				}
//
//			}
//
//			fsi_v2 = fsi_v2->next;
//		}
//		fsi_v1 = fsi_v1->next;
//	}
//
//	return resultL2;
//}

void find_maximal_tree(int minSize, int maxSize, int numOfObj, fsi_set_t **cmp,
		fsi_set_t **result,
		unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>> *insTable2) {

//	for (int i = 0; i < maxSize; i++) {
//		printf("CMP_%d:%d\n", i + 1, cmp[i]->fsi_n);
//		fsi_t *fsi_cur = cmp[i]->head->next;
//		while (fsi_cur != NULL) {
//			comp_support(46, fsi_cur, numOfObj, NULL,
//			NULL);
//			fsi_cur = fsi_cur->next;
//		}
//		print_fsi_set(cmp[i], stdout);
//		printf("\n");
//	}

	int k = maxSize; //current processing size
	while (k >= minSize) {
//		printf("=======================================\n");
		printf("------ size:%d cmp[k]:%d --------\n", k, cmp[k - 1]->fsi_n);
		fsi_t *fsi_cur = cmp[k - 1]->head->next;
		while (fsi_cur != NULL) {
			if (k < maxSize && checked_subset_fsi(result[k], fsi_cur)) {
				fsi_cur = fsi_cur->next;
				continue;
			}

			double sup = calc_sup_tree(fsi_cur, insTable2);

			//=====
//			double sup = comp_support(46, fsi_cur, numOfObj, NULL, NULL);
//			print_fsi(fsi_cur, stdout);
//			exit(-1);
			if (sup >= min_sup) {
				fsi_t *fsi_v2 = copy_fsi(fsi_cur);
				fsi_v2->sup = sup;
//				print_fsi(fsi_v2, stdout);
				add_fsi_set_entry(result[k - 1], fsi_v2);
			} else if (k > minSize) {
				//depth-first for each subset
				//fsi_next[0-(n-1)] = fsi_cur[1-n]
				fsi_t *fsi_next = alloc_fsi(fsi_cur->fea_n - 1);

				for (int i = 0; i < fsi_next->fea_n; i++) {
					fsi_next->feaset[i] = fsi_cur->feaset[i + 1];
				}
//				if (!(checked_subset_fsi(result, maxSize, fsi_next)
//						|| contain_fsi(cmp[k - 2], fsi_next))) {
//					add_fsi_set_entry(cmp[k - 2], copy_fsi(fsi_next));
//				}
				if (contain_fsi(cmp[k - 2], fsi_next) == NULL
						&& !checked_subset_fsi(result, maxSize, fsi_next)) {
					fsi_t *fsi_new = copy_fsi(fsi_next);
					add_fsi_set_entry(cmp[k - 2], fsi_new);
				}

				//----
				//replace fsi_next[i] by fsi_cur[i] in each iteration
				//cur[i+1] is absent in each iteration in this way
				for (int i = 0; i < fsi_next->fea_n; i++) {
					//skipping i-th label to form k-1 subset
					fsi_next->feaset[i] = fsi_cur->feaset[i];
					//inserting to k-1 level/iteration
					if (contain_fsi(cmp[k - 2], fsi_next) == NULL
							&& !checked_subset_fsi(result, maxSize, fsi_next)) {
						fsi_t *fsi_new = copy_fsi(fsi_next);
						add_fsi_set_entry(cmp[k - 2], fsi_new);
					}

				}
			}
			/*s*/
			if (stat_v.memory_v > stat_v.memory_max)
				stat_v.memory_max = stat_v.memory_v;
			/*s*/
			fsi_cur = fsi_cur->next;
		}
//		printf("result[k-1]->fsi_n:%d\n", result[k - 1]->fsi_n);
		//---
		k--;
	}		//end while
}

B_KEY_TYPE calc_sup_tree(fsi_t *fsi_cur,
		unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>> *insTable2) {
	B_KEY_TYPE sup = INFINITY;

	for (int i = 0; i < fsi_cur->fea_n; i++) {
		obj_set_t *obj_set_v = ((bst_node_t*) bst_search(IF_v,
				fsi_cur->feaset[i]))->p_list_obj;
		int cnt = obj_set_v->obj_n;
		double *frac_receive = new double[cnt];
		double frac_receive_UB = 0.0;
		int j = 0;
		if (cost_tag == 2 && fsi_cur->fea_n > 2) {
			obj_node_t *obj_node_v = obj_set_v->head->next;
			// for each object o with the feature f
			while (obj_node_v != NULL) {
				frac_receive[j] = min_frac_receive(fsi_cur, obj_node_v->obj_v);
				frac_receive_UB += frac_receive[j];
				obj_node_v = obj_node_v->next;
				j++;
			}
			if (frac_receive_UB / fea_highest_freq < min_sup) {
				//				print_fsi(fsi_v,stdout);
				//				printf("frac_receive_UB: %.3lf\n",
				//						frac_receive_UB / fea_highest_freq);
				sup = frac_receive_UB;
				sup = sup / fea_highest_freq;

				fsi_cur->sup = sup;
				free(frac_receive);
//				break; //end this label set, since UB cannot >minsup!
				return sup;
			}
		}
		free(frac_receive);
	}
	//=====
//			printf("generating tree for fsi:\n");
//	unordered_set<int> *RI = new unordered_set<int>();

	//generate CIT of fsi_cur
	cit_t *T = const_cit(fsi_cur, insTable2);


//	set_RI(T, RI, fsi_cur->fea_n);
//			printf("calculating sup from tree\n");

	//if T-> depth == k, valid and add sup of each child obj_
//			double sup = add_sup_tree(T, fsi_cur);
//			sup = sup / fea_highest_freq;

	// for each feature f in C
	for (int i = 0; i < fsi_cur->fea_n; i++) {
		B_KEY_TYPE sup_C_f = 0;
		// the corresponding inverted list in IF
		obj_set_t *obj_set_v = ((bst_node_t*) bst_search(IF_v,
				fsi_cur->feaset[i]))->p_list_obj;
		obj_node_t *obj_node_v = obj_set_v->head->next;
		// for each object o with the feature f
		while (obj_node_v != NULL) {
//			if (RI->find(obj_node_v->obj_v->id) != RI->end()) {
			if (is_instance(T, fsi_cur->fea_n, i,obj_node_v->obj_v)) {
				sup_C_f += min_frac_receive(fsi_cur, obj_node_v->obj_v);
			}
			obj_node_v = obj_node_v->next;
		} //end for each obj

//				printf("i:%d sup:%lf\n",i,sup_C_f / fea_highest_freq);
		if (sup_C_f <= sup) {
			sup = sup_C_f;
		}
	} //end for each label
	sup = sup / fea_highest_freq;

	fsi_cur->sup = sup;

//			print_fsi(fsi_cur, stdout);

	cit_release(T);

	return sup;
}

//construct cit-tree for @fsi_v from insTable2
cit_t* const_cit(fsi_t *fsi_v,
		unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>> *insTable2) {
	cit_t *T = cit_ini();
	cit_node_t *cit_node_v1, *cit_node_v2;

	vector<cit_node_t*> size2node;

//	printf("const tree step.1\n");
	//step 1
	FEA_TYPE fea1, fea2;
	fea1 = fsi_v->feaset[0];
	fea2 = fsi_v->feaset[1];

	string key = create_key_for_table(fea1, fea2);
//	printf("key:%s\n", key.c_str());

	auto got = insTable2->find(key);
	if (got == insTable2->end()) //empty tree for this fea pair
		return T;

//	printf("got != end() %d\n", got->second.size());
	//got != end()
//	for (obj_set_t *obj_set_v : got->second) {
	for (pair_t pair_v : got->second) {
//		print_obj_set(obj_set_v,stdout);
//		obj_t *obj_v1 = obj_set_v->head->next->obj_v;
//		obj_t *obj_v2 = obj_set_v->head->next->next->obj_v;
		obj_t *obj_v1 = pair_v.o_1;
		obj_t *obj_v2 = pair_v.o_2;
//		printf("v1:%d %d v2:%d %d\n",obj_v1->id, obj_v1->fea, obj_v2->id, obj_v2->fea);
		cit_node_v1 = is_child(T->root, obj_v1);
		if (cit_node_v1 != NULL) {
			//insert v2 as child of v1
			cit_insert(T, cit_node_v1, obj_v2);
		} else {
			//create two nodes for first two fea of each instance pair in table
			cit_node_v1 = cit_insert(T, T->root, obj_v1);
			cit_node_v2 = cit_insert(T, cit_node_v1, obj_v2);
		}

//		if (fsi_v->fea_n == 2) { //size-2 maximal
//			set_ancestor_flag(cit_node_v2, RI);
//		}
//		printf("inserted\n");
	}

	//===============================================

	//step 2
	int i = 2; //from 3rd fea, i.e., [2]
	while (i < fsi_v->fea_n) {
//		printf("i:%d\n", i);
		//construct node
		fea1 = fsi_v->feaset[i - 1];
		fea2 = fsi_v->feaset[i];
		key = create_key_for_table(fea1, fea2);
		auto got = insTable2->find(key);
		if (got == insTable2->end())
			break;

		//obj_v2 are candidates for create new node at i-th level
		//as child of obj_v1
//		printf(" got->second:%lu\n", got->second.size());
//			for (obj_set_t *obj_set_v : got->second) {
//				obj_t *obj_v1 = obj_set_v->head->next->obj_v;
//				obj_t *obj_v2 = obj_set_v->head->next->next->obj_v;
		for (pair_t pair_v : got->second) {
			obj_t *obj_v1 = pair_v.o_1;
			obj_t *obj_v2 = pair_v.o_2;
//				printf("obj_v1:%d v2:%d\n", obj_v1->id, obj_v2->id);

			//find all possible parents
			//i-1 because searching parent level
			vector<cit_node_t*> list = cit_search(T, i - 1, obj_v1);

//				if(list.size()!=1)
//				printf("list.size:%lu\n",list.size());

			//for each node in the list check whether it is qualified for a new child node
			//by checking all its ancestor

			for (cit_node_t *cit_node_v : list) {
				if (check_clique(T, cit_node_v->parent, obj_v2, insTable2)) {
					//satisfy and add
					cit_node_v2 = cit_insert(T, cit_node_v, obj_v2);
//					if (i == fsi_v->fea_n - 1) {
////							printf("set flag\n");
//						set_ancestor_flag(cit_node_v2, RI);
////							print_ancestor(cit_node_v2);
//
//					}
				}
			}

		}

//		if (!insertedChild)
//			break;
		i++;
	}				//end while
//	printf("check tree2: %d\n", check_tree(T, fsi_v));

	return T;

}

//return the child node if obj_v is child of cit_node_v
//otherwise return null
cit_node_t* is_child(cit_node_t *cit_node_v, obj_t *obj_v) {
	cit_node_t *cit_node_child;

	if (cit_node_v->child_num == 0)
		return NULL;

	cit_node_child = cit_node_v->child->next;
	while (cit_node_child != NULL) {
		if (cit_node_child->obj_v == obj_v)
			return cit_node_child;
		cit_node_child = cit_node_child->next;
	}

	return NULL;
}

void set_ancestor_flag(cit_node_t *cit_node_v, unordered_set<int> *RI) {

	if (cit_node_v->parent != NULL) {

		cit_node_v->flag = true;
		RI->insert(cit_node_v->obj_v->id);
		set_ancestor_flag(cit_node_v->parent, RI);
	}
}

void cit_search_sub(cit_node_t *cit_node_v, int targetDepth,
		vector<cit_node_t*> &list, obj_t *obj_v, int curDepth) {

//	printf("targetDepth:%d\t curDepth:%d\n",targetDepth,curDepth);

	if (curDepth < targetDepth && cit_node_v->child_num > 0) {
		cit_node_t *cit_node_child = cit_node_v->child->next;
		while (cit_node_child != NULL) {
			cit_search_sub(cit_node_child, targetDepth, list, obj_v,
					curDepth + 1);
			cit_node_child = cit_node_child->next;
		}
	}
	if (targetDepth == curDepth) {
//		printf("targetDepth:%d\t ==curDepth:%d\n",targetDepth,curDepth);
//		printf("%d %d\n",cit_node_v->obj_v->id,obj_v->id);
		if (cit_node_v->obj_v->id == obj_v->id) {
			cit_node_t *cit_node_tmp = cit_node_v;
			list.push_back(cit_node_tmp);
		}
	}
	return;
}

//search all nodes with obj_v at level i
vector<cit_node_t*> cit_search(cit_t *T, int targetDepth, obj_t *obj_v) {

	vector<cit_node_t*> list;

	cit_search_sub(T->root, targetDepth, list, obj_v, -1); //root = -1

	return list;
}

bool check_clique(cit_t *T, cit_node_t *cit_node_v, obj_t *obj_v,
		unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>> *insTable2) {

	if (cit_node_v == T->root)
		return true;

	FEA_TYPE fea1 = cit_node_v->obj_v->fea;
	FEA_TYPE fea2 = obj_v->fea;
	string key = create_key_for_table(fea1, fea2);

	auto got = insTable2->find(key);
	if (got != insTable2->end()) {
//		for (auto obj_set_v : got->second) {
//			if (cit_node_v->obj_v == obj_set_v->head->next->obj_v
//					&& obj_v == obj_set_v->head->next->next->obj_v)
//				return check_clique(T, cit_node_v->parent, obj_v, insTable2);
//		}
		pair_t pair_v;
		pair_v.o_1 = cit_node_v->obj_v;
		pair_v.o_2 = obj_v;
		auto got2 = got->second.find(pair_v);
		if (got2 != got->second.end()) {
			return check_clique(T, cit_node_v->parent, obj_v, insTable2);
		}
	}

	return false;
}

void add_sup_tree_sub(cit_node_t *cit_node_v, int depth, fsi_t *fsi_v,
		double *sups) {

	//handle this level node
	if (depth > -1 && cit_node_v->flag)
		sups[depth] += min_frac_receive(fsi_v, cit_node_v->obj_v);

	//recursive
	if (cit_node_v->child_num > 0) {
		cit_node_t *cit_node_child = cit_node_v->child->next;
		while (cit_node_child != NULL) {
			add_sup_tree_sub(cit_node_child, depth + 1, fsi_v, sups);
			cit_node_child = cit_node_child->next;
		}
	}
}

//
B_KEY_TYPE add_sup_tree(cit_t *T, fsi_t *fsi_v) {

	B_KEY_TYPE sup = INFINITY;

	double *sups; //1d array each corresponds to sup of group by fea

	sups = new double[fsi_v->fea_n];

	add_sup_tree_sub(T->root, -1, fsi_v, sups);

	//return the min one
	for (int i = 0; i < fsi_v->fea_n; i++) {
		if (sups[i] < sup)
			sup = sups[i];
	}

	return sup;
}

//set to true to all its ancestor node if it has depth == fsi_n
void set_RI_sub(cit_node_t *cit_node_v, unordered_set<int> *RI, int num,
		int depth) {

	if (num == depth) {
		//if number of fea equal to cur level
		set_ancestor_flag(cit_node_v, RI);
		return;
	}

	//recursive
	if (cit_node_v->child_num > 0) {
		cit_node_t *cit_node_child = cit_node_v->child->next;
		while (cit_node_child != NULL) {
			set_RI_sub(cit_node_child, RI, num, depth + 1);
			cit_node_child = cit_node_child->next;
		}
	}
}

void set_RI(cit_t *T, unordered_set<int> *RI, int num) {
	set_RI_sub(T->root, RI, num, 0);
}

bool is_instance_sub(cit_node_t *cit_node_v, int num, int curDepth) {

	if (curDepth+1==num) {
		//if number of fea equal to cur level
		return true;
	}

	//recursive
	if (cit_node_v->child_num > 0) {
		cit_node_t *cit_node_child = cit_node_v->child->next;
		while (cit_node_child != NULL) {
			if (is_instance_sub(cit_node_child, num, curDepth + 1))
				return true;
			cit_node_child = cit_node_child->next;
		}
	}
	return false;
}

bool is_instance(cit_t *T, int num, int targetDepth, obj_t *obj_v) {

	vector<cit_node_t*> list = cit_search(T, targetDepth, obj_v);

//	printf("list.size:%d\n",list.size());

	if (list.size() == 0)
		return false;
	for (cit_node_t *cit_node_v : list) {
		//true if any one node is true
		if (is_instance_sub(cit_node_v, num, targetDepth))
			return true;
	}
	return false;
}

//check if we can get sth from the insTable2
void test_insTable2(
		unordered_map<string, unordered_set<pair_t, PairHash, PairEqual>> insTable2,
		data_t *data_v) {
	fsi_t *fsi_cur = alloc_fsi(2);
	fsi_cur->feaset[0] = 200;
	fsi_cur->feaset[1] = 202;

	FEA_TYPE fea1 = fsi_cur->feaset[0];
	FEA_TYPE fea2 = fsi_cur->feaset[1];
	string key = create_key_for_table(fea1, fea2);

	double *sups = new double[2];

	auto got = insTable2.find(key);
	if (got != insTable2.end()) {
//		for (auto pair_v : got->second) {
//
//			sups[0] += min_frac_receive(fsi_cur, pair_v.o_1);
//			sups[1] += min_frac_receive(fsi_cur,pair_v.o_2);
//			printf("o1:%d %lf o2:%d %lf\n", pair_v.o_1->id,min_frac_receive(fsi_cur, pair_v.o_1),
//					pair_v.o_2->id,min_frac_receive(fsi_cur, pair_v.o_2));
//		}
		pair_t pair_v;
		pair_v.o_1 = &(data_v->obj_v[31116 - 1]);
		pair_v.o_2 = &(data_v->obj_v[31118 - 1]);
		printf("o1:%d %u o2:%d %u\n", pair_v.o_1->id, pair_v.o_1->fea,
				pair_v.o_2->id, pair_v.o_2->fea);
		auto got2 = got->second.find(pair_v);
		if (got2 != got->second.end()) {
			printf("o1:%d %lf o2:%d %lf\n", pair_v.o_1->id,
					min_frac_receive(fsi_cur, pair_v.o_1), pair_v.o_2->id,
					min_frac_receive(fsi_cur, pair_v.o_2));
		} else {
			printf("not found!\n");
		}
	}

	printf("sups[0]:%lf \t [1]:%lf\n", sups[0] / fea_highest_freq,
			sups[1] / fea_highest_freq);

//
//	printf("generating tree for fsi:\n");
//	print_fsi(fsi_cur, stdout);
//
//	//generate CIT of fsi_cur
//	cit_t *T = const_cit(fsi_cur, insTable2);
//
//	printf("calculating sup from tree\n");
//
//	//if T-> depth == k, valid and add sup of each child obj_
//	double sup = add_sup_tree(T, fsi_cur);
//	sup = sup / fea_highest_freq;
//	fsi_cur->sup = sup;
//
//	print_fsi(fsi_cur, stdout);
//
//	cit_release(T);

}
