/*
 *	Author: Harry
 *	Email: khchanak@cse.ust.hk
 */

#include "yoo_alg.h"

/*
 * The implementation of joinless colocation mining algorithm.
 */
fsi_set_t** joinless_mining(data_t *data_v, int numOfObj, int numOfFea) {
	unordered_map<FEA_TYPE, vector<obj_set_t*>> SN;

	fsi_set_t **result; //storing overall result
	fsi_set_t *fsi_set_cur; //prev = L_{k-1}, cur = L_{k}
	fsi_t *fsi_v1, *fsi_v2, *fsi_v3;
	FEA_TYPE fea_v1;
	B_KEY_TYPE sup;
	int i;

	//Initialize the structure for storing L_1, L_2, ..., L_{|F|}
	result = (fsi_set_t**) malloc(numOfFea * sizeof(fsi_set_t*));
	memset(result, 0, numOfFea * sizeof(fsi_set_t*));

	/*s*/
	stat_v.memory_v += numOfFea * sizeof(fsi_set_t*);
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	SN = gen_star_neighborhoods(data_v);

	/*s*/
	for (auto it = SN.begin(); it != SN.end(); ++it)
		stat_v.memory_v += it->second.size() * sizeof(obj_set_t*);
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

#ifndef WIN32
	struct rusage query_sta, query_end;
	float sys_t, usr_t, usr_t_sum = 0;
	GetCurTime(&query_sta);
#endif

	//L_1
	result[0] = const_L1_apriori();

#ifndef WIN32
	GetCurTime(&query_end);

	GetTime(&query_sta, &query_end, &usr_t, &sys_t);
	printf("L_1 \ttime:%0.5lf result[0]:%d\n", usr_t, result[0]->fsi_n);

	GetCurTime(&query_sta);
#endif
	//L_(i+1)
	for (i = 1; i < numOfFea; i++) {
		fsi_v1 = result[i - 1]->head->next;
		if (fsi_v1 == NULL)
			break;

		fsi_set_cur = alloc_fsi_set();
		while (fsi_v1->next != NULL) {
			fsi_v2 = fsi_v1->next;
			while (fsi_v2 != NULL) {
				//join step.
				if ((fea_v1 = join_check(fsi_v1, fsi_v2)) != -1) {
					fsi_v3 = add_fsi(fsi_v1, fea_v1);

					//prune step.
					if (i > 1 && !all_subesets_exist(result[i - 1], fsi_v3)) {
						/*s*/
						stat_v.memory_v -= sizeof(fsi_t)
								+ fsi_v3->fea_n * sizeof(FEA_TYPE);
						/*s*/
						release_fsi(fsi_v3);
						fsi_v2 = fsi_v2->next;
						continue;
					}

					//----------------------------
					check_star_instance(fsi_v3, SN);

					if (fsi_v3->fea_n > 2) {
						if (cost_tag == 1) //participation based
							sup = comp_PI(fsi_v3, NULL);
						else
							//fraction based
							sup = comp_sup(fsi_v3, NULL);

						if (sup < min_sup) {
							//fsi_v3 is not a frequent pattern
							for (auto i = fsi_v3->obj_set_list_v.begin();
									i != fsi_v3->obj_set_list_v.end(); ++i)
								release_obj_set(*i);

							/*s*/
							stat_v.memory_v -= fsi_v3->obj_set_list_v.size()
									* sizeof(obj_set_t*);
							/*s*/

//                            delete (fsi_v3->obj_set_list_v);
							vector<struct obj_set*>().swap(
									fsi_v3->obj_set_list_v);

							release_fsi(fsi_v3);

							fsi_v2 = fsi_v2->next;
							continue;
						}
						filter_clique_instance(fsi_v3, result[i - 1],
								&fsi_v3->obj_set_list_v);
					}
					//count step.
					if (cost_tag == 1) //participation based
						sup = comp_PI(fsi_v3, NULL);
					else
						//fraction based
						sup = comp_sup(fsi_v3, NULL);

					//----------------------------
					if (sup >= min_sup)
						add_fsi_set_entry(fsi_set_cur, fsi_v3);
					else {
						//fsi_v3 is not a frequent pattern
						for (auto i = fsi_v3->obj_set_list_v.begin();
								i != fsi_v3->obj_set_list_v.end(); ++i)
							release_obj_set(*i);

						/*s*/
						stat_v.memory_v -= fsi_v3->obj_set_list_v.size()
								* sizeof(obj_set_t*);
						stat_v.memory_v -= sizeof(fsi_t)
								+ fsi_v3->fea_n * sizeof(FEA_TYPE);
						/*s*/

//                        delete (fsi_v3->obj_set_list_v);
						vector<struct obj_set*>().swap(fsi_v3->obj_set_list_v);
						release_fsi(fsi_v3);
					}
				}
				fsi_v2 = fsi_v2->next;
			}
			fsi_v1 = fsi_v1->next;
		}
#ifndef WIN32
		GetCurTime(&query_end);

		GetTime(&query_sta, &query_end, &usr_t, &sys_t);
		printf("L_%d\ttime:%0.5lf\t fsi_n:%d\n", i + 1, usr_t,
				fsi_set_cur->fsi_n);
		GetCurTime(&query_sta);
#endif
		/*t*/
		//print out each level for debug
		if (debug_mode) {
			printf("L_%d\n", i + 1);
			int cnt = 0;
			fsi_t *fsi_v = fsi_set_cur->head->next;
			while (fsi_v != NULL) {
				print_fsi(fsi_v, stdout);

				fsi_v = fsi_v->next;
				cnt++;
			}
		}
		//printf("L_%d: %d\n",i+1, fsi_set_cur->fsi_n);
		result[i] = fsi_set_cur;
	}

	return result;
}

//generate a set of feature f_i star neighborhood
unordered_map<FEA_TYPE, vector<obj_set_t*>> gen_star_neighborhoods(
		data_t *data_v) {
	unordered_map<FEA_TYPE, vector<obj_set_t*>> SN;
	vector<obj_set_t*> temp;
	k_node_t *k_head;
	disk_t *disk_v;
	loc_t *loc_v;
	obj_set_t *obj_set_v;
	obj_t *obj_v;

	printf("gen_star_neighborhoods\n");

	k_head = collect_keywords_bst(IF_v);

	//for each object
	for (int i = 0; i < data_v->obj_n; i++) {

		//------
		FEA_TYPE fea = data_v->obj_v[i].fea;
		bst_node_t *bst_node_v = bst_search(IF_v, fea);
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
		stat_v.memory_v += sizeof(vector<obj_set_t*> );
		if (stat_v.memory_v > stat_v.memory_max)
			stat_v.memory_max = stat_v.memory_v;
		/*s*/

		temp = gen_one_star_neighborhood(obj_v, obj_set_v);

		auto got = SN.find(obj_v->fea);
		if (got == SN.end()) //if the obj_set_list for this feature is not exist yet
				{
			SN.insert( { fea, temp });

		} else {
			//copy temp to the end of corresponding list in SN
			got->second.insert(got->second.end(), temp.begin(), temp.end());
//            delete (temp);
			vector<struct obj_set*>().swap(temp);
		}
		release_disk(disk_v);
		release_obj_set(obj_set_v);
		release_loc(loc_v);
	}

	release_k_list(k_head);

	return SN;
}

/*
 *  generate one star neighborhood and insert into temp
 *  obj in @obj_set_v with fea > obj_v->fea are included in the star neighborhood
 *  @obj_v is the first (i.e, head->next) element in S_0
 *
 */
vector<obj_set_t*> gen_one_star_neighborhood(obj_t *obj_v,
		obj_set_t *obj_set_v) {
	vector<obj_set_t*> temp;
	obj_set_t *S_0;
	obj_node_t *obj_node_v;

	//Initialize the S_0.
	S_0 = alloc_obj_set();

	obj_node_v = obj_set_v->head->next;
	while (obj_node_v != NULL) {
		///added: skip if fea must not freqent
		FEA_TYPE fea = obj_node_v->obj_v->fea;
		if (bst_search(IF_v, fea)->totalWeight >= fea_highest_freq * min_sup
				&& fea > obj_v->fea)
			add_obj_set_entry(obj_node_v->obj_v, S_0);
		obj_node_v = obj_node_v->next;
	}

	add_obj_set_entry(obj_v, S_0);
	temp.push_back(S_0);

	return temp;
}

/*
 *  1. Find all row instances of fsi_v
 *  (based on the materialized star instances @SN and by filtering technique)
 *  2. Result stored in fsi_v->obj_set_list_v
 */
void check_star_instance(fsi_t *fsi_v,
		unordered_map<FEA_TYPE, vector<obj_set_t*>> SN) {

	vector<obj_set_t*> obj_set_list_v, obj_set_list_temp;
	obj_set_t *obj_set_v;
	bst_t *inverted_list;
	psi_t *psi_v;

	//    if(fsi_v->fea_n==3)
	//    {
	//        print_fsi(fsi_v, stdout);
	//    }

	//finding star neighborhoods with the first feature in fsi_v
	auto got = SN.find(fsi_v->feaset[0]);
	if (got == SN.end()) {
		printf("feature %d not found!! exiting.\n", fsi_v->feaset[0]);
		exit(-1);
	}

	//all star neigbhorhoods with first object having fsi_v->feaset[0]
	obj_set_list_v = got->second;

	psi_v = fsi_to_psi(fsi_v, fsi_v->feaset[0]); //excluding the first one
//	printf("obj_set_list_v.size:%lu \n",obj_set_list_v.size());
			//-------------------------------------------------------------------------
			//for each star neighborhood
	for (auto i = obj_set_list_v.begin(); i != obj_set_list_v.end(); ++i) {
		obj_set_v = *i;

		inverted_list = const_IF(obj_set_v, psi_v);

		//generate star instances
		obj_set_list_temp = filter_star_instance(obj_set_v->head->next->obj_v,
				inverted_list);
		release_IF(inverted_list);

		/*t*/
		//        for (auto j = obj_set_list_temp->begin(); j != obj_set_list_v->end(); ++j)
		//            print_obj_set(*j, stdout);
		/*t*/

		//append obj_set_list_temp to fsi_v->obj_set_list_v for k+1 to filter clique instance
//        if (fsi_v->obj_set_list_v == NULL) {
		if (fsi_v->obj_set_list_v.size() == 0) {
//			printf("hihi\n");
			fsi_v->obj_set_list_v = obj_set_list_temp;
		} else {
//			printf("qwqw\n");
			//copy the obj_set_list_temp to the end of fsi_v->obj_set_list_v
			fsi_v->obj_set_list_v.insert(fsi_v->obj_set_list_v.end(),
					obj_set_list_temp.begin(), obj_set_list_temp.end());
//            delete (obj_set_list_temp);
			vector<struct obj_set*>().swap(obj_set_list_temp);
		}
	}

	release_psi(psi_v);
	return;
}

/*
 * generate star instances
 * @obj_v: first object in the instance
 * @inverted_list: generated from the corresponding star neighborhood
 * return @obj_set_list_v: start instances each involves @obj_v
 */
vector<obj_set_t*> filter_star_instance(obj_t *obj_v, bst_t *inverted_list) {
	vector<obj_set_t*> obj_set_list_v = { };
	obj_set_t *S_0;

	S_0 = alloc_obj_set();
	add_obj_set_entry(obj_v, S_0);

//    obj_set_list_v = new vector<obj_set_t*>();

	filter_star_instance_sub(inverted_list, S_0, &obj_set_list_v);

	/*s*/
	stat_v.memory_v += obj_set_list_v.size() * sizeof(obj_set_t*)
			+ sizeof(vector<obj_set_t*> );
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	release_obj_set(S_0);

	return obj_set_list_v;
}

/*
 *	The sub-procedure
 *
 */
void filter_star_instance_sub(bst_t *IF_v, obj_set_t *S_0,
		vector<obj_set_t*> *obj_set_list_v) {

	obj_t *obj_v;
	bst_node_t *bst_node_v;
	obj_node_t *obj_node_v;
	bst_node_list_t *bst_node_list_v;

	if (IF_v->node_n == 0) {
		obj_set_t *S;
		S = copy_obj_set(S_0);
		obj_set_list_v->push_back(S);
		return;
	}

	bst_node_v = IF_v->root;
	obj_node_v = bst_node_v->p_list_obj->head->next;
	while (obj_node_v != NULL) {
		//Pick an object.
		obj_v = obj_node_v->obj_v;

		//Update the IF_v.
		bst_node_list_v = update_IF_obj(IF_v, obj_v);

		//Update the S_0.
		//obj_v is added at the first place of S_0.
		add_obj_set_entry(obj_v, S_0);

		//Sub-procedure.
		filter_star_instance_sub(IF_v, S_0, obj_set_list_v);

		//Restore the S_0.
		remove_obj_set_entry(S_0);

		//Restore the IF_v.
		restore_IF_bst_node_list(IF_v, bst_node_list_v);

		release_bst_node_list(bst_node_list_v);

		//Try the next object candidate.
		obj_node_v = obj_node_v->next;
	}
	return;
}

//for each obj_set_v in fsi_cur->obj_set_list_v, check clique instance for objects except the one with fsi_cur->feaset[0]
//remove it from the list if not clique
//checking can be done by searching C_k-1
void filter_clique_instance(fsi_t *fsi_cur, fsi_set_t *fsi_set_v,
		vector<obj_set_t*> *obj_set_list_v) {
	fsi_t *fsi_v;
	obj_node_t *obj_node_temp;
	//-------------------------
	//first we find the corresponding fsi
	// e.g., fsi_cur = {A, B, C}; fsi_v = {B, C}
	fsi_v = fsi_set_v->head->next;
	while (fsi_v != NULL) {
		//if(fsi_cur->fea[1-] == fsi_v->fea[]
		if (memcmp(&fsi_cur->feaset[1], fsi_v->feaset,
				(fsi_v->fea_n) * sizeof(FEA_TYPE)) == 0)
			break;
		fsi_v = fsi_v->next;
	}
	//fsi_v == correspodning fsi || NULL
	if (fsi_v == NULL) {
		printf("fsi_v not found!! exiting\n");
		print_fsi(fsi_v, stdout);
		exit(-1);
	}

	//-------------------
	//second we check whether obj_set_v (i) exist in fsi_v->obj_set_list-v
	for (auto i = obj_set_list_v->begin(); i != obj_set_list_v->end();
			/*++i*/) {
		bool tag = false;
		for (auto j = fsi_v->obj_set_list_v.begin();
				j != fsi_v->obj_set_list_v.end(); ++j) {
			if (check_obj_set_equal(*i, *j, fsi_cur->feaset[0])) {
				tag = true;
				break;
			}
		}

		if (!tag) //not exist
		{
			//remove current obj_set_v from obj_set_list_v
			release_obj_set(*i);
			obj_set_list_v->erase(i);
		} else {
			//the current object set is clique
			i++;
		}
	}
}

/*
 * check whether two object sets contain the same set of objects
 * obj with fea in v1 is not checked
 */
bool check_obj_set_equal(obj_set_t *v1, obj_set_t *v2, FEA_TYPE fea) {
	obj_node_t *obj_node_v1, *obj_node_v2;

	if (v1->obj_n - 1 != v2->obj_n) {
		printf("error!\n");
		print_obj_set(v1, stdout);
		print_obj_set(v2, stdout);
		exit(-1);
	}

	obj_node_v1 = v1->head->next;
	while (obj_node_v1 != NULL) {
		if (obj_node_v1->obj_v->fea == fea) {
			obj_node_v1 = obj_node_v1->next;
			continue;
		}
		bool tag = false;
		obj_node_v2 = v2->head->next;
		while (obj_node_v2 != NULL) {
			if (obj_node_v1->obj_v == obj_node_v2->obj_v) {
				tag = true;
				break;
			}
			obj_node_v2 = obj_node_v2->next;
		}
		if (!tag)
			return false;

		obj_node_v1 = obj_node_v1->next;
	}
	return true;
}

/*
 *  Output: @sup = sup(C)
 */
B_KEY_TYPE comp_PI(fsi_t *fsi_v, obj_set_t *O) {

	double participation_index, participation_ratio;
	obj_node_t *obj_node_v;
	obj_set_t *obj_set_v;

	participation_index = INFINITY;

	//for each feature f in C
	for (int i = 0; i < fsi_v->fea_n; i++) {
		participation_ratio = 0;

		if (O != NULL) //L1
			obj_set_v = O;
		else
			//L_k
			//the corresponding inverted list in IF
			obj_set_v =
					((bst_node_t*) bst_search(IF_v, fsi_v->feaset[i]))->p_list_obj;

		obj_node_v = obj_set_v->head->next;

		//for each object o with the feature f
		while (obj_node_v != NULL) {
			if (obj_exist(fsi_v->obj_set_list_v, obj_node_v->obj_v)) {
				participation_ratio += 1;
			}

			obj_node_v = obj_node_v->next;
		}

		participation_ratio = participation_ratio / (double) obj_set_v->obj_n;

		if (participation_ratio <= participation_index)
			participation_index = participation_ratio;
	}

	fsi_v->sup = participation_index;

	return participation_index;
}

/*
 *  Output: @sup = sup(C)
 */
B_KEY_TYPE comp_sup(fsi_t *fsi_v, obj_set_t *O) {

	B_KEY_TYPE sup, sup_C_f;
	obj_node_t *obj_node_v;
	obj_set_t *obj_set_v;
	loc_t *loc_v;
	disk_t *disk_v;

	sup = INFINITY;
	sup_C_f = 0;

	//for each feature f in C
	for (int i = 0; i < fsi_v->fea_n; i++) {
		sup_C_f = 0;

		if (O != NULL) //L1
			obj_node_v = O->head->next;
		else
			//L_k
			//the corresponding inverted list in IF
			obj_node_v =
					((bst_node_t*) bst_search(IF_v, fsi_v->feaset[i]))->p_list_obj->head->next;

		//for each object o with the feature f
		while (obj_node_v != NULL) {
			if (obj_exist(fsi_v->obj_set_list_v, obj_node_v->obj_v)) {
				sup_C_f += min_frac_receive(fsi_v, obj_node_v->obj_v);
				if (sup_C_f > sup) {
					// stat_v.S4_sum++;
					break;
				}
			}

			obj_node_v = obj_node_v->next;
		}

		if (sup_C_f <= sup)
			sup = sup_C_f;
	}

	sup = sup / fea_highest_freq;
	fsi_v->sup = sup;

	return sup;
}

// check whether the object @obj_v exist in any of the object set in @obj_set_list_v
bool obj_exist(vector<obj_set_t*> obj_set_list_v, obj_t *obj_v) {
	obj_node_t *obj_node_v;
	for (auto i = obj_set_list_v.begin(); i < obj_set_list_v.end(); i++) {
		obj_node_v = (*i)->head->next;
		while (obj_node_v != NULL) {
			if (obj_node_v->obj_v == obj_v)
				return true;
			obj_node_v = obj_node_v->next;
		}
	}
	return false;
}

