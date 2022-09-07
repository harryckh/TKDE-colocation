/*
 *	Author: Harry
 *	Email: khchanak@cse.ust.hk
 */

#include "frac.h"

/*
 *  [Apriori]
 */

fsi_set_t** apriori(int alg_opt, int numOfObj, int numOfFea) {
	int i, j;
	fsi_set_t **result; // storing overall result
	fsi_set_t *fsi_set_cur; // prev = L_{k-1}, cur = L_{k}
	fsi_t *fsi_v1, *fsi_v2, *fsi_v3;

	FEA_TYPE fea_v1;
	B_KEY_TYPE sup;
	double avg_sup;

	i = 0;
	// Initialize the structure for storing L_1, L_2, ..., L_{|F|}
	result = (fsi_set_t**) malloc(numOfFea * sizeof(fsi_set_t*));
	memset(result, 0, numOfFea * sizeof(fsi_set_t*));

	/*s*/
	stat_v.memory_v += numOfFea * sizeof(fsi_set_t*);
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

#ifndef WIN32
	struct rusage query_sta, query_end;
	float sys_t, usr_t, usr_t_sum = 0;
	GetCurTime(&query_sta);
#endif

	// L_1
	result[0] = const_L1_apriori();

#ifndef WIN32
	GetCurTime(&query_end);

	GetTime(&query_sta, &query_end, &usr_t, &sys_t);
	printf("L_1 \ttime:%0.5lf result[0]:%d\n", usr_t, result[0]->fsi_n);
	GetCurTime(&query_sta);
#endif

	// L_(i)
	for (i = 1; i < numOfFea; i++) {
		avg_sup = 0;
		fsi_v1 = result[i - 1]->head->next;
		if (fsi_v1 == NULL)
			break;
		fsi_set_cur = alloc_fsi_set();
		while (fsi_v1->next != NULL) {
			fsi_v2 = fsi_v1->next;
			while (fsi_v2 != NULL) {
				// join step.
				if ((fea_v1 = join_check(fsi_v1, fsi_v2)) != -1) {
					fsi_v3 = add_fsi(fsi_v1, fea_v1);

					// prune step.
					if (i > 1 && !all_subesets_exist(result[i - 1], fsi_v3)) {
						/*s*/
						stat_v.memory_v -= sizeof(fsi_t)
								+ fsi_v3->fea_n * sizeof(FEA_TYPE);
						/*s*/

						release_fsi(fsi_v3);
						fsi_v2 = fsi_v2->next;
						continue;
					}
					// count step.
					sup = comp_support(alg_opt, fsi_v3, numOfObj, NULL, NULL);
					if (sup >= min_sup) {
						add_fsi_set_entry(fsi_set_cur, fsi_v3);
						avg_sup += sup;
					} else {
						/*s*/
						stat_v.memory_v -= sizeof(fsi_t)
								+ fsi_v3->fea_n * sizeof(FEA_TYPE);
						/*s*/
						release_fsi(fsi_v3);
					}
				}
				fsi_v2 = fsi_v2->next;
			}
			fsi_v1 = fsi_v1->next;
		}
		//----
#ifndef WIN32
		GetCurTime(&query_end);

		GetTime(&query_sta, &query_end, &usr_t, &sys_t);
		double tmp = 0;
		if (fsi_set_cur->fsi_n > 0)
			tmp = avg_sup / (double) fsi_set_cur->fsi_n;
//		printf("L_%d \ttime:%0.5lf\t fsi_n:%d avg_sup:%0.5lf\n", i + 1, usr_t,
//				fsi_set_cur->fsi_n, tmp);
		printf("L_%d \ttime:%0.5lf\t size:%d\n", i + 1, usr_t,
						fsi_set_cur->fsi_n);
		GetCurTime(&query_sta);
#endif
		/*t*/
		// print out each level for debug
		if (false) {
			printf("L_%d\n", i + 1);
			int cnt = 0;
			fsi_t *fsi_v = fsi_set_cur->head->next;
			while (fsi_v != NULL) {
				print_fsi(fsi_v, stdout);

				fsi_v = fsi_v->next;
				cnt++;
			}
		}
		// printf("L_%d: %d\n",i+1, fsi_set_cur->fsi_n);

		result[i] = fsi_set_cur;
	}

	return result;
}

/*
 *  Construct the L1 feature set.
 *
 *
 */
fsi_set_t* const_L1_apriori() {
	fsi_set_t *fsi_set_v;
	fsi_t *fsi_v;
	bst_node_t *bst_node_v;
	B_KEY_TYPE sup;

	fsi_set_v = alloc_fsi_set();

	// for each feature
	bst_node_v = bst_get_min(IF_v->root);
	while (bst_node_v != NULL) {
		fsi_v = alloc_fsi(1);
		fsi_v->feaset[0] = bst_node_v->key;

		if (cost_tag == 2) {
			// L1 is number of objects with the feature
//			sup = bst_node_v->p_list_obj->obj_n;
			sup = bst_node_v->totalWeight;
			sup = sup / fea_highest_freq; // normalization
		} else {
			sup = 1; // by definition
		}

		fsi_v->sup = sup;

		if (sup >= min_sup)
			add_fsi_set_entry(fsi_set_v, fsi_v);
		else
			release_fsi(fsi_v);

		bst_node_v = bst_successor(bst_node_v);
	}

	return fsi_set_v;
}

/*
 *	The checking for the join step.
 *
 *	Check whether two sets @fsi_v1 and @fsi_v2
 *	share s-1 features, where s is the number of features in @fsi_v1 or
 *@fsi_v2.
 *
 *	return  -1 if the join check fails,
 *	otherwise, return the (only) feature in fsi_v2 that is not contained in
 *fsi_v1.
 *
 */
FEA_TYPE join_check(fsi_t *fsi_v1, fsi_t *fsi_v2) {
	FEA_TYPE fea_v;

	fea_v = -1;

	if (fsi_v1->fea_n != fsi_v2->fea_n)
		return -1;

	if (memcmp(fsi_v1->feaset, fsi_v2->feaset,
			(fsi_v1->fea_n - 1) * sizeof(FEA_TYPE)) == 0) {
		fea_v = fsi_v2->feaset[fsi_v2->fea_n - 1];
	}

	return fea_v;
}

/*
 *  The checking for the prune step.
 *
 *  return true if all subsets of @fsi_v appear in @fsi_set_v.
 *
 *  Method:
 *  ignore one of the feature one by one to form the k-1 size feature set
 *  thus, number of k-1 size feature set = k
 *
 *  for each such feature set, we check whether it is in @fsi_set_v
 *  if not, return false
 */
bool all_subesets_exist(fsi_set_t *fsi_set_v, fsi_t *fsi_v) {
	int i, x, y;
	fsi_t *fsi_temp;

	// i is the pos to be igore to have k-1 size
	for (i = 0; i < fsi_v->fea_n; i++) {
		// check whether fsi_v\[i] exist in fsi_set_v
		fsi_temp = fsi_set_v->head->next;
		while (fsi_temp != NULL) {
			// check whether fsi_temp == fsi_v\[i]
			for (x = 0; x < fsi_v->fea_n; x++) {
				if (x == i)
					continue;
				for (y = 0; y < fsi_temp->fea_n; y++) {
					if (fsi_temp->feaset[y] == fsi_v->feaset[x])
						break;
				}
				if (y == fsi_temp->fea_n) { // cannot find
					break;
				}
			}
			if (x == fsi_v->fea_n) // fsi_temp matches fsi_v\[i]
					{
				break;
			}
			fsi_temp = fsi_temp->next;
		}
		if (fsi_temp == NULL) // no matches
		{
			return false;
		}
	}

	return true;
}

//---------------------------------------------------------

/*
 cost_tag
 1: participation
 2: fraction-score
 3: partitioning
 4: construdction
 5: enumeration
 */
B_KEY_TYPE comp_support(int alg_opt, fsi_t *fsi_v, int numOfObj,
		unordered_set<int> *RI1, unordered_set<int> *RI2) {
	switch (cost_tag) {
	case 1: //participation-based
	case 2: //fraction-score
		return comp_support_fraction(alg_opt, fsi_v, numOfObj, RI1, RI2);
	case 3: //partitioning-based
		return comp_support_partitioning(alg_opt, fsi_v, numOfObj);
	case 4: //construction-based
		return comp_support_construction(alg_opt, fsi_v, numOfObj);
	case 5: //enumeration-based
		return comp_support_enum(alg_opt, fsi_v, numOfObj);
	default:
		return 0;
	}
}

/*
 * compute the support using partitioning-based approach
 * partition the space into grids implicitiy
 *
 * Method: use the objects with first feature to locate possible cells
 * count the number of objects of each feature in the cell
 * use the min one (possibly 0) as the "number of transactions" in the cell
 * sup+= num. of transactions
 *
 */
B_KEY_TYPE comp_support_partitioning(int alg_opt, fsi_t *fsi_v, int numOfObj) {
	B_KEY_TYPE sup;
	obj_node_t *obj_node_v, *obj_node_cur;
	obj_set_t *obj_set_v, *obj_set_cur;

	sup = 0;
	// use the first feature for finding grid cell
	obj_set_v = ((bst_node_t*) bst_search(IF_v, fsi_v->feaset[0]))->p_list_obj;
	obj_node_v = obj_set_v->head->next;
	// for each object o with the feature f
	while (obj_node_v != NULL) {

		//find the grid x y
		int x = floor(obj_node_v->obj_v->MBR[0].min / dist_thr);
		int y = floor(obj_node_v->obj_v->MBR[1].min / dist_thr);

		//------
		//we need to skip previous counted cell to avoid double count
		bool countedCellFlag = false;
		obj_node_t *obj_node_v2 = obj_set_v->head->next;
		while (obj_node_v2 != obj_node_v) {
			int x2 = floor(obj_node_v2->obj_v->MBR[0].min / dist_thr);
			int y2 = floor(obj_node_v2->obj_v->MBR[1].min / dist_thr);

			if (x == x2 && y == y2) {
				countedCellFlag = true;
				break;
			}
			obj_node_v2 = obj_node_v2->next;
		}
		if (countedCellFlag) {
			obj_node_v = obj_node_v->next;
			continue;
		}
		//------

		//		printf("x:%d\ty:%d\t, obj_v->x:%lf\t y:%lf \n", x,y,obj_node_v->obj_v->MBR[0].min ,obj_node_v->obj_v->MBR[1].min );
		//
		//count the number of objects of each feature in the cell
		//use the min one (possibly 0) as the support
		double cnt_min = INFINITY;
		double weight_min = INFINITY;

		//for each feature
		for (int i = 0; i < fsi_v->fea_n; i++) {
			int cnt = 0;
			obj_set_cur =
					((bst_node_t*) bst_search(IF_v, fsi_v->feaset[i]))->p_list_obj;
			obj_node_cur = obj_set_cur->head->next;
			while (obj_node_cur != NULL) {
				int x2 = floor(obj_node_cur->obj_v->MBR[0].min / dist_thr);
				int y2 = floor(obj_node_cur->obj_v->MBR[1].min / dist_thr);

				if (x == x2 && y == y2) {
					cnt++;
					if (obj_node_cur->obj_v->weight < weight_min)
						weight_min = obj_node_cur->obj_v->weight;
				}
				obj_node_cur = obj_node_cur->next;
			}
			if (cnt < cnt_min)
				cnt_min = cnt;
		}		//end for each feature
		sup += cnt_min * weight_min;

		obj_node_v = obj_node_v->next;
	}		//end for each grid

	sup = sup / (double) fea_highest_freq;
	fsi_v->sup = sup;

	return sup;
}

/*
 * compute the support using construction-based approach
 * Method: arbitrary construct instance of the feautre set
 * until no more instance can be constructed
 * count the number of such instances
 *
 * following the KDD01 paper
 */
B_KEY_TYPE comp_support_construction(int alg_opt, fsi_t *fsi_v, int numOfObj) {
	B_KEY_TYPE sup, sup2;
	obj_node_t *obj_node_v, *obj_node_cur;
	obj_set_t *obj_set_v, *obj_set_cur;

	sup = INFINITY;
	// use the objects in k-th feature as pivot
	for (int k = 0; k < fsi_v->fea_n; k++) {
		//RI keep track of *used* objects
		bool *RI = new bool[numOfObj]();
		for (int i = 0; i < numOfObj; i++)
			RI[i] = false;

		sup2 = 0;
		obj_set_v =
				((bst_node_t*) bst_search(IF_v, fsi_v->feaset[k]))->p_list_obj;
		obj_node_v = obj_set_v->head->next;
		// for each object o with the feature f
		while (obj_node_v != NULL) {
			obj_set_t *S_0 = alloc_obj_set();
			add_obj_set_entry(obj_node_v->obj_v, S_0);

			//for each remaining feature
			for (int i = 0; i < fsi_v->fea_n; i++) {
				if (i == k)
					continue;
				//---
				//calculate new center
				double x = 0, y = 0;
				obj_node_t *obj_node_v2 = S_0->head->next;
				while (obj_node_v2 != NULL) {
					x += obj_node_v2->obj_v->MBR[0].min;
					y += obj_node_v2->obj_v->MBR[1].min;
					obj_node_v2 = obj_node_v2->next;
				}
				x = x / S_0->obj_n;
				y = y / S_0->obj_n;
				//---

				obj_set_cur =
						((bst_node_t*) bst_search(IF_v, fsi_v->feaset[i]))->p_list_obj;
				obj_node_cur = obj_set_cur->head->next;
				while (obj_node_cur != NULL) {
					if (!RI[obj_node_cur->obj_v->id - 1]) {
						//						if(check_dist_constraint(S_0, obj_node_cur->obj_v, obj_node_v->obj_v, dist_thr)) {
						if (sqrt(
								pow(obj_node_cur->obj_v->MBR[0].min - x, 2)
										+ pow(
												obj_node_cur->obj_v->MBR[1].min
														- y, 2))
								<= dist_thr / 2) {
							add_obj_set_entry(obj_node_cur->obj_v, S_0);
							break;
						}
					}
					obj_node_cur = obj_node_cur->next;
				}
			}
			//if S_0 contains all features,
			//            if (S_0->obj_n == fsi_v->fea_n - 1) {
			if (S_0->obj_n == fsi_v->fea_n) {
				double min_weight = obj_node_v->obj_v->weight;
				obj_node_cur = S_0->head->next;
				while (obj_node_cur != NULL) {
					RI[obj_node_cur->obj_v->id - 1] = true;
					//--
					if (obj_node_cur->obj_v->weight < min_weight) {
						min_weight = obj_node_cur->obj_v->weight;
					}
					//--
					obj_node_cur = obj_node_cur->next;
				}
				RI[obj_node_v->obj_v->id - 1] = true;
//				sup2++;
//				sup2+=obj_node_v->obj_v->weight;
				sup2 += min_weight;
			}
			release_obj_set(S_0);

			obj_node_v = obj_node_v->next;
		}			//end for each obj with k-th feature

		delete[] RI;
		//    printf("cnt:%d\t sup2:%lf\n",obj_set_v->obj_n,sup2);

		if (sup2 < sup)
			sup = sup2;
	}			//end for each feature

	sup = sup / fea_highest_freq;
	fsi_v->sup = sup;

	return sup;
}

/*
 * compute the support using enumeration-based approach
 * construct *all* row instances and count the number
 *
 * Not used to compute support
 */
B_KEY_TYPE comp_support_enum(int alg_opt, fsi_t *fsi_v, int numOfObj) {
	B_KEY_TYPE sup, sup_C_f;
	obj_node_t *obj_node_v;
	obj_set_t *obj_set_v;

	bst_t *inverted_file;

	inverted_file = bst_ini();
	for (int i = 1; i < fsi_v->fea_n; i++) {
		bst_node_t *bst_node_v;
		// Insert a new keyword in the binary tree.
		bst_node_v = (bst_node_t*) malloc(sizeof(bst_node_t));
		memset(bst_node_v, 0, sizeof(bst_node_t));

		/*s*/
		stat_v.memory_v += sizeof(bst_node_t);
		if (stat_v.memory_v > stat_v.memory_max)
			stat_v.memory_max = stat_v.memory_v;
		/*s*/

		// Update the posting list.
		bst_node_v->key = fsi_v->feaset[i];
		bst_node_v->p_list_obj = copy_obj_set(
				((bst_node_t*) bst_search(IF_v, fsi_v->feaset[i]))->p_list_obj);
		printf("%d \t", bst_node_v->p_list_obj->obj_n);
		bst_insert(inverted_file, bst_node_v);
	}
	int total = 0;
	//arbitrary use the first feature for finding row instances
	obj_set_v = ((bst_node_t*) bst_search(IF_v, fsi_v->feaset[0]))->p_list_obj;
	obj_node_v = obj_set_v->head->next;
	// for each object o with the feature f
	while (obj_node_v != NULL) {

		// Initialize the S_0.
		obj_set_t *S_0 = alloc_obj_set();
		int cnt = 0;
		enum_sub(inverted_file, S_0, obj_node_v->obj_v, cnt);
		//		printf("cnt:%d\n",cnt);
		total += cnt;
		release_obj_set(S_0);

		obj_node_v = obj_node_v->next;
	}

	printf("%d\n", obj_set_v->obj_n);

	printf("total:%d\n", total);

	bst_release(inverted_file);

	sup = total;
	fsi_v->sup = sup;

	return sup;
}

void enum_sub(bst_t *IF_v, obj_set_t *S_0, obj_t *o, int &cnt) {
	obj_t *obj_v;
	bst_node_t *bst_node_v;
	obj_node_t *obj_node_v;
	bst_node_list_t *bst_node_list_v;

	if (IF_v->node_n == 0) {
		cnt++;
		return; // An empty one.
	}
	bst_node_v = IF_v->root;
	obj_node_v = bst_node_v->p_list_obj->head->next;
	while (obj_node_v != NULL) {
		// Pick an object.
		obj_v = obj_node_v->obj_v;

		// Distance constraint checking.
		if (!check_dist_constraint(S_0, obj_v, o, dist_thr)) {
			obj_node_v = obj_node_v->next;
			continue;
		}

		// Update the IF_v.
		bst_node_list_v = update_IF_obj(IF_v, obj_v);

		// Update the S_0.
		// obj_v is added at the first place of S_0.
		add_obj_set_entry(obj_v, S_0);

		// Sub-procedure.
		enum_sub(IF_v, S_0, o, cnt);

		// Restore the S_0.
		remove_obj_set_entry(S_0);

		// Restore the IF_v.
		restore_IF_bst_node_list(IF_v, bst_node_list_v);

		release_bst_node_list(bst_node_list_v);

		// Try the next object candidate.
		obj_node_v = obj_node_v->next;
	}
	return;
}

/*
 * The implementaiton of "SupportComputation" in the paper.
 *  @fsi_v = current feature set C
 *
 *  Output: @sup = sup(C)
 */
B_KEY_TYPE comp_support_fraction(int alg_opt, fsi_t *fsi_v, int numOfObj,
		unordered_set<int> *RI1, unordered_set<int> *RI2) {
	B_KEY_TYPE sup, sup_C_f;
	obj_node_t *obj_node_v;
	obj_set_t *obj_set_v;

	sup = INFINITY;
	sup_C_f = 0;

	double sumWeight = 0.0;

	//    FEA_TYPE fea_worst;
	//    obj_set_t* obj_set_worst;

	// keep track objects that are involved in row instances that have been found
	// note that obj_id-1 is needed inside the array
//	bool *RI = new bool[numOfObj]();
//	unordered_set<int>* RIno = new bool[numOfObj]();
	unordered_set<int> *RI = new unordered_set<int>();
	unordered_set<int> *RIno = new unordered_set<int>();

	// for each feature f in C
	for (int i = 0; i < fsi_v->fea_n; i++) {
		sup_C_f = 0;

		// the corresponding inverted list in IF
		obj_set_v =
				((bst_node_t*) bst_search(IF_v, fsi_v->feaset[i]))->p_list_obj;

		//--------------------
		//assume all objects have RI
		int cnt = obj_set_v->obj_n;
		double *frac_receive = new double[cnt];
		double frac_receive_UB = 0.0;
		int j = 0;
		if (cost_tag == 2 && fsi_v->fea_n > 2) {
			obj_node_v = obj_set_v->head->next;
			// for each object o with the feature f
			while (obj_node_v != NULL) {
				frac_receive[j] = min_frac_receive(fsi_v, obj_node_v->obj_v);
				frac_receive_UB += frac_receive[j];
				obj_node_v = obj_node_v->next;
				j++;
			}
			if (frac_receive_UB / fea_highest_freq < min_sup) {
//				print_fsi(fsi_v,stdout);
//				printf("frac_receive_UB: %.3lf\n",
//						frac_receive_UB / fea_highest_freq);
				sup = frac_receive_UB;
				free(frac_receive);
				break; //end this label set, since UB cannot >minsup!
			}
		}
		//--------------------
		obj_node_v = obj_set_v->head->next;
		j = 0;
		sumWeight = 0.0;
		// for each object o with the feature f
		while (obj_node_v != NULL) {
			int flag = 0;
			if (cost_tag == 2) {
				flag = filters(alg_opt, fsi_v, obj_node_v->obj_v, i, RI, RIno,
						RI1, RI2);
			}
			if (flag == 1
					|| (flag == 0
							&& check_row_instance(alg_opt, fsi_v,
									obj_node_v->obj_v, RI, RIno))) {
//				RI[obj_node_v->obj_v->id - 1] = 1;
				RI->insert(obj_node_v->obj_v->id);
				if (cost_tag == 1) {
//					sup_C_f += 1; //each group is counted as 1
					sup_C_f += obj_node_v->obj_v->weight;
				} else if (cost_tag == 2) {
					if (fsi_v->fea_n == 2)
						sup_C_f += min_frac_receive(fsi_v, obj_node_v->obj_v);
					else
						sup_C_f += frac_receive[j];
					//early stopping: sup grouped by this fea is higher than by other
					if (sup_C_f > sup)
						break;
				}
			}
			//---
			if (cost_tag == 2 && fsi_v->fea_n > 2) {
				//early stopping: assume all remaining objects contribute
				frac_receive_UB -= frac_receive[j];
				if ((sup_C_f + frac_receive_UB) / fea_highest_freq < min_sup) {
					break;
				}
			}
			//---
			j++;
			sumWeight += obj_node_v->obj_v->weight;

			obj_node_v = obj_node_v->next;
		}			//end loop each object
		//if participation.
		if (cost_tag == 1)
//			sup_C_f = sup_C_f / (double) obj_set_v->obj_n;
			sup_C_f = sup_C_f / (double) sumWeight;

		// maintain the smallest one here
		if (sup_C_f <= sup) {
			sup = sup_C_f;

			// early stopping: sup of this fea set < threshold
			if ((cost_tag == 1 && sup < min_sup)
					|| (cost_tag == 2 && sup / fea_highest_freq < min_sup)) {
				free(frac_receive);
				break;
			}
		}
		free(frac_receive);
	}			//end for each label

	if (cost_tag == 2)
		sup = sup / fea_highest_freq;
	fsi_v->sup = sup;

	/*s*/
	stat_v.memory_v += RI->size() * sizeof(int);
	if (alg_opt < 60)
		stat_v.memory_v += RIno->size() * sizeof(int);
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	//new filter 3
	if (alg_opt == 46)
		fsi_v->RI = RI;
	else {
		stat_v.memory_v -= RI->size() * sizeof(int);
		RI->clear();
		RI = NULL;
	}
	//new verification step 2
//	delete[] RIno;
//	RIno = NULL;

	/*s*/
	if (alg_opt < 60) {
		stat_v.memory_v -= RIno->size() * sizeof(int);
	}
	/*s*/

	RIno->clear();
	RIno = NULL;

//	if (fsi_v->RI != NULL) {
//		int cnt2 = 0;
//		for (int j = 0; j < numOfObj; j++) {
//			if (fsi_v->RI[j] != 0)
//				cnt2++;
//		}
//		printf("cnt2:%d\n", cnt2);
//	}
	return sup;
}

/*
 * @return 	1 = true
 * 			0 = dont know
 * 			-1 = skip / continue
 */
int filters(int alg_opt, fsi_t *fsi_v, obj_t *obj_v, int i,
		unordered_set<int> *RI, unordered_set<int> *RIno,
		unordered_set<int> *RI1, unordered_set<int> *RI2) {
	int flag = 0;

	// ---- F1 ---- //new
	//(only) size 2 can use fraction to infer involve in a row instance
#ifndef WIN32
	float sys_t, S0_t = 0;
	struct rusage S0_sta, S0_end;
	GetCurTime(&S0_sta);
#endif
	if (alg_opt < 60 && cost_tag == 2)
		if (fsi_v->fea_n == 2) {
			int j = i == 0 ? 1 : 0;
			if (obj_v->frac_f[fsi_v->feaset[j] - 1] != 0) {
				flag = 1;
			} else {
				flag = -1;
			}
			//all size 2 are filtered by this step actually
			stat_v.S0_sum++;
#ifndef WIN32
			GetCurTime(&S0_end);
			GetTime(&S0_sta, &S0_end, &S0_t, &sys_t);
			stat_v.S0_time += S0_t;
#endif
			return flag;
		}
#ifndef WIN32
	GetCurTime(&S0_end);
	GetTime(&S0_sta, &S0_end, &S0_t, &sys_t);
	stat_v.S0_time += S0_t;
#endif
	//----F2----- (original F1)
	//i.e., previous found object
#ifndef WIN32
	sys_t = 0;
	float S3_t = 0;
	struct rusage S3_sta, S3_end;
	GetCurTime(&S3_sta);
#endif

//				flag = RI[obj_node_v->obj_v->id - 1];
	if (RI->find(obj_v->id) != RI->end()) {
		flag = 1;
		stat_v.S3_sum++;

#ifndef WIN32
		GetCurTime(&S3_end);
		GetTime(&S3_sta, &S3_end, &S3_t, &sys_t);
		stat_v.S3_time += S3_t;
#endif
		return flag;
	}

#ifndef WIN32
	GetCurTime(&S3_end);
	GetTime(&S3_sta, &S3_end, &S3_t, &sys_t);
	stat_v.S3_time += S3_t;
#endif
	//----------F2'----------- //new
	if (alg_opt == 46) {
		//super set checking
		//if found then yes then return yes =1
		if (RI1 != NULL && RI1->find(obj_v->id) != RI1->end()) {
			flag = 1;
			return flag;
		}
	}
	return 0; //dont know
}

/*
 *  Perform check row instance operation
 *  @alg_opt:
 *  1 = combinational
 *  2 = dia-CoSKQ-adapt
 *  3 = optimization-search
 *  41 = filter and verify + 1
 *  42 = filter and verify + 2
 *  43 = filter and verify + 3
 *
 *  @fsi_v: the current checking feature set C
 *  @obj_v: the current checking object
 */
bool check_row_instance(int alg_opt, fsi_t *fsi_v, obj_t *obj_v,
		unordered_set<int> *RI, unordered_set<int> *RIno) {
#ifndef WIN32
	float sys_t, S4_t = 0;
	struct rusage S4_sta, S4_end;
	GetCurTime(&S4_sta);
#endif
	bool flag = true;
	if (alg_opt == 1) { // combinatorial
		flag = combinatorial(fsi_v, obj_v, NULL, RI, RIno);
	} else if (alg_opt == 2) { // Dia
		flag = dia(fsi_v, obj_v, NULL, RI, RIno);
	} else if (alg_opt == 3) { // mCK
		flag = mck(fsi_v, obj_v, NULL, RI, RIno);
	} else if (alg_opt == 41 || alg_opt == 42 || alg_opt == 43) { // filter-and-verification
		flag = filter_and_verify(alg_opt - 40, fsi_v, obj_v, RI, RIno);
	} else if (alg_opt == 46) { // filter-and-verification
		flag = filter_and_verify(alg_opt - 45, fsi_v, obj_v, RI, RIno);
	} else if (alg_opt == 61) { // filter-and-verification
		flag = filter_and_verify(1, fsi_v, obj_v, RI, RIno);
	} else {
		printf("alg_opt error!\n");
		exit(-1);
	}
#ifndef WIN32
	GetCurTime(&S4_end);
	GetTime(&S4_sta, &S4_end, &S4_t, &sys_t);
	stat_v.S4_time += S4_t;
#endif

//---
//	RIno[obj_v->id - 1] = true;
	if (alg_opt < 60 && flag == false) //explicit if to avoid confusion here.
		RIno->insert(obj_v->id);
//---
	return flag;
}

/*
 *  Method 1 : Combinatorial Approach
 *
 *  if we call this function from method 4, S2 = range_query(disk_v,q) to save
 * some time
 *  else S2 = NULL
 */
bool combinatorial(fsi_t *fsi_v, obj_t *obj_v, obj_set_t *S2,
		unordered_set<int> *RI, unordered_set<int> *RIno) {
	bst_t *inverted_list;
	obj_set_t *S_0, *S = NULL;
	obj_set_t *obj_set_v;
	loc_t *loc_v;
	disk_t *disk_v;
	query_t *q;

// Step 1. Range query.

	loc_v = get_obj_loc(obj_v);
	disk_v = alloc_disk(IRTree_v.dim);
	set_disk(disk_v, loc_v, dist_thr);

	q = alloc_query();
	q->loc_v = loc_v;
	q->psi_v = fsi_to_psi(fsi_v, obj_v->fea);

	if (S2 == NULL) {
		obj_set_v = range_query(disk_v, q);

		// Step 2. Construct an IF structure.
		inverted_list = const_IF(obj_set_v, q->psi_v, RIno);
//		inverted_list = const_IF(obj_set_v, q->psi_v);

		release_obj_set(obj_set_v);
		release_disk(disk_v);
		release_query(q);

		// return false if any key is missing
		if (!bst_check_plist_obj_n(inverted_list->root)) {
			release_IF(inverted_list);
			return false;
		}
	} else {
		inverted_list = const_IF(S2, q->psi_v, RIno);
//		inverted_list = const_IF(S2, q->psi_v);
		release_obj_set(S2);
		release_disk(disk_v);
		release_query(q);
	}
// Initialize the S_0.
	S_0 = alloc_obj_set();

// Invoke the sub-procedure "recursively".
	S = combinatorial_sub(inverted_list, S_0, obj_v, dist_thr);

// Release the resources.
	release_obj_set(S_0);
	release_IF(inverted_list);
	if (S != NULL && S->obj_n != 0) {
		//        add_obj_set_entry(obj_v, S);
		//        printf("S5:%f\t%f\n",comp_diameter(S), dist_thr);
		//---
		// update RI
		{
#ifndef WIN32
			float sys_t, S3_t = 0;
			struct rusage S3_sta, S3_end;
			GetCurTime(&S3_sta);
#endif
			obj_node_t *obj_node_v;
			obj_node_v = S->head->next;
			while (obj_node_v != NULL) {
//				RI[obj_node_v->obj_v->id - 1] = true;
				RI->insert(obj_node_v->obj_v->id);
				obj_node_v = obj_node_v->next;
			}

#ifndef WIN32
			GetCurTime(&S3_end);
			GetTime(&S3_sta, &S3_end, &S3_t, &sys_t);
			stat_v.S3_time += S3_t;
#endif
		}
		//---
		release_obj_set(S);
		return true;
	}

	release_obj_set(S);
	return false;
}

/*
 *	The sub-procedure of combinatorial approach.
 *
 * we process as follows: for each object o in the root,
 *  we update the IF s.t. the keywords in o are removed (note that the root must be removed), then o is added to S_0, then we recursive call the function to find another object in the (new) root and update IF and add to S_0...., until the IF is empty, then S covers all the keywords.
 * note that if we do not consider the distance constraint, we can always find a feaisble set (e.g., by adding all objects in IF into S). When we consider the distance constraint, we may return NULL.
 */
obj_set_t* combinatorial_sub(bst_t *IF_v, obj_set_t *S_0, obj_t *o,
B_KEY_TYPE d) {
	obj_t *obj_v;
	obj_set_t *S;
	bst_node_t *bst_node_v;
	obj_node_t *obj_node_v;
	bst_node_list_t *bst_node_list_v;

	if (IF_v->node_n == 0)
		return alloc_obj_set(); // An empty one.

	bst_node_v = IF_v->root;
	obj_node_v = bst_node_v->p_list_obj->head->next;
	while (obj_node_v != NULL) {
		// Pick an object.
		obj_v = obj_node_v->obj_v;

		// Distance constraint checking.
		if (!check_dist_constraint(S_0, obj_v, o, d)) {
			obj_node_v = obj_node_v->next;
			continue;
		}

		// Update the IF_v.
		bst_node_list_v = update_IF_obj(IF_v, obj_v);

		// Update the S_0.
		// obj_v is added at the first place of S_0.
		add_obj_set_entry(obj_v, S_0);

		// Sub-procedure.
		S = combinatorial_sub(IF_v, S_0, o, d);

		// Restore the S_0.
		remove_obj_set_entry(S_0);

		// Restore the IF_v.
		restore_IF_bst_node_list(IF_v, bst_node_list_v);

		release_bst_node_list(bst_node_list_v);

		// Checking.
		if (S != NULL) {
			// Include obj_v into S.
			add_obj_set_entry(obj_v, S);

			return S;
		}

		// S == NULL.

		// Try the next object candidate.
		obj_node_v = obj_node_v->next;
	}
	return NULL;
}

/*
 * method for inverted file
 * return true if all nodes have p_list->obj_n > 0
 */
bool bst_check_plist_obj_n(bst_node_t *bst_node_v) {
	if (bst_node_v == NULL)
		return true;

	return (bst_node_v->p_list_obj->obj_n > 0)
			&& bst_check_plist_obj_n(bst_node_v->left)
			&& bst_check_plist_obj_n(bst_node_v->right);
}

/*
 * method for inverted file
 * return true if all nodes of fsi_v->feaset[i] exist and have p_list->obj_n > 0
 */
bool bst_check_plist(bst_t *bst_v, fsi_t *fsi_v, FEA_TYPE fea) {
	bst_node_t *bst_node_v;

	if (bst_v == NULL)
		return false;

	for (int i = 0; i < fsi_v->fea_n; i++) {
		if (fsi_v->feaset[i] == fea)
			continue;

		bst_node_v = bst_search(bst_v, fsi_v->feaset[i]);
		if (bst_node_v == NULL || bst_node_v->p_list_obj->obj_n == 0) {
			return false;
		}
	}
	return true;
}

/*
 *  Method 2 : Dia-CoSKQ Approach
 *
 *  if we call this function from method 4, S2 = range_query(disk_v,q) to save
 * some time
 *  else S2 = NULL
 *
 */
bool dia(fsi_t *fsi_v, obj_t *obj_v, obj_set_t *S2, unordered_set<int> *RI,
		unordered_set<int> *RIno) {
	query_t *q;
	obj_set_t *S;
	loc_t *loc_v;
	disk_t *disk_v;

// Constrcut Dia-CoSKQ instance.
	loc_v = get_obj_loc(obj_v);
	disk_v = alloc_disk(IRTree_v.dim);
	set_disk(disk_v, loc_v, dist_thr);

// q.loc = obj_v->loc
// q.psi = fsi_v - o.f
	q = alloc_query();
	q->loc_v = loc_v;
	q->psi_v = fsi_to_psi(fsi_v, obj_v->fea);

// invoke the Dia-Exact for the constructed instance
// note that the UB disk (and LB disk) is formed inside the Dia-Exact
	S = CostEnum(q, 1, 0, S2, disk_v);

	release_query(q);
	release_disk(disk_v);
	release_obj_set(S2);

	if (S != NULL) {
		add_obj_set_entry(obj_v, S);
		if (comp_diameter(S) <= dist_thr) { // update RI
			{
#ifndef WIN32
				float sys_t, S3_t = 0;
				struct rusage S3_sta, S3_end;
				GetCurTime(&S3_sta);
#endif
				obj_node_t *obj_node_v;
				obj_node_v = S->head->next;
				while (obj_node_v != NULL) {
//					RI[obj_node_v->obj_v->id - 1] = true;
					RI->insert(obj_node_v->obj_v->id);
					obj_node_v = obj_node_v->next;
				}

#ifndef WIN32
				GetCurTime(&S3_end);
				GetTime(&S3_sta, &S3_end, &S3_t, &sys_t);
				stat_v.S3_time += S3_t;
#endif
			}
			release_obj_set(S);
			return true;
		}
	}

	release_obj_set(S);
	return false;
}

/*
 *  Method 3 : mCK Approach
 *
 *  if we call this function from method 4, S2 = range_query(disk_v,q) to save
 * some time
 *  else S2 = NULL
 *
 */
bool mck(fsi_t *fsi_v, obj_t *obj_v, obj_set_t *S2, unordered_set<int> *RI,
		unordered_set<int> *RIno) {
	query_t *q;
	obj_set_t *S = NULL;
	loc_t *loc_v;
	disk_t *disk_v;
	obj_set_t *obj_set_v;

// Constrcut mCK instance.
	loc_v = get_obj_loc(obj_v);
	disk_v = alloc_disk(IRTree_v.dim);
	set_disk(disk_v, loc_v, dist_thr);

	q = alloc_query();
	q->loc_v = loc_v;
	q->psi_v = fsi_to_psi(fsi_v, obj_v->fea);

	if (S2 == NULL)
		// retrieve all objects within the disk of dist_thr
		obj_set_v = range_query(disk_v, q);
	else
		obj_set_v = S2;

// we pass the objects inside the disk (formed by disk_thr, q)
// to restict the search space of the method
	if (check_feasibility(obj_set_v, q->psi_v)) {
		S = mck_Exact(q, 0.01, obj_set_v, disk_v);
	}

	release_query(q);
	release_disk(disk_v);
	release_obj_set(obj_set_v);

	if (S != NULL) {
		add_obj_set_entry(obj_v, S);
		if (comp_diameter(S) <= dist_thr) {
			// update RI
			{
#ifndef WIN32
				float sys_t, S3_t = 0;
				struct rusage S3_sta, S3_end;
				GetCurTime(&S3_sta);
#endif
				obj_node_t *obj_node_v;
				obj_node_v = S->head->next;
				while (obj_node_v != NULL) {
//					RI[obj_node_v->obj_v->id - 1] = true;
					RI->insert(obj_node_v->obj_v->id);
					obj_node_v = obj_node_v->next;
				}

#ifndef WIN32
				GetCurTime(&S3_end);
				GetTime(&S3_sta, &S3_end, &S3_t, &sys_t);
				stat_v.S3_time += S3_t;
#endif
			}
			release_obj_set(S);
			return true;
		}
	}

	release_obj_set(S);
	return false;
}

/*
 *  Method 4 : Filtering-and-Verification Approach
 *
 *  @alg_opt2:
 *  1 = (method 1) combinatorial
 *  2 = (method 2) Dia
 *  3 = (method 3) mCK
 */
bool filter_and_verify(int alg_opt2, fsi_t *fsi_v, obj_t *obj_v,
		unordered_set<int> *RI, unordered_set<int> *RIno) {
	loc_t *loc_v;
	obj_set_t *obj_set_v = NULL; // range query Disk(o,d) or NULL
	query_t *q = NULL;
	int temp = -1; // 1:true,0:false,-1:null
	bool feasibleflag1 = false; // true = feasible, false = not feasible
	bool feasibleflag2 = false;
	bool feasibleflag3 = false;

//---------------------------------------------------
// Filter 3. feasiblility check in D(o,d) by N_o_f

#ifndef WIN32
	float sys_t, S1_t = 0;
	struct rusage S1_sta, S1_end;
	GetCurTime(&S1_sta);
#endif
	if (cost_tag == 2) {
		feasibleflag1 = check_Nof_feasibility(fsi_v, obj_v);
#ifndef WIN32
		GetCurTime(&S1_end);
		GetTime(&S1_sta, &S1_end, &S1_t, &sys_t);
		stat_v.S1_time += S1_t;
#endif

		if (!feasibleflag1) {
			stat_v.S1_sum++;
			return false;
		}

		if (fsi_v->fea_n == 2 && feasibleflag1) {
			stat_v.S1_sum++;
			return true;
		}
	}
//---------------------------------------------------
// Filter 4. range query in D(o,d/2)
	loc_v = get_obj_loc(obj_v);
	q = alloc_query();
	q->loc_v = loc_v;
	q->psi_v = fsi_to_psi(fsi_v, obj_v->fea);
	if (cost_tag == 2) {
#ifndef WIN32
		float S2_t = 0;
		struct rusage S2_sta, S2_end;
		GetCurTime(&S2_sta);
#endif

		feasibleflag2 = check_Nof2_feasibility(fsi_v, obj_v);

#ifndef WIN32
		GetCurTime(&S2_end);
		GetTime(&S2_sta, &S2_end, &S2_t, &sys_t);
		stat_v.S2_time += S2_t;
#endif

		if (feasibleflag2) {
			stat_v.S2_sum++;
			release_query(q);
			return true;
		}
	}

//---------------------------------------------------
// verification

	stat_v.S5_sum++;

//---------------------------------------------------

#ifndef WIN32
	float S5_t = 0;
	struct rusage S5_sta, S5_end;
	GetCurTime(&S5_sta);
#endif
// note that obj_set_v is released inside the functions
	if (alg_opt2 == 1)
		feasibleflag3 = combinatorial(fsi_v, obj_v, obj_set_v, RI, RIno);
	else if (alg_opt2 == 2)
		feasibleflag3 = dia(fsi_v, obj_v, obj_set_v, RI, RIno);
	else if (alg_opt2 == 3)
		feasibleflag3 = mck(fsi_v, obj_v, obj_set_v, RI, RIno);
#ifndef WIN32
	GetCurTime(&S5_end);
	GetTime(&S5_sta, &S5_end, &S5_t, &sys_t);
	stat_v.S5_time += S5_t;
#endif

	release_query(q);
	return feasibleflag3;
}

/*
 * Check whether all fea in @fsi_v > 0 in obj_v->N_o_f
 * @return
 * true: yes
 * false: no
 */
bool check_Nof_feasibility(fsi_t *fsi_v, obj_t *obj_v) {
	for (int i = 0; i < fsi_v->fea_n; i++) {
		if (fsi_v->feaset[i] == obj_v->fea)
			continue;

//		if (obj_v->N_o_f[fsi_v->feaset[i] - 1] == 0) {
		if (obj_v->frac_f[fsi_v->feaset[i] - 1] == 0) {
			// not found
			return false;
		}
	}
	return true;
}

/*
 * Check whether all fea in @fsi_v > 0 in obj_v->N_o_f2
 * @return
 * true: yes
 * false: no
 */
bool check_Nof2_feasibility(fsi_t *fsi_v, obj_t *obj_v) {
	for (int i = 0; i < fsi_v->fea_n; i++) {
		if (fsi_v->feaset[i] == obj_v->fea)
			continue;

		if (obj_v->N_o_f2[fsi_v->feaset[i] - 1] != 1) {
			// not found
			return false;
		}
	}
	return true;
}

int pruning_NNset(query_t *q, obj_t *obj_v) {
	obj_set_t *obj_set_v;
	k_node_t *k_node_v;
	B_KEY_TYPE dia;

	obj_set_v = alloc_obj_set();
	k_node_v = q->psi_v->k_head->next;

	while (k_node_v != NULL) {
		add_obj_set_entry(const_NN_key(q->loc_v, k_node_v->key, NULL),
				obj_set_v);
		k_node_v = k_node_v->next;
	}

	if (obj_set_v != NULL) {
		add_obj_set_entry(obj_v, obj_set_v);
		dia = comp_diameter(obj_set_v);
		release_obj_set(obj_set_v);
		if (dia <= dist_thr)
			return 1;
		if (dia / 2 > dist_thr)
			return 0;
	}

	return -1;
}

//---------------------------------------------------------------

/*
 *  The implementation of "FractionAggregation" in the paper.
 *
 *  Input:
 *  @fsi_v = current feature set C
 *  @obj_v = an object o
 *
 *  Output: frac-fs(o|C)
 *
 *  new: computation based on obj_v->frac_f which is pre-computed
 *  return the minimum fraction obj_v receives wrt fsi_v
 */
B_KEY_TYPE min_frac_receive(fsi_t *fsi_v, obj_t *obj_v) {
	B_KEY_TYPE fs = INFINITY;

	for (int i = 0; i < fsi_v->fea_n; i++) {
		if (fsi_v->feaset[i] == obj_v->fea)
			continue;

		if (obj_v->frac_f[fsi_v->feaset[i] - 1] < fs)
			fs = obj_v->frac_f[fsi_v->feaset[i] - 1];
	}

	return fs;
}

//-----------------------------------------

/*
 * The implementation of "FractionComputation" in the paper.
 *
 *  Pre-computation of obj_v[o].N_o_f for all objects in data_v
 *
 */
void precomputation(data_t *data_v, B_KEY_TYPE dist_thr) {
	k_node_t *k_head, *k_node_v;
	disk_t *disk_v;
	loc_t *loc_v;
	obj_set_t *obj_set_v;
	obj_node_t *obj_node_v;
//    std::unordered_map<FEA_TYPE,int>* map, * map2;
//    std::unordered_map<FEA_TYPE,float> * tmp;
	int *map;
	short int *map2;
	float *tmp;

	k_head = collect_keywords_bst(IF_v);

// for each object
	for (int i = 0; i < data_v->obj_n; i++) {

		//--checking to avoid unecessary memory allocation
		FEA_TYPE fea = data_v->obj_v[i].fea;
		if (!fea_freq_check(fea))
			continue;

		//---
		data_v->obj_v[i].frac_f = new float[data_v->key_n]();
//		data_v->obj_v[i].N_o_f = new int[data_v->key_n]();
		data_v->obj_v[i].N_o_f2 = new short int[data_v->key_n]();

		/*s*/
		stat_v.memory_v += (sizeof(float) + sizeof(short int)) * data_v->key_n;
		/*s*/
	}

	/*s*/
//		stat_v.memory_v += (sizeof(float) * data_v->key_n) * data_v->obj_n;
//	//	stat_v.memory_v += (sizeof(int) * data_v->key_n) * data_v->obj_n * 2;
//		stat_v.memory_v += (sizeof(short int) * data_v->key_n) * data_v->obj_n;
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

//	printf("pre-compute size:%lu\n",(sizeof(float) * data_v->key_n) * data_v->obj_n+(sizeof(int) * data_v->key_n) * data_v->obj_n * 2);
//*we iterate each object that contain frequent labels
// for each object
	for (int i = 0; i < data_v->obj_n; i++) {

		FEA_TYPE fea = data_v->obj_v[i].fea;
		double weight = data_v->obj_v[i].weight;

		if (!fea_freq_check(fea))
			continue;

//		printf("precompute: %d\n",i);

//		map = data_v->obj_v[i].N_o_f;
		map = new int[data_v->key_n]();
		map2 = data_v->obj_v[i].N_o_f2;

		// 2. find all objects in range D(o,d)
		loc_v = get_obj_loc(&data_v->obj_v[i]);

		disk_v = alloc_disk(IRTree_v.dim);
		set_disk(disk_v, loc_v, dist_thr);

		obj_set_v = range_query(disk_v);

		// 3. loop each obj in range
		obj_node_v = obj_set_v->head->next;
		while (obj_node_v != NULL) {

			map[obj_node_v->obj_v->fea - 1]++;

			if (calc_dist_obj(obj_node_v->obj_v, &data_v->obj_v[i])
					<= dist_thr / 2.0)
				map2[obj_node_v->obj_v->fea - 1] = 1;

			obj_node_v = obj_node_v->next;
		}

		//-------------
		// assign fraction from this object to each object in range
		obj_node_v = obj_set_v->head->next;
		while (obj_node_v != NULL) {

			tmp = obj_node_v->obj_v->frac_f;
			if (tmp != NULL) {
				//if frac_f ==NULL, this obj's label is not freq

//			tmp[data_v->obj_v[i].fea - 1] += 1
//					/ ((double) map[obj_node_v->obj_v->fea - 1]);
				tmp[data_v->obj_v[i].fea - 1] += weight
						/ ((double) map[obj_node_v->obj_v->fea - 1]);

				if (tmp[data_v->obj_v[i].fea - 1] > 1)
					tmp[data_v->obj_v[i].fea - 1] = 1;
			}
			obj_node_v = obj_node_v->next;
		}

		//-------------
		release_disk(disk_v);
		release_obj_set(obj_set_v);
		release_loc(loc_v);
		delete (map);
	}

	release_k_list(k_head);
}

/**
 * return true if this fea has freq > threshold
 * i.e., need to be considered
 */
bool fea_freq_check(FEA_TYPE fea) {
	bst_node_t *bst_node_v = bst_search(IF_v, fea);
	if (bst_node_v != NULL) {
		double n = bst_node_v->totalWeight;
		if (n < fea_highest_freq * min_sup) {
			return false;
		} else
			return true;
	}
	return false;
}

/*
 *	Construct the IF on a set of objects @obj_set_v for the keywords in @psi_v.
 *
 *	1. The IF structure is indexed by a binary search tree.
 *	2. No ordering is imposed in IF.
 */
bst_t* const_IF(obj_set_t *obj_set_v, psi_t *psi_v, unordered_set<int> *RIno) {
	int i;
	bst_t *IF_v;
	k_node_t *k_node_v;
	obj_node_t *obj_node_v;
	bst_node_t *bst_node_v;

//IF_v = alloc_IF( psi_v->key_n);
	IF_v = bst_ini();

	k_node_v = psi_v->k_head->next;
//for( i=0; i<psi_v->key_n; i++)
	while (k_node_v != NULL) {
		//IF_v->entry_v[ i].key = k_node_v->key;
		bst_node_v = (bst_node_t*) malloc(sizeof(bst_node_t));
		memset(bst_node_v, 0, sizeof(bst_node_t));

		/*s*/
		stat_v.memory_v += sizeof(bst_node_t);
		/*s*/

		bst_node_v->key = k_node_v->key;
		bst_node_v->p_list_obj = alloc_obj_set();

		obj_node_v = obj_set_v->head->next;
		while (obj_node_v != NULL) {
//			if (RIno[obj_node_v->obj_v->id - 1] == 0 //not no found
			if ((RIno == NULL
					|| RIno->find(obj_node_v->obj_v->id) == RIno->end())
					&& has_key_obj(obj_node_v->obj_v, k_node_v->key))
				add_obj_set_entry(obj_node_v->obj_v, bst_node_v->p_list_obj);

			obj_node_v = obj_node_v->next;
		}

		bst_insert(IF_v, bst_node_v);

		k_node_v = k_node_v->next;
	}

	/*s*/
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return IF_v;
}

//============================ maximal =======================

/*
 *  [Maximal]
 */

fsi_set_t** maximal(int alg_opt, int numOfObj, int numOfFea) {
	fsi_set_t **maximalResult; // storing overall result
	fsi_set_t **seedResult; // storing overall result

// Initialize the structure for storing L_1, L_2, ..., L_{|F|}
	maximalResult = (fsi_set_t**) malloc(numOfFea * sizeof(fsi_set_t*));
	memset(maximalResult, 0, numOfFea * sizeof(fsi_set_t*));

	seedResult = (fsi_set_t**) malloc(numOfFea * sizeof(fsi_set_t*));
	memset(seedResult, 0, numOfFea * sizeof(fsi_set_t*));

	/*s*/
	stat_v.memory_v += (numOfFea + numOfFea) * sizeof(fsi_set_t*);
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
	seedResult[0] = const_L1_apriori();

#ifndef WIN32
	GetCurTime(&query_end);

	GetTime(&query_sta, &query_end, &usr_t, &sys_t);
	printf("L_1:%d\n", seedResult[0]->fsi_n);
	GetCurTime(&query_sta);
#endif
	fsi_set_t *fsi_set_v = NULL;
	fsi_set_t **cmp = NULL;
	int maxSize = 0;
	int k = 1; //k=1 = L_2
	while (seedResult[k - 1]->fsi_n > 0) {
		//Step 1
		double sup_avg = 0.0;
		maxSize = 0;
		seedResult[k] = const_L_k_apriori(alg_opt, numOfObj, seedResult[k - 1],
				k - 1, sup_avg);
		printf("L_%d:%d sup_avg:%.5lf\n", k + 1, seedResult[k]->fsi_n, sup_avg);

		///===========================================
		//Step 2
		fsi_set_v = NULL;
		if (seedResult[k] != NULL && seedResult[k]->fsi_n > 0) {

//			printf("generate graph from size %d!\n", k + 1);
			//gen graph from L_k
			adjacency_list.clear();
			adjacency_list = genGraph(numOfFea, seedResult[k]);

			//Step 2.2 find candidate maximal cliques
			set<vector<int>> maximalClqiues = run_BronKerbosch();

			fsi_set_v = candidate_maximal_patterns(maximalClqiues, maxSize);
			cmp = distribute_by_size(fsi_set_v, maxSize);
			release_fsi_set(fsi_set_v);


//			for (int i = 0; i < maxSize; i++) {
//				printf("CMP_%d:%d\n", i + 1, cmp[i]->fsi_n);
//			}

			//check breaking condition here.
			if (sup_avg == 0) {
				break;
			} else {
				for (int i = 0; i < maxSize; i++) {
					release_fsi_set(cmp[i]);
				}
			}

		}
		///===========================================
		k++;
	}

#ifndef WIN32
	GetCurTime(&query_end);

	GetTime(&query_sta, &query_end, &usr_t, &sys_t);
	printf("P1 \ttime:%0.5lf\n", usr_t);

	GetCurTime(&query_sta);
#endif

//	printf("S%d\n", 3);
//	print_fsi_set(fsi_set_v,stdout);

//==================================
//Step 3 compute sup for each candidate
	find_maximal(k, maxSize, numOfObj, cmp, maximalResult);
//============================================

#ifndef WIN32
	GetCurTime(&query_end);

	GetTime(&query_sta, &query_end, &usr_t, &sys_t);
	printf("P2 \ttime:%0.5lf\n", usr_t);

	GetCurTime(&query_sta);
#endif

	//Step 5
	//insert back the remaining size 2 pattern
	for (int i = 1; i < k; i++) {
		if (seedResult[i] != NULL) {
			fsi_t *fsi_v = seedResult[i]->head->next;
			while (fsi_v != NULL) {
				if (seedResult[i + 1] == NULL
						|| is_maximal(fsi_v, seedResult[i + 1])) {
					//					print_fsi(fsi_v, _fp);
					fsi_t *fsi_v2 = copy_fsi(fsi_v);
					fsi_v2->sup = fsi_v->sup;
					add_fsi_set_entry(maximalResult[i], fsi_v2);
				}
				fsi_v = fsi_v->next;
			}
		}
	}

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

//constrcut L2 by apriori
//i=1 = L2
//i=2 = L3
fsi_set_t* const_L_k_apriori(int alg_opt, int numOfObj, fsi_set_t *L, int i,
B_KEY_TYPE &sup_avg) {
// >=L_2
	fsi_t *fsi_v1, *fsi_v2, *fsi_v3;

	FEA_TYPE fea_v1;
	B_KEY_TYPE sup, sup_sum;

	sup_sum = 0;
	fsi_set_t *fsi_set_cur = alloc_fsi_set();

	if (L->fsi_n == 0)
		return fsi_set_cur;

	fsi_v1 = L->head->next;
	while (fsi_v1->next != NULL) {
		fsi_v2 = fsi_v1->next;
		while (fsi_v2 != NULL) {
			// join step.
			if ((fea_v1 = join_check(fsi_v1, fsi_v2)) != -1) {
				fsi_v3 = add_fsi(fsi_v1, fea_v1);
				// prune step.
				if (i > 1 && !all_subesets_exist(L, fsi_v3)) {
					/*s*/
					stat_v.memory_v -= sizeof(fsi_t)
							+ fsi_v3->fea_n * sizeof(FEA_TYPE);
					/*s*/

					release_fsi(fsi_v3);
					fsi_v2 = fsi_v2->next;
					continue;
				}
				// count step.
				sup = comp_support(alg_opt, fsi_v3, numOfObj, NULL,
				NULL);

				if (sup >= min_sup) {
					sup_sum += sup;
					add_fsi_set_entry(fsi_set_cur, fsi_v3);
				} else {
					/*s*/
					stat_v.memory_v -= sizeof(fsi_t)
							+ fsi_v3->fea_n * sizeof(FEA_TYPE);
					/*s*/
					release_fsi(fsi_v3);
				}
			}
			fsi_v2 = fsi_v2->next;
		}
		fsi_v1 = fsi_v1->next;
	}

	if (fsi_set_cur->fsi_n > 0)
		sup_avg = sup_sum / (double) fsi_set_cur->fsi_n;

	/*t*/
// print out each level for debug
	if (false) {
		printf("L_%d\n", 2);
		int cnt = 0;
		fsi_t *fsi_v = fsi_set_cur->head->next;
		while (fsi_v != NULL) {
			print_fsi(fsi_v, stdout);

			fsi_v = fsi_v->next;
			cnt++;
		}
	}

	return fsi_set_cur;
}

//convert to adjacency list

vector<list<int>> genGraph(int numOfFea, fsi_set_t *fsi_set_v) {

	vector<list<int>> adj;

	for (int i = 0; i < numOfFea; i++) {
		list<int> tmp;

		adj.push_back(tmp);
	}

//create edges for each pattern
	fsi_t *fsi_v = fsi_set_v->head->next;
	while (fsi_v != NULL) {

		//create edge for each pair of node
		for (int i = 0; i < fsi_v->fea_n; i++) {
			for (int j = i + 1; j < fsi_v->fea_n; j++) {
				FEA_TYPE fea1 = fsi_v->feaset[i];
				FEA_TYPE fea2 = fsi_v->feaset[j];
				if (std::find(adj[fea1 - 1].begin(), adj[fea1 - 1].end(),
						fea2 - 1) == adj[fea1 - 1].end()) {
					adj[fea1 - 1].push_back(fea2 - 1);
					adj[fea2 - 1].push_back(fea1 - 1);
				}
			}
		}

		fsi_v = fsi_v->next;
	}

	//sort list
	for (int i = 0; i < numOfFea; i++) {

		adj[i].sort();
	}
	return adj;
}

//convert cliques to pattern
//@maxSize = maximum clique size
fsi_set_t* candidate_maximal_patterns(set<vector<int>> cliques, int &maxSize) {

	fsi_set_t *fsi_set_v;

	fsi_set_v = alloc_fsi_set();

	for (auto c : cliques) {
		fsi_t *fsi_v = alloc_fsi(c.size());
		for (int i = 0; i < c.size(); i++) {
			fsi_v->feaset[i] = c[i] + 1; //offset fix!
		}
		add_fsi_set_entry(fsi_set_v, fsi_v);
		//---
		if (fsi_v->fea_n > maxSize)
			maxSize = fsi_v->fea_n;
	}

	return fsi_set_v;

}

//distribute the fsi_v in @fsi_set_v by its size
fsi_set_t** distribute_by_size(fsi_set_t *fsi_set_v, int maxSize) {
	fsi_set_t **cmp; //candidate maximal patterns by size
	cmp = (fsi_set_t**) (malloc(maxSize * sizeof(fsi_set_t*)));
	memset(cmp, 0, maxSize * sizeof(fsi_set_t*));
	/*s*/
	stat_v.memory_v += maxSize * sizeof(fsi_set_t*);
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;

	/*s*/
	for (int i = 0; i < maxSize; i++) {
		cmp[i] = alloc_fsi_set();
	}
	//	printf("maxSize:%d\n", maxSize);
	//distribute CMP to L
	fsi_t *fsi_v = fsi_set_v->head->next;
	while (fsi_v != NULL) {
		add_fsi_set_entry(cmp[fsi_v->fea_n - 1], copy_fsi(fsi_v));
		fsi_v = fsi_v->next;
	}
//		release_fsi_set(fsi_set_v);

	return cmp;
}
/**
 * @minSize = minimum size that need to be considered for maximal pattern
 * @maxSize = maximum size that ...
 *
 */
void find_maximal(int minSize, int maxSize, int numOfObj, fsi_set_t **cmp,
		fsi_set_t **result) {

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
//		printf("------ size:%d cmp[k]:%d --------\n", k, cmp[k - 1]->fsi_n);
		fsi_t *fsi_cur = cmp[k - 1]->head->next;
		while (fsi_cur != NULL) {
			if (k < maxSize && checked_subset_fsi(result[k], fsi_cur)) {
				fsi_cur = fsi_cur->next;
				continue;
			}
			double sup = comp_support(46, fsi_cur, numOfObj, fsi_cur->parentRI,
			NULL);

//			//===
//			print_fsi(fsi_cur, stdout);
//			exit(-1);
////			//===

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
				fsi_t *fsi_tmp = contain_fsi(cmp[k - 2], fsi_next);
				if (fsi_tmp != NULL) {
					//insert to its parentRI
					stat_v.memory_v -= fsi_tmp->parentRI->size() * sizeof(int);
					fsi_tmp->parentRI->insert(fsi_cur->RI->begin(),
							fsi_cur->RI->end());
					stat_v.memory_v += fsi_tmp->parentRI->size() * sizeof(int);

				} else if (!checked_subset_fsi(result, maxSize, fsi_next)) {
					fsi_t *fsi_new = copy_fsi(fsi_next);
					//---
					//				fsi_next->parentRI = fsi_cur->RI;
					fsi_new->parentRI = new unordered_set<int>();
					fsi_new->parentRI->insert(fsi_cur->RI->begin(),
							fsi_cur->RI->end());
					stat_v.memory_v += fsi_cur->RI->size() * sizeof(int);
					//---
					add_fsi_set_entry(cmp[k - 2], fsi_new);
				}

				//----
				//replace fsi_next[i] by fsi_cur[i] in each iteration
				//cur[i+1] is absent in each iteration in this way
				for (int i = 0; i < fsi_next->fea_n; i++) {
					//skipping i-th label to form k-1 subset
					fsi_next->feaset[i] = fsi_cur->feaset[i];
					//inserting to k-1 level/iteration
					fsi_t *fsi_tmp = contain_fsi(cmp[k - 2], fsi_next);
					if (fsi_tmp != NULL) {
						//insert to its parentRI
						stat_v.memory_v -= fsi_tmp->parentRI->size()
								* sizeof(int);
						fsi_tmp->parentRI->insert(fsi_cur->RI->begin(),
								fsi_cur->RI->end());
						stat_v.memory_v += fsi_tmp->parentRI->size()
								* sizeof(int);

					} else if (!checked_subset_fsi(result, maxSize, fsi_next)) {
						fsi_t *fsi_new = copy_fsi(fsi_next);
						//---
						//				fsi_next->parentRI = fsi_cur->RI;
						fsi_new->parentRI = new unordered_set<int>();
						fsi_new->parentRI->insert(fsi_cur->RI->begin(),
								fsi_cur->RI->end());
						stat_v.memory_v += fsi_cur->RI->size() * sizeof(int);
						//---
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
		//memory clear
		//releasing k+1 level RI
		//this (i.e., k-th) level RI need to be used as k-1 parentRI !
		//some fsi_v have NULL RI because skipped in above loop
		if (k < maxSize) {
			fsi_t *fsi_v1 = cmp[k]->head->next;
			while (fsi_v1 != NULL) {
				if (fsi_v1->RI != NULL) {
					/*s*/
					stat_v.memory_v -= fsi_v1->RI->size() * sizeof(int);
					/*s*/
					fsi_v1->RI->clear();
					fsi_v1->RI = NULL;
				}
				fsi_v1 = fsi_v1->next;
			}
			//---
			fsi_v1 = cmp[k - 1]->head->next;
			while (fsi_v1 != NULL) {
				if (fsi_v1->parentRI != NULL) {
					/*s*/
					stat_v.memory_v -= fsi_v1->parentRI->size() * sizeof(int);
					/*s*/
					fsi_v1->parentRI->clear();
					fsi_v1->parentRI = NULL;
				}
				fsi_v1 = fsi_v1->next;
			}
			//---
			release_fsi_set(cmp[k]);
		}
		k--;
	}		//end while
}

//return true if fsi_v is in fsi_set_cv
fsi_t* contain_fsi(fsi_set_t *fsi_set_v, fsi_t *fsi_v) {

	fsi_t *fsi_v1 = fsi_set_v->head->next;

	while (fsi_v1 != NULL) {
		if (fsi_v1->fea_n == fsi_v->fea_n) {
			if (memcmp(fsi_v1->feaset, fsi_v->feaset,
					(fsi_v1->fea_n) * sizeof(FEA_TYPE)) == 0) {
				return fsi_v1;
			}
		}
		//---
		fsi_v1 = fsi_v1->next;
	}

	return NULL;
}

// return true if fsi_v1 is subset of fsi_v2
bool is_subset_fsi(fsi_t *fsi_v2, fsi_t *fsi_v) {

	FEA_TYPE *fset1 = fsi_v->feaset;
	FEA_TYPE *fset2 = fsi_v2->feaset;

	int i = 0;		//for v1
	int j = 0;		//for v2

	//subset check that only work if both fset are sorted
	//in ascending order
	while (i < fsi_v->fea_n && j < fsi_v2->fea_n) {
		if (fset1[i] < fset2[j]) {
			return false;
		} else if (fset1[i] > fset2[j]) {
			j++;
		} else {
			i++;
			j++;
		}
	}
	if (i == fsi_v->fea_n)
		return true;
	else
		return false;
}

//return true if fsi_v is any subset of fsi_set_v
bool checked_subset_fsi(fsi_set_t *fsi_set_v, fsi_t *fsi_v) {

	fsi_t *fsi_v2 = fsi_set_v->head->next;

	while (fsi_v2 != NULL) {
//		if (fsi_v2->fea_n > fsi_v->fea_n) {
//			//checking each label in fsi_v
//			int cnt = 0;
//			int j = 0;
//			for (int i = 0; i < fsi_v2->fea_n; i++) {
//				if (fsi_v2->feaset[i] == fsi_v->feaset[j]) {
//					cnt++;
//					j++;
//					if (cnt == fsi_v->fea_n)	//all matches, it is a subset
//						return true;
//				}
//			}
//		}
		if (is_subset_fsi(fsi_v2, fsi_v))
			return true;
		//---
		fsi_v2 = fsi_v2->next;
	}

	return false;
}

//return true if fsi_v is any subset of fsi_set_v[k]
//where k={fsi_v->fea_n+1, maxSize}
bool checked_subset_fsi(fsi_set_t **fsi_set_v, int maxSize, fsi_t *fsi_v) {

	if (fsi_v->fea_n >= maxSize)
		return false;

	for (int i = fsi_v->fea_n + 2; i < maxSize; i++) {
		if (checked_subset_fsi(fsi_set_v[i - 1], fsi_v))
			return true;
	}

	return false;
}
