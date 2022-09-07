//
//  data_struct_operation.cpp
//  colocation
//
// ** for "basic" operations like create/release object_set_t / fsi_t
//
//  Created by Harry on 26/10/2018.
//  Copyright © 2018 Harry. All rights reserved.
//

#include <stdio.h>
#include "data_struct_operation.h"

/*
 *	Allocate a fsi_t structure.
 */
fsi_t* alloc_fsi(int key_n) {
	fsi_t *fsi_v;

	fsi_v = (fsi_t*) malloc(sizeof(fsi_t));
	memset(fsi_v, 0, sizeof(fsi_t));

	fsi_v->fea_n = key_n;

	fsi_v->feaset = (FEA_TYPE*) malloc(key_n * sizeof(FEA_TYPE));
	memset(fsi_v->feaset, 0, key_n * sizeof(FEA_TYPE));

	/*s*/
	stat_v.memory_v += sizeof(fsi_t) + key_n * sizeof(FEA_TYPE);
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return fsi_v;
}

/*
 *	Allocate a fsi_set_t structure.
 */
fsi_set_t* alloc_fsi_set() {
	fsi_set_t *fsi_set_v;

	fsi_set_v = (fsi_set_t*) malloc(sizeof(fsi_set_t));
	memset(fsi_set_v, 0, sizeof(fsi_set_t));

	fsi_set_v->head = (fsi_t*) malloc(sizeof(fsi_t));
	memset(fsi_set_v->head, 0, sizeof(fsi_t));

	/*s*/
	stat_v.memory_v += sizeof(fsi_set_t) + sizeof(fsi_t);
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return fsi_set_v;
}

/*
 *	Allocate a fsi_set_list_t structure.
 */
/*
 fsi_set_list_t* alloc_fsi_set_list(){
 
 fsi_set_list_t* fsi_set_list_v;
 
 fsi_set_list_v = ( fsi_set_list_t*)malloc( sizeof( fsi_set_list_t));
 memset( fsi_set_list_v, 0, sizeof( fsi_set_list_t));
 
 fsi_set_list_v->head = (fsi_set_t*) malloc( sizeof(fsi_set_t));
 memset(fsi_set_list_v->head, 0, sizeof(fsi_set_t));
 
 
 stat_v.memory_v += sizeof( fsi_set_t) + sizeof( fsi_t);
 if( stat_v.memory_v > stat_v.memory_max)
 stat_v.memory_max = stat_v.memory_v;
 
 
 return fsi_set_list_v;
 }
 */

/*
 *	Release a fsi_t structure.
 */
void release_fsi(fsi_t *fsi_v) {
	/*s*/
	stat_v.memory_v -= sizeof(fsi_t) + fsi_v->fea_n * sizeof(FEA_TYPE);
	/*s*/

	free(fsi_v->feaset);
	free(fsi_v);
}

/*
 *	Release the memory of an fsi_set_t structure.
 */
void release_fsi_set(fsi_set_t *fsi_set_v) {
	fsi_t *fsi_v1, *fsi_v2;

	if (!fsi_set_v)
		return;

	/*s*/
	stat_v.memory_v -= sizeof(fsi_set_t)
			+ (fsi_set_v->fsi_n + 1) * sizeof(fsi_t);
	/*s*/

	fsi_v1 = fsi_set_v->head;
	while (fsi_v1->next != NULL) {
		fsi_v2 = fsi_v1->next;
		release_fsi(fsi_v1);

		fsi_v1 = fsi_v2;
	}
	release_fsi(fsi_v1);

	free(fsi_set_v);
}

/*
 *	Release the memory of an fsi_set_list_t structure.
 */
/*
 void release_fsi_set_list( fsi_set_list_t* fsi_set_list_v)
 {
 fsi_set_t* fsi_set_v1, *fsi_set_v2;
 
 if( !fsi_set_list_v)
 return;
 
 
 stat_v.memory_v -= sizeof( fsi_set_list_t) + sizeof( fsi_set_t);
 
 fsi_set_v1 = fsi_set_list_v->head;
 while( fsi_set_v1->next != NULL)
 {
 fsi_set_v2 = fsi_set_v1->next;
 release_fsi_set(fsi_set_v1);
 
 fsi_set_v1 = fsi_set_v2;
 }
 release_fsi_set(fsi_set_v1);
 
 free( fsi_set_list_v);
 }
 */

void add_fsi_set_entry(fsi_set_t *fsi_set_v, fsi_t *fsi_v) {
	fsi_v->next = fsi_set_v->head->next;
	fsi_set_v->head->next = fsi_v;

	fsi_set_v->fsi_n++;
}

//remove the first item
void remove_fsi_set_entry(fsi_set_t *fsi_set_v) {
	if(fsi_set_v->fsi_n==0)
		return;

	fsi_t* fsi_v = fsi_set_v->head->next;
	fsi_set_v->head->next = fsi_v->next;

	release_fsi(fsi_v);

	fsi_set_v->fsi_n--;
}

/*
 *  returning a new fsi = @fsi_v + fea
 */
fsi_t* add_fsi(fsi_t *fsi_v, FEA_TYPE fea) {
	fsi_t *rtn;

	rtn = alloc_fsi(fsi_v->fea_n + 1);

	memcpy(rtn->feaset, fsi_v->feaset, fsi_v->fea_n * sizeof(FEA_TYPE));

	//--maintain ordering
	// note that rtn->feaset[rtn->fea_n-1] is the LAST element
	if (rtn->feaset[rtn->fea_n - 2] > fea) {
		rtn->feaset[rtn->fea_n - 1] = rtn->feaset[rtn->fea_n - 2];
		rtn->feaset[rtn->fea_n - 2] = fea;
	} else
		rtn->feaset[rtn->fea_n - 1] = fea;

	//--

	return rtn;
}

/*
 *  returning a new fsi = @fsi_v
 */
fsi_t* copy_fsi(fsi_t *fsi_v) {
	fsi_t *rtn;

	rtn = alloc_fsi(fsi_v->fea_n);

	memcpy(rtn->feaset, fsi_v->feaset, fsi_v->fea_n * sizeof(FEA_TYPE));

	rtn->parentRI= fsi_v->parentRI;

	return rtn;
}

/*
 * constrcut an inverted file using obj_v->fea for the objects in @obj_set_v
 *
 */
bst_t* const_IF(obj_set_t *obj_set_v) {
	bst_t *inverted_file;

	bst_node_t *bst_node_v;
	obj_node_t *obj_node_v;

	inverted_file = bst_ini();

	// Insert all the objects to construct the IF
	obj_node_v = obj_set_v->head->next;
	while (obj_node_v != NULL) {
		bst_node_v = bst_search(inverted_file, obj_node_v->obj_v->fea);

		if (bst_node_v != NULL) {
			// Update the posting list.
			add_obj_set_entry(obj_node_v->obj_v, bst_node_v->p_list_obj);
		} else // bst_node_v = NULL
		{
			// Insert a new keyword in the binary tree.
			bst_node_v = (bst_node_t*) malloc(sizeof(bst_node_t));
			memset(bst_node_v, 0, sizeof(bst_node_t));

			/*s*/
			stat_v.memory_v += sizeof(bst_node_t);
			if (stat_v.memory_v > stat_v.memory_max)
				stat_v.memory_max = stat_v.memory_v;
			/*s*/

			// Update the posting list.
			bst_node_v->key = obj_node_v->obj_v->fea;
			bst_node_v->p_list_obj = alloc_obj_set();

			add_obj_set_entry(obj_node_v->obj_v, bst_node_v->p_list_obj);
			bst_insert(inverted_file, bst_node_v);
		}
		obj_node_v = obj_node_v->next;
	}

	return inverted_file;
}

void release_IF_sub(bst_node_t *x) {
	release_obj_set(x->p_list_obj);

	if (x->left != NULL)
		release_IF_sub(x->left);
	if (x->right != NULL)
		release_IF_sub(x->right);

	free(x);

	/*s*/
	stat_v.memory_v -= sizeof(bst_node_t);
	/*s*/
}

/*
 *	Release the IF @T.
 *  Different from bst_release, here we need to release p_list_obj
 */
void release_IF(bst_t *T) {
	if (T != NULL) {
		if (T->root != NULL)
			release_IF_sub(T->root);
		free(T);

		/*s*/
		stat_v.memory_v -= sizeof(bst_t);
		/*s*/
	}
}

//-------------------------------------------------------------

/*
 *	Print the @fsi in @o_fp.
 */
void print_fsi(fsi_t *fsi_v, FILE *o_fp) {
	for (int i = 0; i < fsi_v->fea_n; i++) {
		fprintf(o_fp, "%d  ", fsi_v->feaset[i]);
	}

	fprintf(o_fp, "%lf  ", fsi_v->sup);

	fprintf(o_fp, "\n");
}

/*
 *	Print the @fsi in @o_fp.
 */
void print_fsi_set(fsi_set_t *fsi_set_v, FILE *o_fp) {
	fprintf(o_fp, "%d\n", fsi_set_v->fsi_n);

	fsi_t *fsi_v = fsi_set_v->head->next;
	while (fsi_v != NULL) {
		print_fsi(fsi_v, o_fp);
		fsi_v = fsi_v->next;

	}
}
void print_fsi_set(fsi_set_t **result, int numOfFea, FILE *o_fp) {
	for (int i = 0; i < numOfFea; i++) {
		if (result[i] != NULL && result[i]->fsi_n > 0) {
			fsi_t *fsi_v = result[i]->head->next;

			if (fsi_v != NULL)
				fprintf(o_fp, "%d\n", result[i]->fsi_n);

			while (fsi_v != NULL) {
				print_fsi(fsi_v, o_fp);

				fsi_v = fsi_v->next;
			}
			fprintf(o_fp, "\n");
		}
	}
}

void print_maximal_fsi_set(fsi_set_t **result, int numOfFea, FILE *o_fp) {
	for (int i = 0; i < numOfFea; i++) {
		if (result[i] != NULL) {
			fsi_t *fsi_v = result[i]->head->next;

			if (fsi_v != NULL)
				fprintf(o_fp, "%d\n", result[i]->fsi_n);

			while (fsi_v != NULL) {
				if (i == numOfFea - 1 || is_maximal(fsi_v, result[i + 1]))
					print_fsi(fsi_v, o_fp);

				fsi_v = fsi_v->next;
			}
			fprintf(o_fp, "\n");
		}
	}
}

/*
 * check whether any set in fsi_set_v is superset of fsi_cur
 */

bool is_maximal(fsi_t *fsi_cur, fsi_set_t *fsi_set_v) {
	fsi_t *fsi_v;

	fsi_v = fsi_set_v->head->next;
	while (fsi_v != NULL) {
		if (is_subset(fsi_cur, fsi_v))
			return false;
		fsi_v = fsi_v->next;
	}
	return true;
}

/**
 * Check whether fsi_v1 is a subset of fsi_v2
 */
bool is_subset(fsi_t *fsi_v1, fsi_t *fsi_v2) {
	if (fsi_v1->fea_n > fsi_v2->fea_n)
		return false;

	for (int i = 0; i < fsi_v1->fea_n; i++) {
		bool tag = false;
		FEA_TYPE fea = fsi_v1->feaset[i];
		for (int j = 0; j < fsi_v2->fea_n; j++) {
			if (fea == fsi_v2->feaset[j]) {
				tag = true;
				break;
			}
		}
		if (!tag)
			return false;
	}
	return true;
}

/*
 * Check whether @fea is in @fsi_v
 */
bool has_fea(fsi_t *fsi_v, FEA_TYPE fea) {
	for (int i = 0; i < fsi_v->fea_n; i++) {
		if (fsi_v->feaset[i] == fea)
			return true;
	}
	return false;
}

/*
 * check whether two object sets interest
 */
bool check_obj_set_equal_interest(obj_set_t *v1, obj_set_t *v2) {
	obj_node_t *obj_node_v1, *obj_node_v2;

	obj_node_v1 = v1->head->next;
	while (obj_node_v1 != NULL) {
		obj_node_v2 = v2->head->next;
		while (obj_node_v2 != NULL) {
			if (obj_node_v1->obj_v == obj_node_v2->obj_v) {
				return true;
			}
			obj_node_v2 = obj_node_v2->next;
		}

		obj_node_v1 = obj_node_v1->next;
	}
	return false;
}

/*
 *	Check whether the object set @O_t coveres the feasible set @fsi.
 *  ** possible running time improvement by iterating obj ?
 */
bool check_feasibility(obj_set_t *O_t, fsi_t *fsi_v) {
	int tag;
	obj_node_t *obj_node_v;

	for (int i = 0; i < fsi_v->fea_n; i++) {
		tag = 0;
		obj_node_v = O_t->head->next;
		while (obj_node_v != NULL) {
			if (obj_node_v->obj_v->fea == fsi_v->feaset[i]) {
				tag = 1;
				break;
			}

			obj_node_v = obj_node_v->next;
		}

		if (tag == 0)
			return false;
	}

	return true;
}

/*
 *  convert fsi to psi
 *  @feaToExclude will not be included in the resulting psi
 *
 */
psi_t* fsi_to_psi(fsi_t *fsi_v, FEA_TYPE feaToExclude) {
	psi_t *psi_v;
	psi_v = alloc_psi();

	for (int i = 0; i < fsi_v->fea_n; i++) {
		if (fsi_v->feaset[i] != feaToExclude)
			add_psi_entry(psi_v, fsi_v->feaset[i]);
	}

	return psi_v;
}

//-------------------------------------------------------------

/*
 *	Retrieve all the objects located at the sub-tree rooted at @node_v.
 *	The retrieved objects are stored in obj_set_v.
 */
void retrieve_sub_tree(node_t *node_v, obj_set_t *&obj_set_v) {
	int i;
	BIT_TYPE p_list;

	if (node_v->level == 0) {
		// node_v is a leaf-node.
		// Retrieve all its objects.
		for (i = 0; i < node_v->num; i++) {
			add_obj_set_entry((obj_t*) (node_v->child[i]), obj_set_v);
		}
	} else {
		// node_v is an inner-node.
		// Invoke the function recursively.
		for (i = 0; i < node_v->num; i++) {
			retrieve_sub_tree((node_t*) (node_v->child[i]), obj_set_v);
		}
	}
}

/*
 *	Range query on the sub-tree rooted at @node_v.
 *	@disk_v indicates the range which is a circle.
 *
 *	The results are stored in @obj_set_v.
 */
void range_query_sub(node_t *node_v, disk_t *disk_v, obj_set_t *&obj_set_v) {
	int i;
	BIT_TYPE p_list;
	range *MBR;

	if (node_v->parent == NULL)
		MBR = get_MBR_node(node_v, IRTree_v.dim);
	else
		MBR = node_v->parent->MBRs[node_v->loc];

	// No overlapping.
	if (!is_overlap(MBR, disk_v))
		return;

	// Enclosed entrely.
	if (is_enclosed(MBR, disk_v)) {
		retrieve_sub_tree(node_v, obj_set_v);
		if (node_v->parent == NULL) {
			free(MBR);

			/*s*/
			stat_v.memory_v -= IRTree_v.dim * sizeof(range);
			/*s*/
		}

		return;
	}

	// The remaining cases.
	if (node_v->level == 0) {
		// node_v is a leaf-node.
		for (i = 0; i < node_v->num; i++) {
			if (is_enclosed(((obj_t*) (node_v->child[i]))->MBR, disk_v))
				add_obj_set_entry((obj_t*) (node_v->child[i]), obj_set_v);
		}
	} else {
		// node_v is an inner-node.

		for (i = 0; i < node_v->num; i++) {
			range_query_sub((node_t*) (node_v->child[i]), disk_v, obj_set_v);
		}
	}

	if (node_v->parent == NULL) {
		free(MBR);
		/*s*/
		stat_v.memory_v -= IRTree_v.dim * sizeof(range);
		/*s*/
	}
}

/*
 *	Circle range query.
 *
 *	DFS: recursive implementation.
 */
obj_set_t* range_query(disk_t *disk_v) {
	obj_set_t *obj_set_v;

	obj_set_v = alloc_obj_set();
	range_query_sub(IRTree_v.root, disk_v, obj_set_v);

	return obj_set_v;
}

