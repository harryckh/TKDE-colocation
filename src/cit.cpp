#include "cit.h"

extern colocation_stat_t stat_v;

cit_t* cit_ini() {
	cit_t *T;

	T = (cit_t*) malloc(sizeof(cit_t));
	memset(T, 0, sizeof(cit_t));

	T->root = alloc_cit_node();
	T->root->parent = NULL;

	T->root->obj_v = NULL;
	/*s*/
	stat_v.memory_v += sizeof(cit_t);
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return T;
}

void cit_release_sub(cit_node_t *x) {
	if (x->child != NULL) {
		cit_node_t *tmp = x->child;
		while (tmp != NULL) {
			cit_release_sub(tmp);
			tmp = tmp->next;
		}

	}

	free(x);

	/*s*/
	stat_v.memory_v -= sizeof(cit_node_t);
	/*s*/
}

/*
 *	Release the binary search tree T.
 */
void cit_release(cit_t *T) {
	if (T != NULL) {
		if (T->root != NULL)
			cit_release_sub(T->root);
		free(T);

		/*s*/
		stat_v.memory_v -= sizeof(cit_t);
		/*s*/
	}
}

cit_node_t* alloc_cit_node() {
	cit_node_t *cit_node_v;

	cit_node_v = (cit_node_t*) malloc(sizeof(cit_node_t));
	memset(cit_node_v, 0, sizeof(cit_node_t));

	//an empty head
	cit_node_v->child = (cit_node_t*) malloc(sizeof(cit_node_t));
	memset(cit_node_v->child, 0, sizeof(cit_node_t));

	cit_node_v->child_num = 0;

	/*s*/
	stat_v.memory_v += sizeof(cit_node_t) * 2;
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/
	return cit_node_v;
}

//insert a new child to x
cit_node_t* cit_insert(cit_t *T, cit_node_t *x, obj_t *obj_v) {

	cit_node_t *cit_node_v = alloc_cit_node();
	cit_node_v->obj_v = obj_v;
	cit_node_v->parent = x;

	//insert to x
	cit_node_v->next = x->child->next;
	x->child->next = cit_node_v;

	x->child_num++;

	T->node_n++;

	return cit_node_v;
}

void print_tree(cit_node_t *x) {
	if (x != NULL) {
		if (x->obj_v != NULL)
			printf("%d\t%d\n", x->obj_v->id, x->obj_v->fea);
		if (x->child_num > 0) {
			cit_node_t *cit_node_child = x->child->next;
			while (cit_node_child != NULL) {
				print_tree(cit_node_child);
				cit_node_child = cit_node_child->next;
			}
		}
	}
	return;
}

void print_ancestor_sub(cit_node_t *cit_node_v) {

	if (cit_node_v->parent != NULL) {
//		printf("%d(%d)\t", cit_node_v->obj_v->id, cit_node_v->obj_v->fea);
		printf("%d\t", cit_node_v->obj_v->fea);
		print_ancestor_sub(cit_node_v->parent);
	}
}
//print the (leaf) node and all its ancestors
void print_ancestor(cit_node_t *cit_node_v) {

	printf("print_ancestor: ");
	print_ancestor_sub(cit_node_v);
	printf("\n");
}


//return true if the (sub)-tree is correct
bool check_tree_sub(cit_node_t *cit_node_v, fsi_t *fsi_v, int depth) {

	if (depth > -1 && cit_node_v->obj_v->fea != fsi_v->feaset[depth])
		return false;

	cit_node_t *cit_node_child = cit_node_v->child->next;
	while (cit_node_child != NULL) {
		if (!check_tree_sub(cit_node_child, fsi_v, depth + 1)) {
			return false;
		}
		cit_node_child = cit_node_child->next;
	}
	return true;
}

//check if the tree is constructed as expected
bool check_tree(cit_t *T, fsi_t *fsi_v) {

	return check_tree_sub(T->root, fsi_v, -1);
}
