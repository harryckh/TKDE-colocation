

#ifndef DATA_STRUCT2_H
#define DATA_STRUCT2_H

#include "irtree.h"


//From cao_alg.h.
//The structure of the node for storing an irtree node entry.
typedef struct node_list
{
	node_t*				node_v;
	struct node_list*	next;
}	node_list_t;

//The structure of the list for storing a set of irtree node entries.
typedef struct cns
{
	node_list_t*	list_head;
	int				node_n;

	struct cns*		next;		//For cns_list_t.
}	cns_t;

//The structure of storing an cns entry.
typedef struct cns_list
{
	cns_t*				head;
	int					cns_n;

	struct cns_list*	next;	//For cns_list_set_t.
}	cns_list_t;

//The structure for storing a set of cns_list_t structures.
typedef struct cns_list_set
{
	cns_list_t*	head;
	int			list_n;
} cns_list_set_t;

//The structure for storing a set of object sets (obj_set_t).
typedef struct obj_set_list
{
	int						obj_set_n;
	obj_set_t*				obj_set_v;
	struct obj_set_list*	next;

	struct obj_set_list*	down;		//For cross list.
}	obj_set_list_t;
	

#endif