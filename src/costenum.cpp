
/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */

#include "costenum.h"


/*
 *	Allocate a loc_t structure.
 */
loc_t* alloc_loc( int dim)
{
	loc_t* loc_v;

	loc_v = ( loc_t*)malloc( sizeof( loc_t));
	memset( loc_v, 0, sizeof( loc_t));

	loc_v->coord = ( float*)malloc( dim * sizeof( float));
	memset( loc_v->coord, 0, dim * sizeof( float));

	loc_v->dim = dim;

	return loc_v;
}

/*
 *	Retrieve the loc_t information of an obj_t structure @obj_v.
 */
loc_t* get_obj_loc( obj_t* obj_v)
{
	int i;
	loc_t* loc_v;

	loc_v = alloc_loc( IRTree_v.dim);
	for( i=0; i<IRTree_v.dim; i++)
		loc_v->coord[ i] = obj_v->MBR[ i].min;
	
	return loc_v;
}

/*
 *	Copy a loc_t structure.
 */
loc_t* copy_loc( loc_t* loc_v)
{
	int j;
	loc_t* loc_v1;

	loc_v1 = alloc_loc( loc_v->dim);

	for( j=0; j<loc_v->dim; j++)
		loc_v1->coord[ j] = loc_v->coord[ j];

	return loc_v1;
}

/*
 *	Release a loc_t structure.
 */
void release_loc( loc_t* loc_v)
{
	free( loc_v->coord);
	free( loc_v);
}

/*
 *	Allocate a psi_t structure.
 */
psi_t* alloc_psi( )
{
	psi_t* psi_v;

	psi_v = ( psi_t*)malloc( sizeof( psi_t));
	memset( psi_v, 0, sizeof( psi_t));

	psi_v->k_head = ( k_node_t*)malloc( sizeof( k_node_t));
	memset( psi_v->k_head, 0, sizeof( k_node_t));

    /*s*/
    stat_v.memory_v += sizeof(psi_t) + sizeof( k_node_t);
    if( stat_v.memory_v > stat_v.memory_max)
        stat_v.memory_max = stat_v.memory_v;
    /*s*/

	return psi_v;
}

/*
 *	Add an keyword into the psi_t structure.
 */
void add_psi_entry( psi_t* psi_v, KEY_TYPE key)
{
	k_node_t* k_node_v;

	k_node_v = ( k_node_t*)malloc( sizeof( k_node_t));
	memset( k_node_v, 0, sizeof( k_node_t));

    /*s*/
    stat_v.memory_v += sizeof( k_node_t);
    if( stat_v.memory_v > stat_v.memory_max)
        stat_v.memory_max = stat_v.memory_v;
    /*s*/

    
	k_node_v->key = key;
	k_node_v->next =  psi_v->k_head->next;
	psi_v->k_head->next = k_node_v;
	psi_v->key_n ++;
}

/*
 *	Release a psi_t structure.
 */
void release_psi( psi_t* psi_v)
{
	k_node_t* k_node_v1, *k_node_v2;

    if(psi_v==NULL)
        return;
    
    /*s*/
    stat_v.memory_v -= sizeof( k_node_t)* (psi_v->key_n+1);
    stat_v.memory_v -= sizeof( psi_t);
    /*s*/
    
	k_node_v1 = psi_v->k_head;
	while( k_node_v1->next != NULL)
	{
		k_node_v2 = k_node_v1->next;
		free( k_node_v1);
		k_node_v1 = k_node_v2;
	}
	free( k_node_v1);
    free( psi_v);
    
    
}

/*
 *	Allocate a query_t structure.
 */
query_t* alloc_query( )
{
	query_t* q;

	q = ( query_t*)malloc( sizeof( query_t));
	memset( q, 0, sizeof( query_t));

	//q->loc_v = alloc_loc( dim);

	//q->psi_v = alloc_psi( );

	return q;
}

/*
 *	Print a query_t structure.
 */
void print_query( query_t* q, FILE* o_fp)
{
	int i;
	k_node_t* k_node_v;

	//Location.
	fprintf( o_fp, "%f", q->loc_v->coord[ 0]);
	for( i=1; i<q->loc_v->dim; i++)
	{
		fprintf( o_fp, ",%f", q->loc_v->coord[ i]);
	}
	fprintf( o_fp, "\n");

	//Keywords.
	fprintf( o_fp, "%i", q->psi_v->key_n);
	k_node_v = q->psi_v->k_head->next;
	while( k_node_v != NULL)
	{
		fprintf( o_fp, ",%i", ( int)( k_node_v->key));
		k_node_v = k_node_v->next;
	}
	fprintf( o_fp, "\n\n");

	return;
}

/*
 *	Read a query instance from a file stream @i_fp.
 */
query_t* read_query( FILE* i_fp)
{
	int i;
	char des;
	KEY_TYPE key;
	query_t* q;
	k_node_t* k_node_v;

	q = alloc_query( );

	q->loc_v = alloc_loc( IRTree_v.dim);
	q->psi_v = alloc_psi( );

	//Location.
	if( fscanf( i_fp, "%f", &q->loc_v->coord[ 0]) == EOF)
	{
		release_query( q);
		return NULL;
	}

	for( i=1; i<IRTree_v.dim; i++)
		fscanf( i_fp, "%c%f", &des, &q->loc_v->coord[ i]);

	//Keywords.
	k_node_v = q->psi_v->k_head;
	fscanf( i_fp, "%i", &q->psi_v->key_n);
	for( i=0; i<q->psi_v->key_n; i++)
	{
		fscanf( i_fp, "%c%lf", &des, &key);
		add_keyword_entry( k_node_v, key);
	}

	return q;		
}

/*
 *	Release a query_t structure.
 */
void release_query( query_t* q)
{
    
    if(q==NULL) return;
	release_loc( q->loc_v);
	release_psi( q->psi_v);

	free( q);
}

/*
 *	Allocate the memory for a disk_t structure.
 */
disk_t* alloc_disk( int dim)
{
	disk_t* disk_v;

	disk_v = ( disk_t*)malloc( sizeof( disk_t));
	memset( disk_v, 0, sizeof( disk_t));

	disk_v->center = alloc_loc( dim);

	/*s*/
	stat_v.memory_v += sizeof( disk_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/
	

	return disk_v;
}

/*
 *	Initialize a disk_t structure.
 */
void set_disk( disk_t* disk_v, loc_t* loc_v, B_KEY_TYPE radius)
{
	int i;

	disk_v->radius = radius;
	
	for( i=0; i<loc_v->dim; i++)
		disk_v->center->coord[ i] = loc_v->coord[ i];
}

/*
 *	Construct a disk_t structure with its center of @loc_v and
 *	its radius of @radius.
 */
disk_t* const_disk( loc_t* loc_v, B_KEY_TYPE radius)
{
	int i;
	disk_t* disk_v;


	disk_v = alloc_disk( IRTree_v.dim);
	for( i=0; i<IRTree_v.dim; i++)
		disk_v->center->coord[ i] = loc_v->coord[ i];
	disk_v->center->dim = IRTree_v.dim;

	disk_v->radius = radius;

	return disk_v;
}

/*
 *	Release a disk_t structure.
 */
void release_disk( disk_t* disk_v)
{
    if(disk_v == NULL) return;
    
	release_loc( disk_v->center);
	free( disk_v);

	/*s*/
	stat_v.memory_v -= sizeof( disk_t);
	/*s*/
}

/*
 *	Allocate the memory for an obj_set_t structure.
 */
obj_set_t* alloc_obj_set( )
{
	obj_set_t* obj_set_v;

	obj_set_v = ( obj_set_t*)malloc( sizeof( obj_set_t));
	memset( obj_set_v, 0, sizeof( obj_set_t));

	obj_set_v->head = ( obj_node_t*)malloc( sizeof( obj_node_t));
	memset( obj_set_v->head, 0, sizeof( obj_node_t));

	/*s*/
	stat_v.memory_v += sizeof( obj_set_t) + sizeof( obj_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return obj_set_v;
}

/*
 *	Copy an obj_set_t structure.
 */
obj_set_t* copy_obj_set( obj_set_t* obj_set_v)
{
	obj_set_t* rtn;
	obj_node_t* obj_node_v, *obj_node_iter;

	rtn= alloc_obj_set( );
	obj_node_v = rtn->head;

	obj_node_iter = obj_set_v->head->next;
	while( obj_node_iter != NULL)
	{
		obj_node_v->next = ( obj_node_t*)malloc( sizeof( obj_node_t));
		memset( obj_node_v->next, 0, sizeof( obj_node_t));

		obj_node_v = obj_node_v->next;
		obj_node_v->obj_v = obj_node_iter->obj_v;
		
		obj_node_iter = obj_node_iter->next;
	}

	rtn->obj_n = obj_set_v->obj_n;

	/*s*/
	stat_v.memory_v += rtn->obj_n * sizeof( obj_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return rtn;		
}

/*
 *	Check whether two obj_set_t structures obj_set_v1 and obj_set_v2 have the same content.
 */
bool has_same_content_obj_set( obj_set_t* obj_set_v1, obj_set_t* obj_set_v2)
{
	obj_node_t* iter1, *iter2;

	if( obj_set_v1->obj_n != obj_set_v2->obj_n)
		return false;

	iter1 = obj_set_v1->head->next;
	while( iter1 != NULL)
	{
		iter2 = obj_set_v2->head->next;
		while( iter2 != NULL)
		{
			if( iter2->obj_v == iter1->obj_v)
				break;

			iter2 = iter2->next;
		}

		if( iter2 == NULL)
			return false;

		iter1 = iter1->next;
	}

	return true;
}

/*
 *	Remove the identical objs from @obj_set_v.
 */
void remove_identical_obj( obj_set_t* obj_set_v)
{
	obj_node_t* obj_node_iter1, *obj_node_iter2, *tmp;

	obj_node_iter1 = obj_set_v->head->next;
	while( obj_node_iter1->next != NULL)
	{
		tmp = obj_node_iter1;	
		obj_node_iter2 = obj_node_iter1->next;
		while( obj_node_iter2 != NULL)
		{
			if( obj_node_iter2->obj_v == obj_node_iter1->obj_v)
			{
				//remove.
				tmp->next = obj_node_iter2->next;
				obj_set_v->obj_n --;
				
				free( obj_node_iter2);
				
				/*s*/
				stat_v.memory_v -= sizeof( obj_node_t);
				/*s*/

				obj_node_iter2 = tmp->next;
				continue;
			}

			tmp = obj_node_iter2;
			obj_node_iter2 = obj_node_iter2->next;
		}

		obj_node_iter1 = obj_node_iter1->next;
		if( obj_node_iter1 == NULL)
			break;
	}

	return;
}

/*
 *	Print an obj_set_t structure.
 */
void print_obj_set( obj_set_t* obj_set_v, FILE* o_fp)
{
	int i;
	obj_node_t* obj_node_iter;
	obj_t* obj_v;

	if( !obj_set_v)
		return;

	for( i=0; i<20; i++)
		fprintf( o_fp, "==");
	fprintf( o_fp, "\n");

	fprintf( o_fp, "%i\n\n", obj_set_v->obj_n);
	obj_node_iter = obj_set_v->head->next;
	while( obj_node_iter != NULL)
	{
		obj_v = obj_node_iter->obj_v;
		fprintf( o_fp, "%i:\t", obj_v->id);
        
//        print_k_list( obj_v->k_head, o_fp);
        fprintf( o_fp, "%d \n ", obj_v->fea);

		obj_node_iter = obj_node_iter->next;
	}
	fprintf( o_fp, "\n");
}

/*
 *	Release the memory of an obj_set_t structure.
 */
void release_obj_set( obj_set_t* obj_set_v)
{
	obj_node_t* obj_node_v1, *obj_node_v2;

	if( obj_set_v == NULL)
		return;

	/*s*/
	stat_v.memory_v -= sizeof( obj_set_t) + 
						( obj_set_v->obj_n + 1) * sizeof( obj_node_t);
	/*s*/

	obj_node_v1 = obj_set_v->head;
	while( obj_node_v1->next != NULL)
	{
		obj_node_v2 = obj_node_v1->next;
		free( obj_node_v1);

		obj_node_v1 = obj_node_v2;
	}
	free( obj_node_v1);

	free( obj_set_v);
}

/*
 *	Allocate a IF_entry_t structure.
 *
IF_entry_t* alloc_IF_entry( )
{
	IF_entry_t* IF_entry_v;

	IF_entry_v = ( IF_entry_t*)malloc( sizeof( IF_entry_t));
	memset( IF_entry_v, 0, sizeof( IF_entry_t));

	IF_entry_v->p_list = alloc_obj_set( );

	return IF_entry_v;
}*/

/*
 *	Release a IF_entry_t structure.
 *
void release_IF_entry( IF_entry_t* IF_entry_v)
{
	release_obj_set( IF_entry_v->p_list);
	free( IF_entry_v);
}*/

/*
 *	Allocate a IF_t structure.
 *
IF_t* alloc_IF( )
{
	int i;
	IF_t* IF_v;

	IF_v = ( IF_t*)malloc( sizeof( IF_t));
	memset( IF_v, 0, sizeof( IF_t));

	IF_v->bst_v = bst_ini( );

	//IF_v->key_n = key_n;

	//IF_v->entry_v = ( IF_entry_t*)malloc( key_n * sizeof( IF_entry_t));
	//memset( IF_v->entry_v, 0, key_n * sizeof( IF_entry_t));

	//for( i=0; i<key_n; i++)
		//IF_v->entry_v[ i].p_list = alloc_obj_set( );
		//IF_v->entry_v[ i] = alloc_IF_entry( );

	return IF_v;
}*/

/*
 *	Release a IF_t structure.
 *
void release_IF( IF_t* IF_v)
{
	//int i;
	
	//for( i=0; i<IF_v->key_n; i++)
		//release_obj_set( IF_v->entry_v[ i].p_list);
		//release_IF_entry( IF_v->entry_v[ i]);

	//free( IF_v->entry_v);
	bst_release( IF_v->bst_v);
	free( IF_v);
}*/

/*
 *	Calculate the distance between two locations @loc_v1 and @loc_v2.
 */
B_KEY_TYPE calc_dist_loc( loc_t* loc_v1, loc_t* loc_v2)
{
	int i;
	B_KEY_TYPE dist;

	dist = 0;
	for( i=0; i<loc_v1->dim; i++)
		dist += pow( loc_v1->coord[ i] - loc_v2->coord[ i], 2);

	return sqrt( dist);
}

/*
 *	Calculate the distance between two objects @obj_v1 and @obj_v2.
 */
B_KEY_TYPE calc_dist_obj( obj_t* obj_v1, obj_t* obj_v2)
{
	B_KEY_TYPE dist;
	loc_t* loc_v1, *loc_v2;

	loc_v1 = get_obj_loc( obj_v1);
	loc_v2 = get_obj_loc( obj_v2);

	dist = calc_dist_loc( loc_v1, loc_v2);

	release_loc( loc_v1);
	release_loc( loc_v2);

	return dist;	
}

/*
 *	Check whether an object @obj_v contains the keyword @key.
 */
bool has_key_obj( obj_t* obj_v, KEY_TYPE key)
{
//	k_node_t* k_node_v;

//	k_node_v = obj_v->k_head->next;
//	while( k_node_v != NULL)
//	{
//		if( k_node_v->key == key)
//			return true;
//
//		k_node_v = k_node_v->next;
//	}

    if(obj_v->fea==key)return true;
    
	return false;
}

/*
 *	Check whether an object @obj_v is "relevant" to the query @q.
 *	That is, whether @obj_v contains a keyword in the query @q.
 */
bool is_relevant_obj( obj_t* obj_v, query_t* q)
{
	KEY_TYPE key;
	k_node_t* k_node_v;

	k_node_v = q->psi_v->k_head->next;
	while( k_node_v != NULL)
	{
		key = k_node_v->key;
		if( has_key_obj( obj_v, key))
			return true;

		k_node_v = k_node_v->next;
	}

	return false;		
}

/*
 *	Check whether the keywords in the query @q are covered by a set of objs in @obj_set_v.
 */
bool is_covered_obj_set( obj_set_t* obj_set_v, query_t* q)
{
	KEY_TYPE key;
	k_node_t* k_node_iter;
	obj_node_t* obj_node_iter;

	k_node_iter = q->psi_v->k_head->next;
	while( k_node_iter != NULL)
	{
		key = k_node_iter->key;

		obj_node_iter = obj_set_v->head->next;
		while( obj_node_iter != NULL)
		{
			if( has_key_obj( obj_node_iter->obj_v, key))
				break;

			obj_node_iter = obj_node_iter->next;
		}

		if( obj_node_iter == NULL)
			return false;

		k_node_iter = k_node_iter->next;
	}
	
	return true;
}

/*
 *	Check whether the sub-tree rooted at a node @node_v
 *	contains a keyword @key.
 *
 */
BIT_TYPE has_key_node( node_t* node_v, KEY_TYPE key)
{
	bst_node_t* bst_node_v;

	if( ( bst_node_v = bst_search( node_v->bst_v, key)))
		return bst_node_v->p_list;
	
	return 0;
}

/*
 *	Check whether a node @node_v is "relevant" or not.
 */
BIT_TYPE is_relevant_node( node_t* node_v, query_t* q)
{
	BIT_TYPE res, res_t;

	k_node_t* k_node_v;

	res = 0;
	k_node_v = q->psi_v->k_head->next;
	while( k_node_v != NULL)
	{
		res_t = has_key_node( node_v, k_node_v->key);
		union_bit( res, res_t);

		if( res == UINT_MAX)
			return res;

		k_node_v = k_node_v->next;
	}

	return res;
}

/*
 *	Calculate the minimum distance between the MBR of a node @node_v 
 *	and a location @loc_v.
 */
B_KEY_TYPE calc_minDist( range* MBR, loc_t* loc_v)
{
	int i;
	B_KEY_TYPE m_dist;
	
	//Calculate the minimum distance between the MBR and the location.
	m_dist = 0;
	for( i=0; i<loc_v->dim; i++)
	{
		if( loc_v->coord[ i] < MBR[ i].min)
			m_dist += pow( MBR[ i].min - loc_v->coord[ i], 2);
		else if( loc_v->coord[ i] > MBR[ i].max)
			m_dist += pow( loc_v->coord[ i] - MBR[ i].max, 2);
	}

	return sqrt( m_dist);
}

/*
 *	Calculate the maxDist between a MBR @MBR and 
 */
B_KEY_TYPE calc_maxDist( range* MBR, loc_t* loc_v)
{
	int i;
	B_KEY_TYPE maxDist;

	maxDist = 0;
	for( i=0; i<loc_v->dim; i++)
	{
		if( loc_v->coord[ i] <= ( MBR[ i].min + MBR[ i].max) / 2)
			maxDist += pow( MBR[ i].max - loc_v->coord[ i], 2);
		else
			maxDist += pow( loc_v->coord[ i] - MBR[ i].min, 2);
	}

	return sqrt( maxDist);
}

/*
 *	Calculate the minimum distance between the MBR of a node @node_v and a location @loc_v.
 */
B_KEY_TYPE calc_minDist_node( node_t* node_v, loc_t* loc_v)
{
	B_KEY_TYPE dist;
	range* MBR;

	if( node_v->parent == NULL)
		MBR = get_MBR_node( node_v, IRTree_v.dim);
	else
		MBR = node_v->parent->MBRs[ node_v->loc];

	dist = calc_minDist( MBR, loc_v);

	if( node_v->parent == NULL)
	{
		free( MBR);

		/*s*/
		stat_v.memory_v -= sizeof( IRTree_v.dim * sizeof( range));
		/*s*/
	}

	return dist;
}

/*
 *	Check whether a MBR @MBR overlaps with a disk @disk_v.
 */
bool is_overlap( range* MBR, disk_t* disk_v)
{
	B_KEY_TYPE min_dist;

	min_dist = calc_minDist( MBR, disk_v->center);
	if(  min_dist <= disk_v->radius)
		return true;
	else
		return false;
}

/*
 *	Constrained NN search.
 *	Find the NN that contains the keyword @key and
 *	is located in the dist @disk_v.
 *
 *	disk_v == NULL indicates that the disk is the whole space,
 *	thus, in this case, this constraint is invalid (interface consideration).
 *
 *	Method: best-first search.
 *	
 */
obj_t* const_NN_key( loc_t* loc_v, KEY_TYPE key, disk_t* disk_v)
{	
	int size, top, i, rear;
	BIT_TYPE p_list;

	B_KEY_TYPE min_dist, c_dist;
	range* MBR;
	obj_t* obj_v;
	node_t* node_v;
	b_heap_t* b_h;

	//Checking 1.
	if( IRTree_v.root->num == 0)
		return NULL;

	//Checking 2: involved in the search process implicitly.
	//if( !( bst_node_v = bst_search( IRTree.root->bst_v, key)))
		//return NULL;

	//Checking 3:
	MBR = get_MBR_node( IRTree_v.root, IRTree_v.dim);
	if( disk_v != NULL && !is_overlap( MBR, disk_v))
		return NULL;

	size = IRTree_v.obj_n / M + 1;
	b_h = alloc_b_heap( size);

	rear = 1;
	b_h->obj_arr[ rear].node_v = IRTree_v.root;
	b_h->obj_arr[ rear].key = calc_minDist( MBR, loc_v);
	free( MBR);

	/*s*/
	stat_v.memory_v -= sizeof( IRTree_v.dim * sizeof( range));
	/*s*/

	b_h_insert( b_h, rear++);

	min_dist = DBL_MAX;
	obj_v = NULL;

	while( !b_h_is_empty( b_h))
	{
		top = b_h_get_top( b_h);
		if( b_h->obj_arr[ top].key >= min_dist)
			break;

		node_v = b_h->obj_arr[ top].node_v;

		//Keyword constraint (for all its entries).
		//if( !( bst_node_v = bst_search( node_v->bst_v, key)))
			//continue;
		if( ( p_list = has_key_node( node_v, key)) == 0)
			continue;
		
		//bst_node_v != NULL.
		for( i=0; i<node_v->num; i++)
		{
			//Keyword constraint (for a specific entry).
			if( !get_k_bit( p_list, i))
				continue;

			//the i'th entry contains the keyword.
			if( node_v->level == 0)
			{
				//node_v is a leaf-node.
				MBR = ( ( obj_t*)( node_v->child[ i]))->MBR;
				
				//Disk constraint (if any, i.e., disk_v != NULL).
				if( disk_v != NULL && !is_overlap( MBR, disk_v))
					continue;

				c_dist = calc_minDist( MBR, loc_v);
				if( c_dist < min_dist)
				{
					//Update the current best solution.
					min_dist = c_dist;
					obj_v = ( obj_t*)( node_v->child[ i]);
				}
			}
			else
			{
				//node_v is an inner-node.
				MBR = node_v->MBRs[ ( ( node_t*)( node_v->child[ i]))->loc];

				//Disk constraint (if any).
				if( disk_v != NULL && !is_overlap( MBR, disk_v))
					continue;

				c_dist = calc_minDist( MBR, loc_v);
				if( c_dist < min_dist)
				{
					//Enqueue the child.
					b_h->obj_arr[ rear].node_v = ( node_t*)( node_v->child[ i]);
					b_h->obj_arr[ rear].key = c_dist;

					b_h_insert( b_h, rear++);
				}					
			}
		}//for
	}//while

	release_b_heap( b_h);
		
	return obj_v;
}

/*
 *	The kNN version of "const_NN_key".
 *
 *	Simply use a list for storing the kNN results.
 */
obj_set_t* const_k_NN_key( loc_t* loc_v, KEY_TYPE key, disk_t* disk_v, int k)
{	
	int size, top, i, rear, cnt;
	BIT_TYPE p_list;
	B_KEY_TYPE min_dist, c_dist;
	range* MBR;
	//obj_t* obj_v;
	node_t* node_v;
	b_heap_t* b_h;
	obj_node_t* obj_node_v, *obj_node_v1, *obj_node_rear;
	obj_set_t* obj_set_v;

	//Checking 1.
	if( IRTree_v.root->num == 0)
		return NULL;

	//Checking 2: involved in the search process implicitly.
	//if( !( bst_node_v = bst_search( IRTree.root->bst_v, key)))
		//return NULL;

	//Checking 3:
	MBR = get_MBR_node( IRTree_v.root, IRTree_v.dim);
	if( disk_v != NULL && !is_overlap( MBR, disk_v))
		return NULL;

	//Initialize a heap.
	size = IRTree_v.obj_n / M + 1;
	b_h = alloc_b_heap( size);

	rear = 1;
	b_h->obj_arr[ rear].node_v = IRTree_v.root;
	b_h->obj_arr[ rear].key = calc_minDist( MBR, loc_v);
	free( MBR);

	/*s*/
	stat_v.memory_v -= sizeof( IRTree_v.dim * sizeof( range));
	/*s*/

	b_h_insert( b_h, rear++);


	obj_set_v = alloc_obj_set( );
	min_dist = DBL_MAX;
	//obj_v = NULL;

	while( !b_h_is_empty( b_h))
	{
		top = b_h_get_top( b_h);
		if( b_h->obj_arr[ top].key >= min_dist)
			break;

		node_v = b_h->obj_arr[ top].node_v;

		//Keyword constraint (for all its entries).
		//if( !( bst_node_v = bst_search( node_v->bst_v, key)))
			//continue;
		if( ( p_list = has_key_node( node_v, key)) == 0)
			continue;
		
		//bst_node_v != NULL.
		for( i=0; i<node_v->num; i++)
		{
			//Keyword constraint (for a specific entry).
			if( !get_k_bit( p_list, i))
				continue;

			//the i'th entry contains the keyword.
			if( node_v->level == 0)
			{
				//node_v is a leaf-node.
				MBR = ( ( obj_t*)( node_v->child[ i]))->MBR;
				
				//Disk constraint (if any, i.e., disk_v != NULL).
				if( disk_v != NULL && !is_overlap( MBR, disk_v))
					continue;

				c_dist = calc_minDist( MBR, loc_v);
				if( c_dist < min_dist)
				{
					//Update the current best solution.
					/*
					min_dist = c_dist;
					obj_v = ( obj_t*)( node_v->child[ i]);
					*/
					obj_node_v = ( obj_node_t*)malloc( sizeof( obj_node_t));
					memset( obj_node_v, 0, sizeof( obj_node_t));

					obj_node_v->obj_v = ( obj_t*)( node_v->child[ i]);
					obj_node_v->dist = c_dist;

					obj_node_v1 = obj_set_v->head;
					while( obj_node_v1->next != NULL && obj_node_v1->next->dist < c_dist)
						obj_node_v1 = obj_node_v1->next;

					obj_node_v->next = obj_node_v1->next;
					obj_node_v1->next = obj_node_v;
					obj_set_v->obj_n ++;

					//Update min_dist.
					cnt = 0;
					obj_node_rear = obj_set_v->head;
					while( obj_node_rear->next != NULL && cnt < k)
					{
						obj_node_rear = obj_node_rear->next;

						cnt++;
					}
					min_dist = obj_node_rear->dist;

					if( obj_node_rear->next != NULL)
					{
						free( obj_node_rear->next);
						obj_node_rear->next = NULL;
						obj_set_v->obj_n --;
					}
				}
			}
			else
			{
				//node_v is an inner-node.
				MBR = node_v->MBRs[ ( ( node_t*)( node_v->child[ i]))->loc];

				//Disk constraint (if any).
				if( disk_v != NULL && !is_overlap( MBR, disk_v))
					continue;

				c_dist = calc_minDist( MBR, loc_v);
				if( c_dist < min_dist)
				{
					//Enqueue the child.
					b_h->obj_arr[ rear].node_v = ( node_t*)( node_v->child[ i]);
					b_h->obj_arr[ rear].key = c_dist;

					b_h_insert( b_h, rear++);
				}					
			}
		}//for
	}//while

	release_b_heap( b_h);
		
	return obj_set_v;
}

/*
 *	Check whether a MBR @MBR is enclosed entirely by a disk @disk_v.
 */
bool is_enclosed( range* MBR, disk_t* disk_v)
{
   // printf("a.%f\tb.%f\n", calc_maxDist( MBR, disk_v->center), disk_v->radius);
	if( calc_maxDist( MBR, disk_v->center) <= disk_v->radius)
		return true;
	else
		return false;
}

/*
 *	Add an object entry @obj_v to @obj_set_v.
 */
void add_obj_set_entry( obj_t* obj_v, obj_set_t* obj_set_v)
{
	obj_node_t* obj_node_v;

	obj_node_v = ( obj_node_t*)malloc( sizeof( obj_node_t));
	memset( obj_node_v, 0, sizeof( obj_node_t));

	obj_node_v->obj_v = obj_v;
	obj_node_v->next =  obj_set_v->head->next;
	obj_set_v->head->next = obj_node_v;
	obj_set_v->obj_n ++;

	/*s*/
	stat_v.memory_v += sizeof( obj_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/
}

/*
 *	Remove the first entry from the list of objects @obj_set_v.
 */	
void remove_obj_set_entry( obj_set_t* obj_set_v)
{
	obj_node_t* obj_node_v;

	obj_node_v = obj_set_v->head->next;
	obj_set_v->head->next = obj_node_v->next;
	obj_set_v->obj_n --;

	free( obj_node_v);

	/*s*/
	stat_v.memory_v -= sizeof( obj_node_t);
	/*s*/
}

/*
 *	Retrieve all the objects located at the sub-tree rooted at @node_v.
 *	The retrieved objects are stored in obj_set_v.
 */
void retrieve_sub_tree( node_t* node_v, obj_set_t* &obj_set_v, query_t* q)
{
	int i;
	BIT_TYPE p_list;

	if( node_v->level == 0)
	{
		//node_v is a leaf-node.
		//Retrieve all its objects.
		for( i=0; i<node_v->num; i++)
		{
			if( is_relevant_obj( ( obj_t*)( node_v->child[ i]), q))
				add_obj_set_entry( ( obj_t*)( node_v->child[ i]), obj_set_v);			
		}
	}
	else
	{
		//node_v is an inner-node.
		//Invoke the function recursively.
		p_list = is_relevant_node( node_v, q);
		for( i=0; i<node_v->num; i++)
		{
			if( get_k_bit( p_list, i))
				retrieve_sub_tree( ( node_t*)( node_v->child[ i]), obj_set_v, q);
		}
	}
}

/*
 *	Range query on the sub-tree rooted at @node_v.
 *	@disk_v indicates the range which is a circle.
 *
 *	The results are stored in @obj_set_v.
 */
void range_query_sub( node_t* node_v, disk_t* disk_v, obj_set_t* &obj_set_v, query_t* q)
{
	int i;
	BIT_TYPE p_list;
	range* MBR;

	if( node_v->parent == NULL)
		MBR = get_MBR_node( node_v, IRTree_v.dim);
	else
		MBR = node_v->parent->MBRs[ node_v->loc];

	//No overlapping.
	if( !is_overlap( MBR, disk_v))
		return;

	//Enclosed entrely.
	if( is_enclosed( MBR, disk_v))
	{
		retrieve_sub_tree( node_v, obj_set_v, q);
		if( node_v->parent == NULL)
		{
			free( MBR);

			/*s*/
			stat_v.memory_v -= IRTree_v.dim * sizeof( range);
			/*s*/
		}

		return;
	}

	//The remaining cases.
	if( node_v->level == 0)
	{
		//node_v is a leaf-node.
		for( i=0; i<node_v->num; i++)
		{
			if( is_enclosed( ( ( obj_t*)( node_v->child[ i]))->MBR, disk_v) &&
				is_relevant_obj( ( obj_t*)( node_v->child[ i]), q))
				add_obj_set_entry( ( obj_t*)( node_v->child[ i]), obj_set_v);
		}
	}
	else
	{
		//node_v is an inner-node.
		p_list = is_relevant_node( node_v, q);
		
		for( i=0; i<node_v->num; i++)
		{
			if( get_k_bit( p_list, i))
				range_query_sub( ( node_t*)( node_v->child[ i]), disk_v, obj_set_v, q);
		}
	}

	if( node_v->parent == NULL)
	{
		free( MBR);
		/*s*/
		stat_v.memory_v -= IRTree_v.dim * sizeof( range);
		/*s*/
	}
}

/*
 *	Circle range query.
 *
 *	DFS: recursive implementation.
 */
obj_set_t* range_query( disk_t* disk_v, query_t* q)
{
	obj_set_t* obj_set_v;

	obj_set_v = alloc_obj_set( );
	range_query_sub( IRTree_v.root, disk_v, obj_set_v, q);

	return obj_set_v;
}

/*
 *	Exclusion the objects on the boundary of disk @disk_v from @obj_set_v.
 */
void refine_region( obj_set_t* obj_set_v, disk_t* disk_v)
{
	B_KEY_TYPE dist;
	loc_t* loc_v;
	obj_node_t* obj_node_v1, *obj_node_v2;

	obj_node_v1 = obj_set_v->head;
	obj_node_v2 = obj_node_v1->next;
	while( obj_node_v2 != NULL)
	{
		loc_v = get_obj_loc( obj_node_v2->obj_v);
		dist = calc_dist_loc( loc_v, disk_v->center);
		
		release_loc( loc_v);
		
		if( dist == disk_v->radius)
		{
			obj_node_v1->next = obj_node_v2->next;
			free( obj_node_v2);
			obj_node_v2 = obj_node_v1->next;

			obj_set_v->obj_n --;

			/*s*/
			stat_v.memory_v -= sizeof( obj_node_t);
			/*s*/
			
			continue;
		}

		obj_node_v1 = obj_node_v2;
		obj_node_v2 = obj_node_v1->next;
	}
}

/*
 *	construct a psi_t structure based on an object @obj_v.
 */
psi_t* get_psi_obj( obj_t* obj_v)
{
	psi_t* psi_v;
	k_node_t* k_node_v, *k_node_v1;

	psi_v = alloc_psi( );
	k_node_v = psi_v->k_head;

//	k_node_v1 = obj_v->k_head->next;
//	while( k_node_v1 != NULL)
//	{	
//		add_keyword_entry( k_node_v, k_node_v1->key);
//
//		k_node_v1 = k_node_v1->next;
//	}

    add_keyword_entry(k_node_v, obj_v->fea);
    
	return psi_v;
}

/*
 *	Construct a psi_t structure based on a list of keywords @k_head.
 */
psi_t* const_psi( k_node_t* k_head)
{
	psi_t* psi_v;
	k_node_t* k_node_v;

	psi_v = ( psi_t*)malloc( sizeof( psi_t));
	memset( psi_v, 0, sizeof( psi_t));

	/*s*/
	stat_v.memory_v += sizeof( psi_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	psi_v->k_head = k_head;
	
	k_node_v = k_head->next;
	while( k_node_v != NULL)
	{
		psi_v->key_n ++;
		
		k_node_v = k_node_v->next;
	}

	return psi_v;
}

/*
 *	Exclude the keywords that occur in @k_head2 from @k_head1.
 *
 *	Return the resulting keywords in @k_head1.
 */
k_node_t* key_exclusion( k_node_t* k_head1, k_node_t* k_head2)
{
	int tag;
	k_node_t* k_node_v, *k_node_v1, *k_node_v2, *k_node_v3;

	k_node_v = ( k_node_t*)malloc( sizeof( k_node_t));
	memset( k_node_v, 0, sizeof( k_node_t));
	k_node_v3 = k_node_v;

	/*s*/
	stat_v.memory_v += sizeof( k_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	k_node_v1 = k_head1->next;
	while( k_node_v1 != NULL)
	{
		tag = 0;
		k_node_v2 = k_head2->next;
		while( k_node_v2 != NULL)
		{
			if( k_node_v2->key == k_node_v1->key)
			{
				tag = 1;
				break;
			}

			k_node_v2 = k_node_v2->next;
		}

		if( tag == 0)
		{
			//The current keyword should not excluded.
			add_keyword_entry( k_node_v3, k_node_v1->key);
		}

		k_node_v1 = k_node_v1->next;
	}

	return k_node_v;
}


/*
 *	Exclude the keywords @key from @k_head1.
 *
 *	Return the resulting keywords in @k_head1.
 */
k_node_t* key_exclusion( k_node_t* k_head1, KEY_TYPE key)
{
    int tag;
    k_node_t* k_node_v, *k_node_v1, *k_node_v2, *k_node_v3;
    
    k_node_v = ( k_node_t*)malloc( sizeof( k_node_t));
    memset( k_node_v, 0, sizeof( k_node_t));
    k_node_v3 = k_node_v;
    
    /*s*/
    stat_v.memory_v += sizeof( k_node_t);
    if( stat_v.memory_v > stat_v.memory_max)
        stat_v.memory_max = stat_v.memory_v;
    /*s*/
    
    k_node_v1 = k_head1->next;
    while( k_node_v1 != NULL)
    {
        tag = 0;
//        k_node_v2 = k_head2->next;
//        while( k_node_v2 != NULL)
//        {
//            if( k_node_v2->key == k_node_v1->key)
//            {
//                tag = 1;
//                break;
//            }
//            
//            k_node_v2 = k_node_v2->next;
//        }
        if( k_node_v1->key != key)
//        if( tag == 0)
        {
            //The current keyword should not excluded.
            add_keyword_entry( k_node_v3, k_node_v1->key);
        }
        
        k_node_v1 = k_node_v1->next;
    }
    
    return k_node_v;
}

/*
 *	Exclude the keywords that occur in @psi_v2 from @psi_v1.
 *
 *	Return the resulting keywords in @psi_v1.
 */
psi_t* psi_exclusion( psi_t* psi_v1, psi_t* psi_v2)
{
	int tag;
	k_node_t* k_node_v, *k_node_v1, *k_node_v2;
	psi_t* psi_v;

	psi_v = alloc_psi( );

	k_node_v = psi_v->k_head;
	k_node_v1 = psi_v1->k_head->next;
	while( k_node_v1 != NULL)
	{
		tag = 0;
		k_node_v2 = psi_v2->k_head->next;
		while( k_node_v2 != NULL)
		{
			if( k_node_v2->key == k_node_v1->key)
			{
				tag = 1;
				break;
			}

			k_node_v2 = k_node_v2->next;
		}

		if( tag == 0)
		{
			//The current keyword should not excluded.
			add_keyword_entry( k_node_v, k_node_v1->key);
			psi_v->key_n ++;
		}

		k_node_v1 = k_node_v1->next;
	}

	return psi_v;
}

/*
 *	Exclude the keywords @fae from @psi_v1.
 *
 *
 */
void psi_exclusion( psi_t* psi_v, FEA_TYPE fea)
{
 
    k_node_t* k_node_v1, *k_node_v2;
    
    k_node_v1 = psi_v->k_head;
    while( k_node_v1->next != NULL)
    {
        k_node_v2 = k_node_v1->next;
        
        if( k_node_v2->key == fea)
        {
            //The current keyword in *k_node_v2* should be excluded.
            k_node_v1->next = k_node_v2->next;
            free(k_node_v2);
            
            /*s*/
            stat_v.memory_v -= sizeof( k_node_t);
            /*s*/

            
            psi_v->key_n--;
        }
        else
            k_node_v1 = k_node_v1->next;
    }
    
    return;
}

/*
 *	Check whether the object set @O_t coveres the keyword set @psi.
 */
bool check_feasibility( obj_set_t* O_t, psi_t* psi)
{
	int tag;
	KEY_TYPE key;
	k_node_t* k_node_v, *k_node_v1;
	obj_node_t* obj_node_v;

    
    if(O_t->obj_n < psi->key_n)
        return false;
    
	k_node_v = psi->k_head->next;
	while( k_node_v != NULL)
	{
		key = k_node_v->key;

		tag = 0;
		obj_node_v =  O_t->head->next;
		while( obj_node_v != NULL)
		{
//			k_node_v1 = obj_node_v->obj_v->k_head->next;
//			while( k_node_v1 != 0)
//			{
//				if( k_node_v1->key == key)
//				{
//					tag = 1;
//					break;
//				}
//
//				k_node_v1 = k_node_v1->next;
//			}
//            
//			if( tag == 1)
//				break;

            if(obj_node_v->obj_v->fea == key)
            {
                tag = 1;
                break;
            }
            
			obj_node_v = obj_node_v->next;
		}

		if( tag == 0)
			return false;
		
		k_node_v = k_node_v->next;
	}

	return true;
}

/*
 *	Filter out the objects that are not located in range @disk_v from @O_t.
 */
void obj_filter_range( obj_set_t* &O_t, disk_t* disk_v)
{
	B_KEY_TYPE dist;
	loc_t* loc_v;
	obj_node_t* obj_node_v1, *obj_node_v2;

	obj_node_v1 = O_t->head;
	obj_node_v2 = obj_node_v1->next;
	while( obj_node_v2 != NULL)
	{
		loc_v = get_obj_loc( obj_node_v2->obj_v);
		dist = calc_dist_loc( loc_v, disk_v->center);
		if( dist > disk_v->radius)
		{
			//the obj is not located in the disk.
			obj_node_v1->next = obj_node_v2->next;
			free( obj_node_v2);
			obj_node_v2 = obj_node_v1->next;

			O_t->obj_n --;

			/*s*/
			stat_v.memory_v -= sizeof( obj_node_t);
			/*s*/
		}
		else
		{
			obj_node_v1 = obj_node_v2;
			obj_node_v2 = obj_node_v2->next;
		}

		release_loc( loc_v);
	}		
}

/*
 *	Construct the IF on a set of objects @obj_set_v for the keywords in @psi_v.
 *
 *	1. The IF structure is indexed by a binary search tree.
 *	2. No ordering is imposed in IF.
 */
bst_t* const_IF( obj_set_t* obj_set_v, psi_t* psi_v)
{
	int i;
	bst_t* IF_v;
	k_node_t* k_node_v;
	obj_node_t* obj_node_v;
	bst_node_t* bst_node_v;

	//IF_v = alloc_IF( psi_v->key_n);
	IF_v = bst_ini( );

	k_node_v = psi_v->k_head->next;
	//for( i=0; i<psi_v->key_n; i++)
    while(k_node_v !=NULL)
	{
		//IF_v->entry_v[ i].key = k_node_v->key;
		bst_node_v = ( bst_node_t*)malloc( sizeof( bst_node_t));
		memset( bst_node_v, 0, sizeof( bst_node_t));

		/*s*/
		stat_v.memory_v += sizeof( bst_node_t);
		/*s*/

		bst_node_v->key = k_node_v->key;
		bst_node_v->p_list_obj = alloc_obj_set( );
		
		obj_node_v = obj_set_v->head->next;
		while( obj_node_v != NULL)
		{
			if( has_key_obj( obj_node_v->obj_v, k_node_v->key))
				add_obj_set_entry( obj_node_v->obj_v, bst_node_v->p_list_obj);

			obj_node_v = obj_node_v->next;
		}

		bst_insert( IF_v, bst_node_v);

		k_node_v = k_node_v->next;
	}
	
	/*s*/
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return IF_v;
}

/*
 *	Construct a obj_set_t structure containing 3 objects
 *	from a tri_t structure @triplet_v.
 */
obj_set_t* const_obj_set( tri_t* triplet_v)
{
	obj_set_t* obj_set_v;

	obj_set_v = alloc_obj_set( );
	obj_set_v->obj_n = 3;

	//Including the objects.
	add_obj_set_entry( triplet_v->o, obj_set_v);
	add_obj_set_entry( triplet_v->o_1, obj_set_v);
	add_obj_set_entry( triplet_v->o_2, obj_set_v);

	return obj_set_v;
}

/*
 *	Check the distance constraint.
 */
bool check_dist_constraint( obj_set_t* obj_set_v, obj_t* obj_v, obj_t* o, B_KEY_TYPE d)
{
	obj_node_t* obj_node_iter;

	if( calc_dist_obj( o, obj_v) > d)
		return false;

	obj_node_iter = obj_set_v->head->next;
	while( obj_node_iter != NULL)
	{
		if( calc_dist_obj( obj_node_iter->obj_v, obj_v) > d)
			return false;

		obj_node_iter = obj_node_iter->next;
	}

	return true;
}

/*
 *	Update the IF structure @IF_v by removing the keywords 
 *	that have been covered by an object @obj_v.
 */
bst_node_list_t* update_IF_obj( bst_t* IF_v, obj_t* obj_v)
{
	k_node_t* k_node_v;
	bst_node_t* bst_node_v;
	bst_node_list_t* bst_node_list_v, *tmp;

	bst_node_list_v = ( bst_node_list_t*)malloc( sizeof( bst_node_list_t));
	memset( bst_node_list_v, 0, sizeof( bst_node_list_t));

	/*s*/
	stat_v.memory_v += sizeof( bst_node_list_t);
	/*s*/

//	k_node_v = obj_v->k_head->next;
//	while( k_node_v != NULL)
//	{
//      bst_node_v = bst_search( IF_v, k_node_v->key);

		bst_node_v = bst_search( IF_v, obj_v->fea);

		if( bst_node_v != NULL)
		{
		

		bst_delete( IF_v, bst_node_v);
		//bug.
		bst_node_v->p = NULL;
		bst_node_v->left = NULL;
		bst_node_v->right = NULL;
		
		tmp = ( bst_node_list_t*)malloc( sizeof( bst_node_list_t));
		memset( tmp, 0, sizeof( bst_node_list_t));

		/*s*/
		stat_v.memory_v += sizeof( k_node_t);
		/*s*/

		tmp->bst_node_v = bst_node_v;
		tmp->next = bst_node_list_v->next;
		bst_node_list_v->next = tmp;
        }
//		k_node_v = k_node_v->next;
//	}

	/*s*/
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return bst_node_list_v;
}

/*
 *	Release a bst_node_list_t structure.
 */
void release_bst_node_list( bst_node_list_t* bst_node_list_v)
{
	bst_node_list_t* tmp;
	
	while( bst_node_list_v != NULL)
	{
		tmp = bst_node_list_v->next;
		free( bst_node_list_v);
		bst_node_list_v = tmp;

		/*s*/
		stat_v.memory_v -= sizeof( bst_node_list_t);
		/*s*/
	}
}

/*
 *	Restore the IF_v structure @IF_v by re-including the bst_nodes of @bst_node_list_v.
 */
void restore_IF_bst_node_list( bst_t* IF_v, bst_node_list_t* bst_node_list_v)
{
	bst_node_list_t* bst_node_list_iter;

	bst_node_list_iter = bst_node_list_v->next;
	while( bst_node_list_iter != NULL)
	{
		bst_insert( IF_v, bst_node_list_iter->bst_node_v);

		bst_node_list_iter = bst_node_list_iter->next;
	}
}

/*
 *	Combine two obj_set_t structures @obj_set_v1 and @obj_set_v2.
 *		1. The combined result is stored in @obj_set_v1.
 *		2. The obj_set_v2 is released.
 *
 */
void combine_obj_set( obj_set_t* obj_set_v1, obj_set_t* obj_set_v2)
{
	obj_node_t* obj_node_v;

	//Locate the rear of obj_set_v1.
	obj_node_v = obj_set_v1->head;
	while( obj_node_v->next != NULL)
		obj_node_v = obj_node_v->next;

	obj_node_v->next = obj_set_v2->head->next;
	obj_set_v2->head->next = NULL;
	obj_set_v1->obj_n += obj_set_v2->obj_n;

	release_obj_set( obj_set_v2);
}

/*
 *	The sub-procedure of function "const_feasible_set_sub".
 */
obj_set_t* const_feasible_set_sub( bst_t* IF_v, obj_set_t* S_0, obj_t* o, B_KEY_TYPE d)
{
	obj_t* obj_v;
	obj_set_t* S;
	bst_node_t* bst_node_v;
	obj_node_t* obj_node_v;
	bst_node_list_t* bst_node_list_v;

	if( IF_v->node_n == 0)
		return alloc_obj_set( );			//An empty one.

	bst_node_v = IF_v->root;
	obj_node_v = bst_node_v->p_list_obj->head->next;
	while( obj_node_v != NULL)
	{
		//Pick an object.
		obj_v = obj_node_v->obj_v;

		//Distance constraint checking.
		if( !check_dist_constraint( S_0, obj_v, o, d))
		{
			obj_node_v = obj_node_v->next;
			continue;
		}
		
		//Update the IF_v.
		bst_node_list_v = update_IF_obj( IF_v, obj_v);

		//Update the S_0.
		//obj_v is added at the first place of S_0.
		add_obj_set_entry( obj_v, S_0);

		//Sub-procedure.
		S = const_feasible_set_sub( IF_v, S_0, o, d);

		//Restore the S_0.
		remove_obj_set_entry( S_0);

		//Restore the IF_v.
		restore_IF_bst_node_list( IF_v, bst_node_list_v);

		release_bst_node_list( bst_node_list_v);
		
		//Checking.
		if( S != NULL)
		{	
			//Include obj_v into S.
			add_obj_set_entry( obj_v, S);

			return S;
		}

		//S == NULL.

		//Try the next object candidate.
		obj_node_v = obj_node_v->next;
	}

	return NULL;
}

/*
 *	Check whether there exists a subset of @O_t such that
 *	1. it covers @psi and 2. its diameter is d.
 *
 *	Method: recursive implementation.
 *
 *	If so, return such a subset; otherwise, return NULL.
 */
obj_set_t* const_feasible_set( obj_set_t* O_t, psi_t* psi, obj_t* o, B_KEY_TYPE d)
{
	bst_t* IF_v;
	obj_set_t* S_0, *S;

	//Construct an IF structure.
	IF_v = const_IF( O_t, psi);

	//Initialize the S_0.
	S_0 = alloc_obj_set( );

	//Invoke the sub-procedure "recursively".
	S = const_feasible_set_sub( IF_v, S_0, o, d);

	//Release the resources.
	//bst_release( IF_v);
    release_IF(IF_v);
	release_obj_set( S_0);

	return S;
}

/*
 *	The implementation of the "ConstructFeasibleSet-Appro" procedure in the paper.
 */
obj_set_t* ConstructFeasibleSet_Appro( obj_t* o, query_t* q)
{
	B_KEY_TYPE radius;

	obj_t* obj_v;
	obj_set_t* S;
	k_node_t* k_head1, *k_node_v;
	loc_t* loc_v, *loc_v1;
	disk_t* disk_v;
	
	S = alloc_obj_set( );

	//Location of object o.
	loc_v = get_obj_loc( o);

	//Disk(q, d(o, q)).
	loc_v1 = copy_loc( q->loc_v);
	radius = calc_dist_loc( loc_v, q->loc_v);
	disk_v = const_disk( loc_v1, radius);

	//Include object o in S.
	add_obj_set_entry( o, S);

	//Obtain the "un-covered" keywords by S.
	k_head1 = key_exclusion( q->psi_v->k_head, o->fea);

	//Retrieve the disk-constrained NNs that contain the un-covered keywords.
	k_node_v = k_head1->next;
	while( k_node_v != NULL)
	{
		obj_v = const_NN_key( loc_v, k_node_v->key, disk_v);
		
		if( obj_v == NULL)
		{
			release_obj_set( S);
			S = NULL;
			break;
		}

		add_obj_set_entry( obj_v, S);

		k_node_v = k_node_v->next;
	}

	//Release the memory.
	release_loc( loc_v);
    release_loc( loc_v1);
	release_disk( disk_v);
	release_k_list( k_head1);
	
	return S;
}

/*
 *	Compute the farthest distance between a set of objects @obj_set_v and @q.
 */
B_KEY_TYPE comp_farthest( obj_set_t* obj_set_v, query_t* q)
{
	B_KEY_TYPE far_dist, dist;
	loc_t* loc_v;
	obj_node_t* obj_node_v;

	far_dist = 0;
	obj_node_v = obj_set_v->head->next;
	while( obj_node_v != NULL)
	{
		loc_v = get_obj_loc( obj_node_v->obj_v);
		
		dist = calc_dist_loc( loc_v, q->loc_v);
		if( dist > far_dist)
			far_dist = dist;

		release_loc( loc_v);
		obj_node_v = obj_node_v->next;
	}

	return far_dist;
}

/*
 *	Compute the diameter of a set of objects @obj_set_v.
 */
B_KEY_TYPE comp_diameter( obj_set_t* obj_set_v)
{
	B_KEY_TYPE dia, dist;
	obj_node_t* obj_node_v1, *obj_node_v2;

	if( obj_set_v->obj_n <= 1)
		return 0;

	//There exist at least two objects.
	dia = 0;
	obj_node_v1 = obj_set_v->head->next;
	while( obj_node_v1->next != NULL)
	{
		obj_node_v2 = obj_node_v1->next;
		while( obj_node_v2 != NULL)
		{			
			dist = calc_dist_obj( obj_node_v1->obj_v, obj_node_v2->obj_v);
			if( dist > dia)
				dia = dist;

			obj_node_v2 = obj_node_v2->next;
		}

		obj_node_v1 = obj_node_v1->next;
	}

	return dia;		
}


/*
 *	Compute the cost of set @obj_set_v wrt @q.
 */
B_KEY_TYPE comp_cost(int cost_tag, obj_set_t* obj_set_v, query_t* q)
{
	if( cost_tag == 1) 
		return comp_farthest( obj_set_v, q) + comp_diameter( obj_set_v);
	else
		return max( comp_farthest( obj_set_v, q), comp_diameter( obj_set_v));
}

/*
 *	Compute the lower bound @LB and the upper bound @UB of cost(S*, q).c_1,
 *	where S* is the optimal solution.
 *
 *	return NN(q, q.\psi) if successful; otherwise, return NULL.
 */
obj_set_t* comp_bounds( query_t* q, B_KEY_TYPE &LB, B_KEY_TYPE &UB, disk_t* disk_v)
{
	B_KEY_TYPE dist;
	k_node_t* k_node_v;
	obj_set_t* obj_set_v;
	loc_t* loc_v;
	obj_t* obj_v;

	obj_set_v = alloc_obj_set( );

	//Compute LB.
	LB = 0;
	k_node_v = q->psi_v->k_head->next;
	while( k_node_v != NULL)
	{
		obj_v = const_NN_key( q->loc_v, k_node_v->key, disk_v);
		
		if( obj_v == NULL)
		{
			release_obj_set( obj_set_v);
			return NULL;
		}

		loc_v = get_obj_loc( obj_v);
		dist = calc_dist_loc( loc_v, q->loc_v);
		if( dist > LB)
			LB = dist;
		
        release_loc(loc_v);
		add_obj_set_entry( obj_v, obj_set_v);

		k_node_v = k_node_v->next;
	}

	//Compute the UB.
	UB = comp_cost( cost_tag, obj_set_v, q);

	return obj_set_v;
	//release_obj_set( obj_set_v);
}

/*
 *	Exclude from a set of objects @obj_set_v
 *	those objects that are located in a specified disk @disk_v.
 *
 *	Note that those objects on the boundary of the disk are kept.
 */
void obj_exclusion_disk( obj_set_t* obj_set_v, disk_t* disk_v)
{
	B_KEY_TYPE far_dist;
	obj_node_t* obj_node_v1, *obj_node_v2;

	obj_node_v1 = obj_set_v->head;
	obj_node_v2 = obj_node_v1->next;
	while( obj_node_v2 != NULL)
	{
		far_dist = calc_maxDist( obj_node_v2->obj_v->MBR, disk_v->center);
		if( far_dist < disk_v->radius)
		{
			//The object should be excluded.
			obj_node_v1->next = obj_node_v2->next;
			free( obj_node_v2);
			obj_node_v2 = obj_node_v1->next;

			obj_set_v->obj_n --;

			/*s*/
			stat_v.memory_v -= sizeof( obj_node_t);
			/*s*/

			continue;
		}

		obj_node_v1 = obj_node_v2;
		obj_node_v2 = obj_node_v2->next;
	}
}


/*
 *	Exclude from a set of objects @obj_set_v
 *	those objects that are located *outside* a specified disk @disk_v.
 *
 *	Note that those objects on the boundary of the disk are kept.
 */
void obj_inclusion_disk( obj_set_t* obj_set_v, disk_t* disk_v)
{
    B_KEY_TYPE far_dist;
    obj_node_t* obj_node_v1, *obj_node_v2;
    
    obj_node_v1 = obj_set_v->head;
    obj_node_v2 = obj_node_v1->next;
    while( obj_node_v2 != NULL)
    {
        far_dist = calc_maxDist( obj_node_v2->obj_v->MBR, disk_v->center);
        if( far_dist > disk_v->radius)
        {
            //The object should be excluded.
            obj_node_v1->next = obj_node_v2->next;
            free( obj_node_v2);
            obj_node_v2 = obj_node_v1->next;
            
            obj_set_v->obj_n --;
            
            /*s*/
            stat_v.memory_v -= sizeof( obj_node_t);
            /*s*/
            
            continue;
        }
        
        obj_node_v1 = obj_node_v2;
        obj_node_v2 = obj_node_v2->next;
    }
}

/*
 *	Sort the objects in @obj_set_v by their distances to @q.
 *
 *	Method: heap-sort.
 *	Alternative method: based on the binary search tree.
 */
b_heap_t* heap_sort_obj_set( obj_set_t* obj_set_v, query_t* q)
{
	int cur;
	B_KEY_TYPE dist;
	b_heap_t* b_h;
	obj_node_t* obj_node_v;
	loc_t* loc_v;

	b_h = alloc_b_heap( obj_set_v->obj_n + 1);
	
	cur = 1;
	obj_node_v = obj_set_v->head->next;
	while( obj_node_v != NULL)
	{
		loc_v = get_obj_loc( obj_node_v->obj_v);
		dist = calc_dist_loc( loc_v, q->loc_v);
		release_loc( loc_v);

		b_h->obj_arr[ cur].key = dist;
		b_h->obj_arr[ cur].obj_v = obj_node_v->obj_v;

		b_h_insert( b_h, cur);
		
		cur ++;

		obj_node_v = obj_node_v->next;
	}

	return b_h;
}

/*
 *	The implementation of the "AchievabilityCheck" function in the paper.
 */
obj_set_t* AchievabilityCheck( tri_t triplet_v, query_t* q)
{
	//Disks.
	B_KEY_TYPE radius_1, radius_2;
	loc_t* loc_v1, *loc_v2, *loc_v3;
	disk_t* disk_v1, *disk_v2, *disk_v3, *disk_v_tmp;
	obj_set_t* O_t, *S;
	k_node_t* k_head_1, *k_head_2, *k_head_3;
	psi_t* psi_v;

	//Pre-checking.
	//Exclude the keywords covered by the triplet.
    
	k_head_1 = key_exclusion( q->psi_v->k_head, triplet_v.o->fea);
	k_head_2 = key_exclusion( k_head_1, triplet_v.o_1->fea);
	k_head_3 = key_exclusion( k_head_2, triplet_v.o_2->fea);

	psi_v = const_psi( k_head_3);

	release_k_list( k_head_1);
	release_k_list( k_head_2);
	//Note that k_head_3 would be released when psi_v is released.

	if( psi_v->key_n == 0)
	{
		// <o, o_1, o_2> is achievable.
		//S = {o, o_1, o_2} is a verying feasible set.
		S = const_obj_set( &triplet_v);

		release_psi( psi_v);
		return S;		
	}

	//Normal-checking.
	//Construct the disks.
	//1. Disk(q, d(o, q)).
	loc_v1 = get_obj_loc( triplet_v.o);
	radius_1 = calc_dist_loc( loc_v1, q->loc_v);
	disk_v1 = const_disk( q->loc_v, radius_1);

	//2. Disk(o_1, d(o_1, o_2)) and Disk(o_2, d(o_1, o_2)).
	loc_v2 = get_obj_loc( triplet_v.o_1);
	loc_v3 = get_obj_loc( triplet_v.o_2);
	radius_2 = calc_dist_loc( loc_v2, loc_v3);
	
	disk_v2 = const_disk( loc_v2, radius_2);
	disk_v3 = const_disk( loc_v3, radius_2);

	//Decide the disk for range query.
	//Use disk_v1 at default.
	if( radius_2 < radius_1)
	{
		//Change to use disk_v2 or disk_v3.
		disk_v_tmp = disk_v1;
		disk_v1 = disk_v2;
		disk_v2 = disk_v_tmp;
	}

	//Range query on disk_v1.
	O_t = range_query( disk_v1, q);
	
	//Filter the objs that are not located in disk_v2 or disk_v3.
	obj_filter_range( O_t, disk_v2);
	obj_filter_range( O_t, disk_v3);

	release_loc( loc_v1);
	release_loc( loc_v2);
	release_loc( loc_v3);
	release_disk( disk_v1);
	release_disk( disk_v2);
	release_disk( disk_v3);
	
	//Check whether O_t covers the keywords in the query.
	if( !check_feasibility( O_t, q->psi_v))
	{
		S = NULL;
		goto E;
	}

	/*s*/
	stat_v.O_t_size_sum += O_t->obj_n;
	/*s*/

	//Find a sub-set of O_t which covers psi_v.
	S = const_feasible_set( O_t, psi_v, triplet_v.o, radius_2);

	if( S == NULL)
		goto E;

	//Combine the sub-set with the critical objects.
	add_obj_set_entry( triplet_v.o, S);
	add_obj_set_entry( triplet_v.o_1, S);
	add_obj_set_entry( triplet_v.o_2, S);

	/*s*/
	stat_v.psi_n_sum += psi_v->key_n;
	/*s*/

E:	
	//Release the resources.
	release_psi( psi_v);
	release_obj_set( O_t);

	return S;
}

/*
 *	Handle one specific candidate of <o_1, o_2>.
 */
void process_obj_pair( obj_t* o_1, obj_t* o_2, B_KEY_TYPE dist, bst_t* obj_pair)
{
	bst_node_t* bst_node_v;	

	//Add the pair into the bst.
	bst_node_v = ( bst_node_t*)malloc( sizeof( bst_node_t));
	memset( bst_node_v, 0, sizeof( bst_node_t));

	/*s*/
	stat_v.memory_v += sizeof( bst_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	bst_node_v->key = dist;
	bst_node_v->obj_v1 = o_1;
	bst_node_v->obj_v2 = o_2;

	bst_insert( obj_pair, bst_node_v);
}

/*
 *	Initialize the set of <o_1, o_2> pairs.
 */
bst_t* ini_obj_pair( obj_set_t* region, obj_t* o, B_KEY_TYPE d_l, B_KEY_TYPE d_u, query_t* q)
{
	B_KEY_TYPE d_low, d_high, dist, d_o_q;
	obj_t* o_1, *o_2;
	bst_t* obj_pair;
	obj_node_t* obj_node_v1, *obj_node_v2;
	loc_t* loc_v, *loc_v1, *loc_v2;

	obj_pair = bst_ini( );

	if( region->obj_n == 0)
		return obj_pair;

	//Compute the upper bound of d(o_1, o_2): independent on <o_1, o_2>.
	loc_v = get_obj_loc( o);
	d_o_q = calc_dist_loc( loc_v, q->loc_v);

	if( cost_tag == 1)
		d_high = d_u - d_o_q;
	else
		d_high = d_u;
	release_loc( loc_v);

	obj_node_v1 = region->head->next;
	while( obj_node_v1 != NULL)
	{
		o_1 = obj_node_v1->obj_v;

		obj_node_v2 = obj_node_v1;
		while( obj_node_v2 != NULL)
		{
			o_2 = obj_node_v2->obj_v;

			//Compute the lower bound of d(o_1, o_2): dependent on <o_1, o_2>.
			if( cost_tag == 1)
			{
				loc_v1 = get_obj_loc( o_1);
				loc_v2 = get_obj_loc( o_2);
				d_low = d_l - min( calc_dist_loc( loc_v1, q->loc_v),
					calc_dist_loc( loc_v2, q->loc_v));	
				
				release_loc( loc_v1);
				release_loc( loc_v2);
			}
			else
			{
				d_low = d_o_q;
			}

			//process.
			dist = calc_dist_obj( o_1, o_2);
			if( cost_tag == 1)
			{
				if( !( dist >= d_high || dist < d_low))
					process_obj_pair( o_1, o_2, dist, obj_pair);
			}
			else
			{
				if( !( dist >= d_high || dist <= d_low))
					process_obj_pair( o_1, o_2, dist, obj_pair);
			}

			obj_node_v2 = obj_node_v2->next;
		}

		obj_node_v1 = obj_node_v1->next;
	}

	return obj_pair;
}

/*
 *	Update the set of <o_1, o_2> pairs.
 *		1. region should also be updated.
 */
void update_obj_pair( bst_t* obj_pair, obj_set_t* region, obj_t* o, B_KEY_TYPE d_l, B_KEY_TYPE d_u, query_t* q)
{
	B_KEY_TYPE dist, d_low, d_high, d_o_q;
	obj_node_t* obj_node_v;
	loc_t* loc_v;
	
	//Compute the upper bound of d(o_1, o_2): independent on <o_1, o_2>.
	loc_v = get_obj_loc( o);
	d_o_q = calc_dist_loc( loc_v, q->loc_v);
	if( cost_tag == 1)
		d_high = d_u - d_o_q;
	else
		d_high = d_u;
	release_loc( loc_v);

	//Remove the unsatisfied obj_pairs in the original obj_pair_v.
	if( d_high <= obj_pair->max)
		bst_trim( obj_pair, d_high);	

	//Include o into region. bug.
	add_obj_set_entry( o, region);
	
	obj_node_v = region->head->next;
	while( obj_node_v != NULL)
	{
		//Compute the lower bound of d(o_1, o_2): dependent on <o_1, o_2>.
		//The lower bound is equal to 0.

		if( cost_tag == 1)
			d_low = 0;
		else
			d_low = d_o_q;
		
		//process.
		dist = calc_dist_obj( obj_node_v->obj_v, o);
		if( cost_tag == 1)
		{
			if( !( dist >= d_high || dist < d_low))
				process_obj_pair( obj_node_v->obj_v, o, dist, obj_pair);
		}
		else
		{
			if( !( dist >= d_high || dist <= d_low))
				process_obj_pair( obj_node_v->obj_v, o, dist, obj_pair);
		}

		obj_node_v = obj_node_v->next;
	}
}

/*
 *	The implementation of the "ConstructFeasibleSet-Exact" procedure in the paper CoSK.
 */
obj_set_t* ConstructFeasibleSet_Exact( obj_t* o, query_t* q, bst_t* obj_pair)
{
	int j;
	obj_set_t* S;
	tri_t triplet_v;
	bst_node_t* cur_pair;
	obj_t* o_q;

	if( cost_tag == 2)
	{
		//Pre-checking for the Cost2 measurement.
		//Create o_q.
		o_q = ( obj_t*)malloc( sizeof( obj_t));
		alloc_obj( o_q, q->loc_v->dim);
		for( j=0; j<q->loc_v->dim; j++)
		{
			o_q->MBR[ j].min = q->loc_v->coord[ j];
			o_q->MBR[ j].max = o_q->MBR[ j].min;
		}
        //bug fix.
        o_q->fea = -1;
        //--
        
		S = NULL;
		triplet_v.o = o;
		triplet_v.o_1 = o;
		triplet_v.o_2 = o_q;
		S =  AchievabilityCheck( triplet_v, q);


		/*s*/
		stat_v.achi_sum ++;
		/*s*/

        if( S != NULL){

            return S;
        }
	}

	//Upper bound pruning: involved in CostEnum.

	S = NULL;
	triplet_v.o = o;
	cur_pair = bst_get_min( obj_pair->root);
	while( cur_pair != NULL)
	{
		triplet_v.o_1 = cur_pair->obj_v1;
		triplet_v.o_2 = cur_pair->obj_v2;

		//Lower bound pruning:
		if( calc_dist_obj( triplet_v.o_1, triplet_v.o_2) < 
			max( calc_dist_obj( triplet_v.o, triplet_v.o_1),
			calc_dist_obj( triplet_v.o, triplet_v.o_2)))
		{
			cur_pair = bst_successor( cur_pair);
			continue;
		}

		S = AchievabilityCheck( triplet_v, q);

		/*s*/
		stat_v.achi_sum ++;
		/*s*/
		
		if( S != NULL)
        {
			return S;
        }
		//S == NULL.
		cur_pair = bst_successor( cur_pair);
	}
    
	return S;
}

/*
 *	The implementation of "CostEnum".
 *
 *	s_tag = 1: CostEnum-Exact.
 *	s_tag = 2: CostEnum-Appro.
 *
 *  @disk_v: the disk D(o,d)
 *  @obj_set_v: the object in the disk @disk_v OR NULL
 */
obj_set_t* CostEnum( query_t* q, int s_tag, int prune_tag, obj_set_t* obj_set_v, disk_t* disk_v)
{
	int top;
	B_KEY_TYPE LB, UB, cost_c, cost, dist;
	obj_set_t* S_a;				//the current best solution.
	obj_set_t* S;				//the newly constructed feasible set.
	obj_set_t* R;				//region R.
	disk_t* disk_u, *disk_l;	//the outer and inner disks.
	obj_set_t* region_u, *region_l;
	obj_t* o, *o_next;
	loc_t* loc_v;
	bst_t* obj_pair = NULL;
	b_heap_t* R_heap = NULL;
    
	//Compute the LB and UB of cost(S*, q).c_1;
	S_a = comp_bounds( q, LB, UB, disk_v);
	
	if( S_a == NULL)
		return NULL;

    cost_c = UB;
    if( cost_c <= dist_thr)
    {
        return S_a;
    }
    

    //Initialize region R.
	//An alternative implementation is possible here.
	//Direct range query with the range be a "ring".
	disk_u = alloc_disk( IRTree_v.dim);
	disk_l = alloc_disk( IRTree_v.dim);

	set_disk( disk_u, q->loc_v, UB);
	set_disk( disk_l, q->loc_v, LB);

    //---
    if(obj_set_v == NULL)
    {
        //use the smaller disk
        if(disk_u->radius < disk_v->radius)
            region_u = range_query(disk_u, q);
        else
            region_u = range_query(disk_v, q);
    }else
        region_u = obj_set_v;
    
    //---
    //we must have LB < d here
    //otherwise S_a == NULL
    
	region_l = range_query( disk_l, q);
    refine_region( region_l, disk_l);

	/*t*/
	//print_obj_set( region_l, stdout);
	/*t*/

	obj_exclusion_disk( region_u, disk_l);
    
    UB = UB > disk_v->radius? disk_v->radius : UB;
	R = region_u;

	/*t*/
	//print_obj_set( R, stdout);
	/*t*/

	/*s*/
	//printf( "#cands: %i\n", R->obj_n);
	/*s*/

    
	//Pre-checking.
	if( R->obj_n == 0)
		goto E;

    
	//Sort the objects in R by their distances to q.
	R_heap = heap_sort_obj_set( R, q);
	top = b_h_get_top( R_heap);

	o_next = R_heap->obj_arr[ top].obj_v;

	//Construct the base set of the <o_1, o_2> pairs.
	//region_l should be strict.
	if( s_tag == 1)
	{
		obj_pair = ini_obj_pair( region_l, o_next, LB, UB, q);

		/*s*/
		//printf( "#pairs: %i\n", obj_pair->node_n);
		/*s*/
	}

	//Search.
	while( true)
	{
		if( o_next == NULL)
			break;

		o = o_next;

		top = b_h_get_top( R_heap);
		if( top == 0)
			o_next = NULL;
		else
			o_next = R_heap->obj_arr[ top].obj_v;

		loc_v = get_obj_loc( o);
		dist = calc_dist_loc( loc_v, q->loc_v);
		release_loc( loc_v);

		if( dist > cost_c)
			break;

		/*s*/
		stat_v.n_1_sum ++;
		/*s*/

		//The "ConstructFeasibleSet" procedure.
		if( s_tag == 1)
		{
			//Pre-checking (for the boundary case that |S| = 1).
			S = alloc_obj_set( );
			add_obj_set_entry( o, S);
			if( !is_covered_obj_set( S, q))
			{
				release_obj_set( S);

				if( prune_tag == 1)
				{
					//Update UB by using ConstructFeasibleSet-Appro.
					S = ConstructFeasibleSet_Appro( o, q);
					if( !S)
						continue;
					
					cost = comp_cost(cost_tag, S, q);
					if( cost < cost_c)
					{
						release_obj_set( S_a);
						S_a = S;
						cost_c = cost;
						
						UB = cost_c;
					}
					else
						release_obj_set( S);
				}

				//Update the <o_1, o_2> pairs.
				update_obj_pair( obj_pair, region_l, o, LB, UB, q);

				/*s*/
				//printf( "#pairs: %i\n", obj_pair->node_n);
				/*s*/

				/*t*/
				//in_order_walk( obj_pair->root);
				/*t*/

				S = ConstructFeasibleSet_Exact( o, q, obj_pair);
			}
		}
		else	//s_tag == 2
		{
			S = ConstructFeasibleSet_Appro( o, q);

			/*t*/
			//printf( "base object: %i\n", o->id);
			//print_obj_set( S, stdout);
			/*t*/
		}

		if( !S)
			continue;

		cost = comp_cost(cost_tag, S, q);
        //---
//        printf("cost:%f\n",cost);

        if( cost <= dist_thr)
        {
            //return the set directly without continue checking
            release_obj_set( S_a);
            S_a = S;
            cost_c = cost;
            
            goto E;
        }
        //---
		if( cost < cost_c)
		{
			release_obj_set( S_a);
			S_a = S;
			cost_c = cost;
		}
		else
			release_obj_set( S);
	}//while
E:
	remove_identical_obj( S_a);

	//Release the resource.
	release_disk( disk_u);
	release_disk( disk_l);
    if(obj_set_v == NULL)
        release_obj_set( R); //R==region_u, if obj_set_v not null means to be release in caller
	release_obj_set( region_l);
  	if( s_tag == 1)
		bst_release( obj_pair);
    if(R_heap != NULL)
        release_b_heap( R_heap);

    return S_a;
}

/*
 *	The implementation of the "CostEnum-Exact" algorithm.
 */
obj_set_t* CostEnum_Exact( query_t* q, int prune_tag)
{
	return CostEnum( q, 1, prune_tag, NULL, NULL);
}

/*
 *	The implementation of the "CostEnum-Appro" algorithm.*	
 */
obj_set_t* CostEnum_Appro( query_t* q)
{
	return CostEnum( q, 2, 0, NULL, NULL);
}


