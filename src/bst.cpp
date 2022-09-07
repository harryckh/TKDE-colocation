/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */

#include "bst.h"

extern colocation_stat_t stat_v;


bst_t* bst_ini( )
{
	bst_t* T;

	T = ( bst_t*)malloc( sizeof( bst_t));
	memset( T, 0, sizeof( bst_t));

	//Problem specific.
	T->max = - INT_MAX;
	T->min = INT_MAX;

	/*s*/
	stat_v.memory_v += sizeof( bst_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return T;
}

void bst_release_sub( bst_node_t* x)
{
	if( x->left != NULL)
		bst_release_sub( x->left);
	if( x->right != NULL)
		bst_release_sub( x->right);
    
	free( x);

	/*s*/
	stat_v.memory_v -= sizeof( bst_node_t);
	/*s*/
}

/*
 *	Release the binary search tree T.
 */
void bst_release( bst_t* T)
{
	if( T != NULL)
	{
		if( T->root != NULL)
			bst_release_sub( T->root);
		free( T);

		/*s*/
		stat_v.memory_v -= sizeof( bst_t);
		/*s*/
	}
}

void bst_insert( bst_t* T, bst_node_t* z)
{
	bst_node_t* y, *x;

	y = NULL;
	x = T->root;

	while( x != NULL)
	{
		y = x;
		if( z->key < x->key)
			x = x->left;
		else
			x = x->right;
	}

	z->p = y;
	if( y == NULL)
		T->root = z;	//T was empty.
	else if( z->key < y->key)
		y->left = z;	
	else
		y->right = z;

	//
	T->node_n ++;
	if( T->max < z->key)
		T->max = z->key;
	if( T->min > z->key)
		T->min = z->key;
}

void bst_transplant( bst_t* T, bst_node_t* u, bst_node_t* v)
{
	if( u->p == NULL)
		T->root = v;
	else if( u == u->p->left)
		u->p->left = v;
	else
		u->p->right = v;
	if( v != NULL)
		v->p = u->p;
}


bst_node_t* bst_get_min( bst_node_t* x)
{
	if( x == NULL)
		return NULL;

	while( x->left != NULL)
		x = x->left;
	return x;
}

bst_node_t* bst_get_max( bst_node_t* x)
{
	if( x == NULL)
		return NULL;

	while( x->right != NULL)
		x = x->right;
	return x;
}

void bst_delete( bst_t* T, bst_node_t* x)
{
	bst_node_t* y;

	if( x->left == NULL)
		bst_transplant( T, x, x->right);
	else if( x->right == NULL)
		bst_transplant( T, x, x->left);
	else
	{
		y = bst_get_min( x->right);
		if( y->p != x)
		{
			bst_transplant( T, y, y->right);
			y->right = x->right;
			y->right->p = y;
		}

		bst_transplant( T, x, y);
		y->left = x->left;
		y->left->p = y;
	}

	//Problem specific.
	T->node_n --;

	if( T->node_n == 0)
	{
		T->max = -DBL_MAX;
		T->min =  DBL_MAX;

		return;
	}

	if( x->key == T->max)
		T->max = bst_get_max( T->root)->key;
	if( x->key == T->min)
		T->min = bst_get_min( T->root)->key;
}

/*
 *	Update the key of a bst_node @x.
 */
void bst_update( bst_t* T, bst_node_t* x)
{
	bst_delete( T, x);
/*t/
	in_order_walk_non_recur( T->root);
/*t*/
	x->left = NULL;
	x->right = NULL;		//bug.
	bst_insert( T, x);
/*t/
	in_order_walk_non_recur( T->root);
/*t*/
}

void in_order_walk( bst_node_t* x)
{
	if( x != NULL)
	{
		in_order_walk( x->left);
		printf( "%.2lf\t%d\t%d\n", x->key,x->obj_v1->id,x->obj_v2->id);
		in_order_walk( x->right);
	}
}

/*
 *	
 *	@return null if no successor exists.
 */
bst_node_t*  bst_successor( bst_node_t* x)
{
	bst_node_t* y;

	//y = NULL;

	if( x->right != NULL)
		return bst_get_min( x->right);

	y = x->p;
	while( y != NULL && x == y->right)
	{
		x = y;
		y = y->p;
	}

	return y;
}

/*
 * @return null if no predecessor exists.
 */
bst_node_t* bst_predecessor( bst_node_t* x)
{
	bst_node_t* y;

	if( x->left != NULL)
		return bst_get_max( x->left);

	y = x->p;
	while( y != NULL && x == y->left)
	{
		x = y;
		y = y->p;
	}

	return y;
}

void print_bst( bst_t* T)
{
	//
}

void test_bst( )
{
	int i, node_n;
	bst_t* T;
	bst_node_t* x;

	T = bst_ini( );

	printf( "Testing the binary search tree...\n");
	printf( "Please input the number of nodes: ");
	scanf( "%i", &node_n);

	//srand( time( NULL));
	for( i=0; i<node_n; i++)
	{
		x = ( bst_node_t*)malloc( sizeof( bst_node_t));
		memset( x, 0, sizeof( bst_node_t));

		x->key = rand( )%100;
		bst_insert( T, x);
	}

	printf( "Tree constructed...\n");
	printf( "In-order walking...\n");
	in_order_walk( T->root);
	printf( "\n");

	//printf( "Tree constructed...\n");
	printf( "In-order walking (no-recursive) ...\n");
	in_order_walk_non_recur( T->root);
	printf( "\n");

	//test the "get_next_in_order".
	printf( "In-order walking (get_next) ...\n");

	//test bst_search_range.
	int cnt;
	KEY_TYPE low, high;

	low = 900;
	high = 1000;
	
	bst_node_t* low_n, *high_n;
	cnt = bst_search_range( T, low, high, low_n, high_n);

	if( cnt != 0)
	{
		printf( "low_n->key: %.2lf\n", low_n->key); 
		printf( "high_n->key: %.2lf\n", high_n->key);
	}
	printf( "cnt: %i\n\n", cnt);

	int tag;
	x = T->root;
	tag = 0;
	while( get_next_in_order( x, tag))
	{
		printf( "%.2lf  ", x->key);
		if( !in_order_sub( x, tag))
			break;
	}

	//test "delete".
	printf( "\nDelete the root...\n");
	bst_delete( T, T->root);

	printf( "In-order walking...\n");
	in_order_walk( T->root);
	printf( "\n");

	printf( "Retrieving the successor of the root...\n");
	x = bst_successor( T->root);
	printf( "The key of the root is: %.2lf\n", T->root->key);
	if( x != NULL)
		printf( "The key of the successor of the root is: %.2lf\n", x->key);
	else
		printf( "No predecessors!\n");

	printf( "Retrieving the predecessor of the root...\n");
	x = bst_predecessor( T->root);
	printf( "The key of the root is: %.2lf\n", T->root->key);
	if( x != NULL)
		printf( "The key of the predecessor of the root is: %.2lf\n", x->key);
	else
		printf( "No successors!\n");

	printf( "Testing finished!\n");
	
	bst_release( T);
}


bool in_order_sub( bst_node_t* &x, int &tag)
{
	if( x->right != NULL)
	{
		x = x->right;
		tag = 0;
	}
	else
	{
		if( x->p != NULL)
		{
			if( x == x->p->left)
				tag = -1;
			else
				tag = 1;

			x = x->p;
		}
		else
			return false;
	}

	return true;
}

void in_order_walk_non_recur( bst_node_t* x)
{
	int tag;

	if( x == NULL)
		return;

	//Traverse the binary search tree: pre-order, DFS.
	tag = 0;
	while( true)
	{
		if( tag == 0)
		{
			while( x->left != NULL)
				x = x->left;
			
			//
			printf( "%.2lf  ", x->key);

			if( !in_order_sub( x, tag))
				break;
		}
		else if( tag == -1)
		{
			//Backtrack from the left child.
			//Read x itsef.
			//
			printf( "%.2lf  ", x->key);

			if( !in_order_sub( x, tag))
				break;
		}
		else	//tag ==1
		{
			//Backtrack from the right child.
			//Further backtrack the parent.
			if( x->p != NULL)
			{
				if( x == x->p->left)
					tag = -1;
				else
					tag = 1;
				x = x->p;
			}
			else 
				break;
		}
	}

	printf( "\n\n");
}

/*
 *	Locate the "next" bst_node wrt the "in-order" from the current one @x.
 *
 *	return true if successful, otherwise, return false,
 *	in which case, the current bst_node @x is the last bst_node.
 */
bool get_next_in_order( bst_node_t* &x, int &tag)
{
	//tag should be set to 0 when x is the root.
	while( true)
	{
		if( tag == 0)
		{
			while( x->left != NULL)
				x = x->left;
			
			//printf( "%i  ", x->key);
			return true;
		}
		else if( tag == -1)
		{
			//Backtrack from the left child.
			//Read x itsef.
			//printf( "%i  ", x->key);

			return true;
		}
		else	//tag ==1
		{
			//Backtrack from the right child.
			//Further backtrack the parent.
			if( x->p != NULL)
			{
				if( x == x->p->left)
					tag = -1;
				else
					tag = 1;
				x = x->p;
			}
			else 
				return false;
		}
	}
}

/*
 *	Search the bst_node that contans a specific key @key.
 */
bst_node_t* bst_search( bst_t* bst_v, KEY_TYPE key)
{
	bst_node_t* x;

	x = bst_v->root;
	while( true)
	{
		if( x == NULL)
			return x;
		
		if( x->key == key)
			return x;
		else if( key < x->key)
			x = x->left;
		else //key > x->key
			x = x->right;
	}

	return NULL;
}

/*
 *	Trim the bst @bst_v by removing all nodes with the key larger than or equal to @u_bound.
 */
void bst_trim( bst_t* bst_v, KEY_TYPE u_bound)
{
	bst_node_t* x;

	x = bst_get_max( bst_v->root);
	while( x != NULL && x->key >= u_bound)
	{
		bst_delete( bst_v, x);

		free( x);

		/*s*/
		stat_v.memory_v -= sizeof( bst_node_t);
		/*s*/

		x = bst_get_max( bst_v->root);
	}
}

/*
 *	Search the bst based on a key range [@low, @high].
 *
 *	@low_n retrieves the first node that has key >= low.
 *	@high_n retrieves the last node that has key <= high.
 *
 *	In case that no key values in bst fall in range [low, high], return 0,
 *	otherwise set low_n and high_n properly
 *	and return the number of key values falling in range [low, high].
 */
int bst_search_range( bst_t* bst_v, KEY_TYPE low, KEY_TYPE high, bst_node_t* &low_n, bst_node_t* &high_n)
{
	int cnt;
	bst_node_t* x_1, *x_2;

	if( low > high)
		return 0;

	x_1 = bst_get_min( bst_v->root);
	if( x_1->key > high)
		return 0;
	
	x_2 = bst_get_max( bst_v->root);
	if( x_2->key < low)
		return 0;
		
	//Set low_n.
	if( x_1->key >= low && x_1->key <= high)
		low_n = x_1;
	else
	{
		while( x_1->key < low)
			x_1 = bst_successor( x_1);
		if( x_1->key > high)		//bug.
			return 0;

		low_n = x_1;
	}

	//Set high_n.
	if( x_2->key <= high)
		high_n = x_2;
	else
	{
		while( x_2->key > high)
			x_2 = bst_predecessor( x_2);
		high_n = x_2;
	}

	//Counting.
	cnt = 1;
	while( x_1 != x_2)
	{
		x_1 = bst_successor( x_1);
		cnt ++;
	}

	return cnt;
}

/*
 *	Search those nodes with the "perecent ranks" falling in range [low, high].
 */
int bst_search_percentile_range( bst_t* bst_v, int low, int high, bst_node_t* &low_n, bst_node_t* &high_n)
{
	if( low < 0 || low > 100)
		return 0;
	if( high < 0 || high > 100)
		return 0;
	if( low > high)
		return 0;

	int cnt_i;
	float cnt_f, low_p, high_p;
	bst_node_t* x;

	low_p = low / 100.0;
	high_p = high / 100.0;

	x = bst_get_min( bst_v->root);
	if( x == NULL)
		return 0;

	cnt_f = 1;
	while( cnt_f / bst_v->node_n < low_p && x != NULL)
	{
		x = bst_successor( x);
		cnt_f ++;
	}

	if( x == NULL)
		return 0;

	low_n = x;
	cnt_i = 1;

	while( cnt_f / bst_v->node_n < high_p && x != NULL)
	{
		x =  bst_successor( x);
		cnt_f ++;
		cnt_i ++;
	}

	high_n = x;

	return cnt_i;
}

/*
 *	Locate the bst_node that contains the specific keyword @key_id.
 *
 *	Note that bst_locate is different from bst_search.
 *		-bst_search is based on the key, O(log(n)).
 *		-bst_locate is based on a non-key, O(n).
 */
bst_node_t* bst_locate( bst_t* bst_v, KEY_TYPE key_id)
{
	int tag;
	bst_node_t* x;

	if( bst_v->root == NULL)
		return NULL;

	x = bst_v->root;
	tag = 0;
	while( get_next_in_order( x, tag))
	{
		if( x->key_id == key_id)
			return x;

		if( !in_order_sub( x, tag))
			return NULL;
	}

	return NULL;
}

