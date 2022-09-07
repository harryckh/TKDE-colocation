
/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */


/*
 *	1. The range of k is [0, 31].
 */
#include "bit_operator.h"

#define	BIT_LENGTH	32

/*
 *	Retrieve the k^th bit in the bit array represented by an BIT_TYPE number v.
 */	
int get_k_bit( BIT_TYPE v, int k)
{
	k++;

	if( !( k >= 1 && k <= 32))
		return -1;
	
	//Divide v by 2^(32-k).
	v = v >> ( 32 - k);
	
	return v % 2;
}

/*
 *	Set the k^th bit of the bit array represented by v to be 1.
 */
void insert_k_bit( BIT_TYPE &v, int k)
{
	if( !( get_k_bit( v, k)))
		v += ( int)pow( (float)2, int( 32 - k - 1));		
}

/*
 *	Set the k^th bit of the bit array represented by v to be 0.
 */
void delete_k_bit( BIT_TYPE& v, int k)
{
	if( ( get_k_bit( v, k)))
		v -= ( int)pow( ( float)2, int( 32 - k - 1));
}

/*
 *	Union the bits of two unsigned integers.
 */
void union_bit( BIT_TYPE &v1, BIT_TYPE v2)
{
	int i;

	for( i=0; i<BIT_LENGTH; i++)
	{
		if( get_k_bit( v2, i))
			insert_k_bit( v1, i);
	}
}

