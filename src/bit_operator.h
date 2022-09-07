
/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */

/*
 *	1. The bitwise operations are based on the "BIT_TYPE type",
 *	which has 4 bytes.
 *
 *	2. It is possible to extend the number of bits to 64 bits.
 */


#ifndef BIT_OPERATOR_H
#define BIT_OPERATOR_H


#include <math.h>
//#include <stdint.h>

#define BIT_TYPE	unsigned int


int get_k_bit( BIT_TYPE v, int k);

void insert_k_bit( BIT_TYPE &v, int k);

void delete_k_bit( BIT_TYPE &v, int k);

void union_bit( BIT_TYPE &v1, BIT_TYPE v2);

#endif
