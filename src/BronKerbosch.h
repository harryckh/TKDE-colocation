/*
 * BronKerbosch.h
 *
 *  Created on: 9 Feb 2022
 *      Author: kai-ho
 */

#ifndef BRONKERBOSCH_H_
#define BRONKERBOSCH_H_



#include <vector>
#include <set>
#include <list>
using namespace std;


void BronKerbosch(vector<int> R, vector<int> P, vector<int> X,
		set<vector<int>>* result_cliques, int depth) ;

set<vector<int>> run_BronKerbosch() ;

#endif /* BRONKERBOSCH_H_ */
