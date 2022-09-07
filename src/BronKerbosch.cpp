/*
 * BronKerbosch.cpp
 *
 *  Created on: 9 Feb 2022
 *      Author: kai-ho
 */

#include "BronKerbosch.h"

vector<list<int>> adjacency_list;

void BronKerbosch(vector<int> R, vector<int> P, vector<int> X,
		set<vector<int>> *result_cliques, int depth) {
	/*
	 printf(" %d BronKerbosch\n", depth);

	 printf("R:");
	 for (auto i : R)
	 printf("%d ", i);
	 printf("\n");

	 printf("P:");
	 for (auto i : P)
	 printf("%d ", i);
	 printf("\n");

	 printf("X:");
	 for (auto i : X)
	 printf("%d ", i);
	 printf("\n");
	 */
	if (P.empty() && X.empty()) {
		if (R.size() > 1)
		{
			result_cliques->insert(R);

//		printf("insert R:");
//		for (auto i : R)
//			printf("%d ", i);
//		printf("\n");
		}
	}

	vector<int> P2 = P;
	for (auto node : P2) {
//	for(std::vector<int>::size_type i = 0; i != P.size(); i++)  {
		vector<int> intersectionP = { }, intersectionX = { };

		//N(P)
		for (int nodeP : adjacency_list[node]) {
			for (int node2 : P) {
				if (nodeP == node2) {
					intersectionP.push_back(nodeP);
				}
			}

			//N(X)
			for (int node3 : X) {
				if (nodeP == node3) {
					intersectionX.push_back(nodeP);
				}
			}
		}
		R.push_back(node);
		BronKerbosch(R, intersectionP, intersectionX, result_cliques,
				depth + 1);
		R.pop_back();
//		P.erase(remove(P.begin(), P.end(), node), P.end());
//		vector<int>::iterator position = find(P.begin(), P.end(), node);
//		if (position != P.end())
//			P.erase(position);
		for( vector<int>::iterator iter = P.begin(); iter != P.end(); ++iter )
		{
		    if( *iter == node )
		    {
		        P.erase( iter );
		        break;
		    }
		}
		X.push_back(node);
	}
}

set<vector<int>> run_BronKerbosch() {

	set<vector<int>> result_cliques;
	vector<int> R, P, X;
//	printf("run_BronKerbosch\n");

	//--------
	//gen test
//	for (int i = 0; i < 10; i++) {
//		list<int> tmp;
//		for (int j = 0; j < 10; j++) {
//			if (i != j)
//				tmp.push_back(j);
//		}
//		adjacency_list.push_back(tmp);
//	}
	//--------

	for (int i = 0; i < adjacency_list.size(); i++) {
		P.push_back(i);
	}

	BronKerbosch(R, P, X, &result_cliques, 0);
//	printf("result_cliques.size():%lu\n", result_cliques.size());
//
//	for (auto clique : result_cliques) {
//		for (int node : clique) {
//			printf("%d ", node);
//		}
//		printf("\n");
//	}

	return result_cliques;

}

