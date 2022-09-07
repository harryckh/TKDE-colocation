
/*
 * mCK in SIGMOD15
 */

#ifndef MCK_H
#define	MCK_H

#include "costenum.h"
#include "frac.h"
#include "SmallestEnclosingCircle.h"
#include <limits>


extern int	cost_tag;
extern float dist_thr;

extern IRTree_t IRTree_v;
extern bst_t* IF_v;

//==================== algorithms =================================================

obj_set_t* mck_Exact( query_t* q, B_KEY_TYPE eps, obj_set_t * region_u, disk_t* disk_v);

void search( query_t* q, obj_t* o, obj_t* o1, obj_set_t* G, obj_set_t*& curGroup, B_KEY_TYPE& cost_curGroup);

void search_sub( query_t* q, obj_set_t* selectedSet, obj_set_t* candidateSet, int maxId, B_KEY_TYPE& cost_selectedSet, psi_t* psi_uncover_selectedSet, obj_set_t*& curGroup, B_KEY_TYPE& cost_curGroup);

//=============================================================================

obj_set_t* GKG( query_t* q, obj_set_t* obj_set_v, disk_t* disk_v);

disk_t* MCC1(obj_set_t* obj_set_v);

disk_t* MCC2(obj_set_t* obj_set_v);

disk_t*  makeCircleOnePoint(obj_set_t* obj_set_v, obj_node_t* p) ;
disk_t* makeDiameter(obj_node_t* p1, obj_node_t* p2) ;
disk_t* makeCircleTwoPoints(obj_set_t* obj_set_v, obj_node_t* p, obj_node_t* q);
loc_t* cir_subtract(loc_t* loc_p, obj_node_t* q);
loc_t* subtract(obj_node_t* p, obj_node_t* q);
B_KEY_TYPE calc_cross(loc_t* p, loc_t* q);
disk_t* makeCircumcircle(obj_node_t* p1, obj_node_t* p2, obj_node_t* p3);

bool is_inside(obj_t* obj_v, disk_t* disk_v);
//===============================================================================

obj_set_t* SKEC( query_t* q, int s_tag, B_KEY_TYPE eps, obj_set_t* region_u, disk_t* disk_v);

void findOSKEC(obj_set_t* O_t, obj_node_t* obj_node_v, disk_t*& c_cur, query_t* q, obj_set_t*& S_cur);

void findAppOSKEC(obj_set_t* O_t, obj_node_t* o, B_KEY_TYPE& c_cur_diam, query_t* q, B_KEY_TYPE eps, B_KEY_TYPE LB, obj_set_t*& cur_group, B_KEY_TYPE& cost_curGroup);

//===============================================================================

obj_set_t* SKECaplus( query_t* q, B_KEY_TYPE eps, bool callByExact, obj_set_t*& forExact, disk_t* disk_v);

bool circleScan(obj_set_t* region_u, obj_t* o, B_KEY_TYPE diam, query_t* q, int s_tag, obj_set_t*& curGroup, B_KEY_TYPE& cost_curGroup);

void getAngle( obj_node_t* o, obj_node_t* o_j, B_KEY_TYPE r, B_KEY_TYPE& o_in, B_KEY_TYPE& o_out);

void getAngle2( obj_t* o, obj_t* o_j, B_KEY_TYPE r, B_KEY_TYPE& o_in, B_KEY_TYPE& o_out);

//=============== Tab related ==============================================

bst_t* const_tab( psi_t* psi_v);

bool is_full(bst_node_t* bst_node_v);

void add_Tab_key(bst_t* Tab, k_node_t* k_head);

void add_Tab_key(bst_t* Tab, FEA_TYPE fea);

void remove_Tab_key(bst_t* Tab, k_node_t* k_head);

void remove_Tab_key(bst_t* Tab, FEA_TYPE fea);

void print_Tab(bst_t* Tab);

//================ other functions =====================================================

void add_obj_set_entry_extended2( obj_t* obj_v, obj_set_t* obj_set_v, B_KEY_TYPE angle, bool type);

b_heap_t* heap_sort_obj_set( obj_set_t* obj_set_v);

bool is_contained_obj_set( obj_set_t* obj_set_v, obj_t* obj_v);


bool is_covered_obj_set( obj_set_t* obj_set_v, query_t* q, B_KEY_TYPE d);

//================ query related =====================================================


obj_set_t* range_query3();

void retrieve_sub_tree3( node_t* node_v, obj_set_t* &obj_set_v);

query_t** gen_query_set3( int query_n, int key_n, range* MBR, data_t* data_v, int low, int high, B_KEY_TYPE maxDia);

//from other files

obj_set_t* Cao_Appro1( query_t* q);

obj_set_t* Cao_Appro1( query_t* q, disk_t* disk_v);

void remove_obj_set_entry( obj_set_t* obj_set_v, obj_t* obj_v );

psi_t*  uncover_keyword(query_t* q, obj_set_t* obj_set_v);

int is_relevant_obj( obj_t* obj_v, psi_t* psi_v);

void psi_exclusion( psi_t* psi_v1, obj_t* obj_v);

#endif
