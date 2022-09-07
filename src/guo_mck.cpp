/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */

#include "guo_mck.h"

//the disk will enlarge radius by this value
double disk_r_tolerance = 0.001;


//================================================================================
//============== Exact mCK =======================================================
//================================================================================


obj_set_t* mck_Exact( query_t* q, B_KEY_TYPE eps, obj_set_t * region_u, disk_t* disk_v2)
{
    
    obj_set_t* bestGroup, * G_skeca;
    B_KEY_TYPE diam, cost_curGroup;
    disk_t* disk_v;
    obj_node_t* o;
    
    G_skeca = SKECaplus( q, eps, true, region_u, disk_v2);
    
    if(G_skeca == NULL)
    {
        return NULL;
    }
    
    bestGroup = G_skeca;
    cost_curGroup = comp_diameter(bestGroup);
    
    if(cost_curGroup <= dist_thr)
    {
        return bestGroup;
    }
    
    disk_v = MCC(G_skeca);
    diam =(2/sqrt(3)) * disk_v->radius * 2;
    
    // printf("|region_u|:%d\n",region_u->obj_n);
    
    o = region_u->head->next;
    while ( o!= NULL){
        if (o->maxInvalidRange < diam)
        {
            circleScan(region_u, o->obj_v, diam, q, 1, bestGroup, cost_curGroup);
            if(cost_curGroup <= dist_thr)
                break;
        }
        o=o->next;
    }
    
    //release resources
    //    o = region_u->head->next;
    //    while ( o!= NULL){
    //        release_obj_set(o->range);
    //        o = o->next;
    //    }
    
    release_disk(disk_v);
    
    return bestGroup;
}


//==================================================================================
//==================================================================================
//==================================================================================

// GKG implemtnation
// for each object contain t_inf, run Cao-Appro1 to find NN set
obj_set_t* GKG( query_t* q, obj_set_t* obj_set_v, disk_t* disk_v)
{
    obj_set_t* S_c;				//the current best solution.
    obj_set_t* S;               //the newly constructed feasible set.
    B_KEY_TYPE cost_c, cost;
    //bst_t* IF_v;
    bst_node_t* bst_node_v, *bst_node_min;
    obj_node_t* obj_node_v;
    int min;
    obj_t* obj_v;
    query_t* q_new;
    k_node_t* k_node_v;
    bst_t* inverted_list;
    
    // ----------------------------------------------
    //find the most infrequent query keyword
    inverted_list = const_IF(obj_set_v);
    
    min = INT_MAX;
    bst_node_min = NULL;
    bst_node_v = bst_get_min(inverted_list->root);
    
    while (bst_node_v != NULL){
        
        if( bst_node_v->p_list_obj->obj_n < min)
        {
            min = bst_node_v->p_list_obj->obj_n;
            bst_node_min = bst_node_v;
        }
        bst_node_v = bst_successor( bst_node_v );
    }
    // bst_node_min = node contain min number of obj
    // ----------------------------------------------
    //for each object contain t, construct a new query with different cost_tag
    cost_c = INFINITY;
    S_c = NULL;
    
    obj_node_v = bst_node_min->p_list_obj->head->next;
    while( obj_node_v != NULL){
        
        obj_v = obj_node_v->obj_v;
        
        //    printf("obj_v id:%d\n",obj_v->id);
        
        //construct new query
        q_new = alloc_query();
        q_new->loc_v = alloc_loc(q->loc_v->dim);
        for( int j=0; j<q_new->loc_v->dim; j++)
        {
            q_new->loc_v->coord[j] = obj_v->MBR[j].min ;
        }
        //        k_node_v = key_exclusion(q->psi_v->k_head, obj_v->k_head);
        k_node_v = key_exclusion(q->psi_v->k_head, obj_v->fea);
        q_new->psi_v = const_psi( k_node_v );
        
        //perform new query
        S = Cao_Appro1(q_new, disk_v);
        
        release_query(q_new);
        
        //add the object into the set
        add_obj_set_entry(obj_v, S);
        
        //update current best

        cost = comp_cost(cost_tag, S, q);
        //--
        if( cost <= dist_thr)
        {
            //return the set directly without continue checking
            release_obj_set( S_c);
            S_c = S;
            cost_c = cost;
            goto E;
        }
        //--
        if ( cost < cost_c){
            release_obj_set(S_c);
            S_c = S;
            cost_c = cost;
            S = NULL;
        }else
            release_obj_set(S);
        
       
        obj_node_v = obj_node_v->next;
    }
    
E:
    //bst_release(inverted_list);
    release_IF(inverted_list);
    
    return S_c;
}




/*
 input: an object set
 output: the smallest circle that cover all object in the object set
 * O(n^3) implementation: try all possible circle formed by 2 or 3 objects
 */
disk_t* MCC1(obj_set_t* obj_set_v)
{
    
    disk_t* circle_v;
    double x1,y1,x2,y2,x3,y3,x,y,r;
    loc_t* loc_v;
    disk_t* disk_v;
    
    bool flag;
    double x_min, y_min, r_min = DBL_MAX;
    
    
    //-----------------
    //minus all coordinate  by the first object before we start
//    double x0 = obj_set_v->head->next->obj_v->MBR[0].min;
//    double y0 = obj_set_v->head->next->obj_v->MBR[1].min;
//    for(obj_node_t * p1=obj_set_v->head->next;p1!=NULL;p1=p1->next)
//    {
//        p1->obj_v->MBR[0].min -= x0;
//        p1->obj_v->MBR[0].max -= x0;
//        p1->obj_v->MBR[1].min -= y0;
//        p1->obj_v->MBR[1].max -= y0;
//        
//    }
    //-----------------
    
    
    // printf("obj_set_v:%d\n",obj_set_v->obj_n);
    
    //try all triple / double objs and return the min.
    for(obj_node_t * p1=obj_set_v->head->next;p1->next!=NULL;p1=p1->next)
    {
        for(obj_node_t * p2=p1->next;p2->next!=NULL;p2=p2->next)
        {
            for(obj_node_t * p3=p2->next;;p3=p3->next)
            {
                if(p3==NULL)break;
                
                
                x1=p1->obj_v->MBR[0].min;
                y1=p1->obj_v->MBR[1].min;
                x2=p2->obj_v->MBR[0].min;
                y2=p2->obj_v->MBR[1].min;
                x3=p3->obj_v->MBR[0].min;
                y3=p3->obj_v->MBR[1].min;
                
                
              
//                printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x1,y1,x2,y2,x3,y3);
                //checked, correct
                x =   ((pow(x1,(double)2.0) + pow(y1,(double)2.0)) * (y2-y3)
                       + (pow(x2,(double)2.0) + pow(y2,(double)2.0)) * (y3-y1)
                       + (pow(x3,(double)2.0) + pow(y3,(double)2.0)) * (y1-y2))
                / ( 2* (x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
                
                y =  ((pow(x1,(double)2.0)+pow(y1,(double)2.0))*(x3-x2)
                      + (pow(x2,(double)2.0)+pow(y2,(double)2.0))*(x1-x3)
                      + (pow(x3,(double)2.0)+pow(y3,(double)2.0))*(x2-x1))
                / ( 2 * (x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
                
                r = fmax(fmax(sqrt(pow(x-x1,(double)2.0)+pow(y-y1,(double)2.0)),
                              sqrt(pow(x-x2,(double)2.0)+pow(y-y2,(double)2.0))),
                              sqrt(pow(x-x3,(double)2.0)+pow(y-y3,(double)2.0)));

                //enlarge the tolerance proportional to x,y
                r += fmax(fmax( x, y)/100000, 1.0) * disk_r_tolerance;
     
//                                printf("%f\t%f\n",x1,y1);
//                                printf("%f\t%f\n",x2,y2);
//                                printf("%f\t%f\n",x3,y3);
//                                printf("%f\t%f\t%f\n",x,y, r);
                
                //----------
                loc_v = alloc_loc(2);
                loc_v->coord[0] = x;
                loc_v->coord[1] = y;
                disk_v = const_disk(loc_v, r);
                
                release_loc(loc_v);
                
                flag = true;
                //check whether all object insides the disk
                for(obj_node_t * p4=obj_set_v->head->next;p4!=NULL;p4=p4->next)
                {
                    if(!is_inside(p4->obj_v, disk_v))
                    {
                        flag = false;
                        break;
                    }
                }
                //flag == true now mean yes
                release_disk(disk_v);
                if(flag && r < r_min){
                    //the circle enclose all obj
                    r_min = r;
                    x_min = x;
                    y_min = y;
                    
                }
                //-------------
            }
        }
    }
    
    //-----------------------------
    for(obj_node_t * p1=obj_set_v->head->next;p1->next!=NULL;p1=p1->next){
        for(obj_node_t * p2=p1->next;;p2=p2->next){
            if(p2==NULL)break;
            
            x1=p1->obj_v->MBR[0].min;
            y1=p1->obj_v->MBR[1].min;
            x2=p2->obj_v->MBR[0].min;
            y2=p2->obj_v->MBR[1].min;
            
//            x=(x1+x2)/2.0;
//            y=(y1+y2)/2.0;
//            r = fmax(sqrt(pow(x-x1,(double)2.0)+pow(y-y1,(double)2.0)),
//                     sqrt(pow(x-x2,(double)2.0)+pow(y-y2,(double)2.0)));
            
            x= (p1->obj_v->MBR[0].min + p2->obj_v->MBR[0].min)/2.0;
            y= (p1->obj_v->MBR[1].min + p2->obj_v->MBR[1].min)/2.0;
            r = fmax(sqrt(pow(x-p1->obj_v->MBR[0].min,(double)2.0)+pow(y-p1->obj_v->MBR[1].min,(double)2.0)),
                   sqrt(pow(x-p2->obj_v->MBR[0].min,(double)2.0)+pow(y-p2->obj_v->MBR[1].min,(double)2.0)));
            
            //enlarge the tolerance proportional to x,y
            r += fmax(fmax( x, y)/100000, 1.0) * disk_r_tolerance;
            
//            printf("%f\t%f\n",x1,y1);
//            printf("%f\t%f\n",x2,y2);
//            printf("%f\t%f\t%f\n",x,y, r);

            //----------
            loc_v = alloc_loc(2);
            loc_v->coord[0] = x;
            loc_v->coord[1] = y;
            disk_v = const_disk(loc_v, r);
            
            release_loc(loc_v);
            
            flag = true;
            
            //check whether all object insides the disk
            for(obj_node_t * p4=obj_set_v->head->next;p4!=NULL;p4=p4->next){
                if(!is_inside(p4->obj_v, disk_v))
                {
                    flag = false;
                    break;
                }
            }
            //flag == true now mean yes
            release_disk(disk_v);
            if(flag && r < r_min){
                //the circle enclose all obj
                r_min = r;
                x_min = x;
                y_min = y;
                
            }
            //------------
        }
    }
    //-------------------
    
    circle_v = alloc_disk(2);
    circle_v->center->coord[0] = x_min;
    circle_v->center->coord[1] = y_min;
//    circle_v->center->coord[0] = x_min + x0;
//    circle_v->center->coord[1] = y_min + y0;
    circle_v->radius = r_min;
    
    //printf("cir r:%f\n",r_min);
    
    return circle_v;
}



/*
 input: an object set
 output: the smallest circle that cover all object in the object set
 */
disk_t* MCC2(obj_set_t* obj_set_v){
    
    disk_t* circle_v = NULL;
    obj_node_t* obj_node_v;
    
    //-------
    //minus all coordinate  by the first object before we start
    double x0 = obj_set_v->head->next->obj_v->MBR[0].min;
    double y0 = obj_set_v->head->next->obj_v->MBR[1].min;
    for(obj_node_t * p1=obj_set_v->head->next;p1!=NULL;p1=p1->next)
    {
        p1->obj_v->MBR[0].min -= x0;
        p1->obj_v->MBR[0].max -= x0;
        p1->obj_v->MBR[1].min -= y0;
        p1->obj_v->MBR[1].max -= y0;
        
    }
    //-----------------
    

    
    
    //if only 2 objs, the cir is center = mid-point
    if(obj_set_v->obj_n == 2)
    {
        circle_v = makeDiameter(obj_set_v->head->next, obj_set_v->head->next->next);
        return circle_v;
    }
    
    
    //insert obj one by one
    obj_node_v = obj_set_v->head->next;
    while (obj_node_v != NULL){
        //if the cur circle cannot enclose obj, a new cir is needed
        //obj must be on the boundary of the new cir
        if ( circle_v == NULL || !is_inside(obj_node_v->obj_v, circle_v))
        {
            release_disk(circle_v);
            circle_v = makeCircleOnePoint(obj_set_v, obj_node_v);
        }
        
        obj_node_v=obj_node_v->next;
    }
    
    return circle_v;
}

// One boundary point known
disk_t*  makeCircleOnePoint(obj_set_t* obj_set_v, obj_node_t* p) {
    
    disk_t* cir;
    obj_node_t* obj_node_v;
    
    cir = alloc_disk(2);
    cir->center->coord[0] = p->obj_v->MBR[0].min;
    cir->center->coord[1] = p->obj_v->MBR[1].min;
    cir->radius = 0.0;
    
    
    //   printf("\t 1:%d\n",p->obj_v->id);
    
    //for each obj inserted before p
    obj_node_v = obj_set_v->head->next;
    while (obj_node_v != NULL && obj_node_v != p->next){
        //     printf("\t \t 2:%d\n",obj_node_v->obj_v->id);
        
        //if the cur cir (that with p in boundary) cannot enclose obj, a new cir is needed
        //p and obj must be on the boundary
        
        //if ( !is_enclosed(obj_node_v->obj_v->MBR, cir))
        if ( !is_inside(obj_node_v->obj_v, cir))
        {
            if (cir->radius==0.0){
                release_disk(cir);
                cir = makeDiameter(p, obj_node_v);
            }else{
                release_disk(cir);
                cir = makeCircleTwoPoints(obj_set_v, p, obj_node_v);
            }
        }
        obj_node_v=obj_node_v->next;
    }
    return cir;
}

//construct a circle with two points (diameter)
disk_t* makeDiameter(obj_node_t* p1, obj_node_t* p2) {
    
    disk_t* cir;
    
    
//    double x1 = p1->obj_v->MBR[0].min-p1->obj_v->MBR[0].min;
//    double y1 = p1->obj_v->MBR[1].min-p1->obj_v->MBR[1].min;
//    double x2 = p2->obj_v->MBR[0].min-p1->obj_v->MBR[0].min;
//    double y2 = p2->obj_v->MBR[1].min-p1->obj_v->MBR[1].min;
    
    
    cir = alloc_disk(2);
    cir->center->coord[0] = (p1->obj_v->MBR[0].min + p2->obj_v->MBR[0].min)/2;
    cir->center->coord[1] = (p1->obj_v->MBR[1].min + p2->obj_v->MBR[1].min)/2;
    
//    cir->center->coord[0] = (x1+x2)/2.0;
//    cir->center->coord[1] = (y1+y2)/2.0;

//    cir->radius = fmax(sqrt(pow(cir->center->coord[0]-x1,(double)2.0)+
//                            pow(cir->center->coord[1]-y1,(double)2.0)),
//                       sqrt(pow(cir->center->coord[0]-x2,(double)2.0)+
//                            pow(cir->center->coord[1]-y2,(double)2.0)));

  
//    printf("dia r:%lf\t%lf\n",sqrt(pow(cir->center->coord[0]-p1->obj_v->MBR[0].min,(double)2.0)+
//                pow(cir->center->coord[1]-p1->obj_v->MBR[1].min,(double)2.0)),
//           sqrt(pow(cir->center->coord[0]-p2->obj_v->MBR[0].min,(double)2.0)+
//                pow(cir->center->coord[1]-p2->obj_v->MBR[1].min,(double)2.0)));
    cir->radius = fmax(sqrt(pow(cir->center->coord[0]-p1->obj_v->MBR[0].min,(double)2.0)+
                            pow(cir->center->coord[1]-p1->obj_v->MBR[1].min,(double)2.0)),
                       sqrt(pow(cir->center->coord[0]-p2->obj_v->MBR[0].min,(double)2.0)+
                            pow(cir->center->coord[1]-p2->obj_v->MBR[1].min,(double)2.0)));
    
    
    //enlarge the tolerance proportional to x,y
    cir->radius += fmax(fmax( cir->center->coord[0], cir->center->coord[1])/100000, 1.0) * disk_r_tolerance;
    
//    cir->center->coord[0] += p1->obj_v->MBR[0].min;
//    cir->center->coord[1] += p1->obj_v->MBR[1].min;
//    printf("%lf\t%lf\t%lf\n",cir->center->coord[0],cir->center->coord[1], cir->radius);
    
    return cir;
}


// Two boundary points known
disk_t* makeCircleTwoPoints(obj_set_t* obj_set_v, obj_node_t* p, obj_node_t* q) {
    
    obj_node_t* obj_node_v;
    disk_t *temp, * cir, * left = NULL, * right=NULL;
    bool flag = true;
    
    // printf("\t\tp:%d\t q:%d\t \n",p->obj_v->id,q->obj_v->id);
    
    //---try whether cir by p, q can cover all objs
    temp = makeDiameter(p, q);
    obj_node_v = obj_set_v->head->next;
    while (obj_node_v != NULL && obj_node_v != q->next)
    {
         if (!is_inside(obj_node_v->obj_v, temp))
        {
            flag = false;
            break;
        }
        obj_node_v=obj_node_v->next;
    }
    if (flag)
        return temp;
    release_disk(temp);
    //---
    //for each obj inserted before q
    obj_node_v = obj_set_v->head->next;
    while (obj_node_v != NULL && obj_node_v != q->next){
        loc_t* cp,* lp = NULL,* rp=NULL, * qp, * loc_v;
        //      printf("p:%d\t q:%d\t r:%d\n",p->obj_v->id,q->obj_v->id,obj_node_v->obj_v->id);
        qp = subtract(q,p);
        loc_v = subtract(obj_node_v, p);
        B_KEY_TYPE cross = calc_cross(qp, loc_v);
        
        
        disk_t* c = makeCircumcircle(p, q, obj_node_v);
        
        
        if ( c == NULL){
            obj_node_v=obj_node_v->next;
            release_loc(qp);
            release_loc(loc_v);
            
            continue;
        }
        
        cp = cir_subtract(c->center, p);
        if(left != NULL)
            lp = cir_subtract(left->center, p);
        if(right != NULL)
            rp = cir_subtract(right->center, p);
        
        if (cross > 0 && (left == NULL || calc_cross(qp, cp) > calc_cross(qp, lp))){
            release_disk(left);
            left = c;
        }else if (cross < 0 && (right == NULL || calc_cross(qp, cp)  < calc_cross(qp, rp))){
            release_disk(right);
            right = c;
        }
        
        release_loc(cp);
        if( lp != NULL)
            release_loc(lp);
        if( rp != NULL)
            release_loc(rp);
        
        release_loc(qp);
        release_loc(loc_v);
        
        obj_node_v=obj_node_v->next;
    }
    
    if (right == NULL || (left != NULL && left->radius <= right->radius)){
        
        cir = left;
        release_disk(right);
    }else{
        
        cir = right;
        release_disk(left);
    }
    return cir;
}

//p-q
loc_t* cir_subtract(loc_t* loc_p, obj_node_t* q){
    
    loc_t* loc_v, * loc_q;
    loc_v = alloc_loc(2);
    loc_q = get_obj_loc(q->obj_v);
    loc_v->coord[0] = loc_p->coord[0] - loc_q->coord[0];
    loc_v->coord[1] = loc_p->coord[1] - loc_q->coord[1];
    
    release_loc(loc_q);
    
    return loc_v;
}

//p-q
loc_t* subtract(obj_node_t* p, obj_node_t* q){
    
    loc_t* loc_v, * loc_p, * loc_q;
    
    loc_v = alloc_loc(2);
    loc_p = get_obj_loc(p->obj_v);
    loc_q = get_obj_loc(q->obj_v);
    
    loc_v->coord[0] = loc_p->coord[0] - loc_q->coord[0];
    loc_v->coord[1] = loc_p->coord[1] - loc_q->coord[1];
    
    release_loc(loc_p);
    release_loc(loc_q);
    
    return loc_v;
}

B_KEY_TYPE calc_cross(loc_t* p, loc_t* q){
    return (p->coord[0] * q->coord[1]) - (p->coord[1] * q->coord[0]);
}


//construct circle with three points
disk_t* makeCircumcircle(obj_node_t* p1, obj_node_t* p2, obj_node_t* p3) {
    
    disk_t* cir;
    B_KEY_TYPE x,y,x1,x2,x3,y1,y2,y3,r,d;
    
    x1=p1->obj_v->MBR[0].min-p1->obj_v->MBR[0].min;
    y1=p1->obj_v->MBR[1].min-p1->obj_v->MBR[1].min;
    x2=p2->obj_v->MBR[0].min-p1->obj_v->MBR[0].min;
    y2=p2->obj_v->MBR[1].min-p1->obj_v->MBR[1].min;
    x3=p3->obj_v->MBR[0].min-p1->obj_v->MBR[0].min;
    y3=p3->obj_v->MBR[1].min-p1->obj_v->MBR[1].min;
    
    d = (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) * 2;
    if (d == 0)
        return NULL;
    x = ((x1*x1 + y1*y1) * (y2 - y3) + (x2*x2 + y2*y2) * (y3 - y1) + (x3*x3 + y3*y3) * (y1 - y2))/d;
    y = ((x1*x1 + y1*y1) * (x3 - x2) + (x2*x2 + y2*y2) * (x1 - x3) + (x3*x3 + y3*y3) * (x2 - x1))/d;
    
    r = fmax(fmax(sqrt(pow(x-x1,(double)2.0)+pow(y-y1,(double)2.0)),
                  sqrt(pow(x-x2,(double)2.0)+pow(y-y2,(double)2.0))),
             sqrt(pow(x-x3,(double)2.0)+pow(y-y3,(double)2.0)));
    
    //enlarge the tolerance proportional to x,y
    r += fmax(fmax( x, y)/100000, 1.0) * disk_r_tolerance;

    //--------------------
    
    cir = alloc_disk(2);
    cir->center->coord[0] = x+p1->obj_v->MBR[0].min;
    cir->center->coord[1] = y+p1->obj_v->MBR[1].min;
    cir->radius = r;
    
    return cir;
    
}

//return true if the @obj_v is located inside the disk
bool is_inside(obj_t* obj_v, disk_t* disk_v)
{
//    printf("is_inside:%lf\n",sqrt(pow(obj_v->MBR[0].min-disk_v->center->coord[0],(double)2.0)+
//                                  pow(obj_v->MBR[1].min-disk_v->center->coord[1],(double)2.0)));
    return sqrt(pow(obj_v->MBR[0].min-disk_v->center->coord[0],(double)2.0)+
                pow(obj_v->MBR[1].min-disk_v->center->coord[1],(double)2.0)) <= disk_v->radius;
    

}

//===============================================================================
//===============================================================================
//===============================================================================
//SKEC
//s_tag==1 SKEC
//s_tag==2 SKECa
obj_set_t* SKEC( query_t* q, int s_tag, B_KEY_TYPE eps, obj_set_t* region_u, disk_t* disk_v)
{
    obj_set_t* S_cur;
    obj_node_t *o;
    disk_t* c_cur;
    obj_set_t* G_gkg, * G_skec = NULL;
    B_KEY_TYPE LB, c_cur_dia;
    double cost_curGroup = INFINITY;
    
    G_gkg = GKG(q, region_u, disk_v);
    c_cur = MCC(G_gkg);
    S_cur = copy_obj_set(G_gkg);
    
    LB = comp_cost(cost_tag, G_gkg, q) / (double) 2.0;  //S dia / 2
    c_cur_dia = c_cur->radius * 2;
    
    if (s_tag==2)
        G_skec = alloc_obj_set();
    
    //  region_u = range_query(q);
    
    //  printf("|region_u|:%d\n",region_u->obj_n);
    
    //---------------------------------------------------
    ///findOSKEC
    o = region_u->head->next;
    while ( o!= NULL){
        if(s_tag==1){
            
            G_skec = alloc_obj_set();
            add_obj_set_entry(o->obj_v, G_skec);
            if (is_covered_obj_set(G_skec, q))
                goto E;
            release_obj_set(G_skec);
            
            //c_cur is updated in this function
            findOSKEC(region_u, o, c_cur, q, S_cur);
        }
        else
        {//s_tag==2
            
            findAppOSKEC(region_u, o, c_cur_dia, q, eps, LB, G_skec, cost_curGroup);
        }
        o=o->next;
    }
    //---------------------------------------------------
    //G_skec <- objects in c_cur
    if (s_tag==1)
        //G_skec = range_query(c_cur, q);
        G_skec = S_cur;
E:
    release_obj_set(G_gkg);
    release_obj_set( region_u );
    release_disk(c_cur);
    
    return G_skec;
}


/*
 find the SKEC
 c_cur will update if better is found
 */
void findOSKEC(obj_set_t* O_t, obj_node_t* obj_node_v, disk_t*& c_cur, query_t* q, obj_set_t*& S_cur){
    
    obj_t* obj_v1,* obj_v2;
    loc_t* loc_v, * loc_v1, * loc_v2;
    B_KEY_TYPE dist_o_1, dist_o_2, dist_1_2;
    disk_t* c_can;
    obj_set_t* obj_set_v, * olist, * S;
    obj_node_t* obj_node_v1, * obj_node_v2, * obj_node_new;
    obj_node_t* obj_node_temp;
    psi_t* psi_v;
    
    
    loc_v = get_obj_loc( obj_node_v->obj_v);
    
    olist = alloc_obj_set();
    //-----------------------------
    obj_node_v1 = O_t->head->next;
    while( obj_node_v1 != NULL)
    {
        loc_v1 = get_obj_loc( obj_node_v1->obj_v);
        dist_o_1 = calc_dist_loc( loc_v1, loc_v);
        release_loc(loc_v1);
        
        if (obj_node_v1->obj_v == obj_node_v->obj_v ||
            dist_o_1 > 2 * c_cur->radius){
            obj_node_v1 = obj_node_v1->next;
            continue;
        }
        
        //add to olist while keeping it sorted
        //find the correct position
        obj_node_v2 = olist->head;
        while (obj_node_v2->next != NULL){
            
            if ( dist_o_1 < obj_node_v2->next->dist )
                break;
            obj_node_v2 = obj_node_v2->next;
            
        }
        //obj_node_v2 pointing to the node before insert position
        //-----
        obj_node_new = ( obj_node_t*)malloc( sizeof( obj_node_t));
        memset( obj_node_new, 0, sizeof( obj_node_t));
        
        /*s*/
        stat_v.memory_v += sizeof( obj_node_t);
        if( stat_v.memory_v > stat_v.memory_max)
            stat_v.memory_max = stat_v.memory_v;
        /*s*/
        
        obj_node_new->obj_v = obj_node_v1->obj_v;
        obj_node_new->dist = dist_o_1;
        
        obj_node_new->next = obj_node_v2->next;
        obj_node_v2->next = obj_node_new;
        
        olist->obj_n++;
        //-----
        
        obj_node_v1 = obj_node_v1->next;
    }
    
    //    printf("olist:\n");
    //    print_obj_set(olist, stdout);
    
    add_obj_set_entry(obj_node_v->obj_v, olist);
    if (!is_covered_obj_set(olist, q))
        goto J;
    
    remove_obj_set_entry(olist, obj_node_v->obj_v);
    
    //-----------------------------
    obj_node_v1 = olist->head->next;
    while( obj_node_v1 != NULL)
    {
        obj_v1 = obj_node_v1->obj_v;
        
        loc_v1 = get_obj_loc( obj_v1);
        dist_o_1 = calc_dist_loc( loc_v1, loc_v);
        /*
         printf("o:%d, o1:%d\n",obj_node_v->obj_v->id, obj_v1->id);
         printf("dist:%f, cur:%f\n", dist_o_1,c_cur->radius);
         */
        
        if( dist_o_1 > 2 * c_cur->radius){
            release_loc( loc_v1);
            break;
        }
        //  printf("v1:%d\n",obj_v1->id);
        
        //2
        obj_node_v2 = olist->head->next;
        while( obj_node_v2 != NULL)
        {
            obj_v2 = obj_node_v2->obj_v;
            
            loc_v2 = get_obj_loc( obj_v2);
            dist_o_2 = calc_dist_loc( loc_v2, loc_v);
            dist_1_2 = calc_dist_loc( loc_v1, loc_v2);
            release_loc( loc_v2);
            
            /*
             printf("o:%d, o1:%d, o2:%d\n",obj_node_v->obj_v->id, obj_v1->id,obj_v2->id);
             printf("dist:%f, dist2:%f, dist3:%f, cur:%f\n", dist_o_1,dist_o_2,dist_1_2,c_cur->radius);
             */
            if( dist_o_2 > dist_o_1 ){
                break;
            }
            
            // * this is different from the code in the paper
            if (  dist_1_2 > 2 * c_cur->radius){
                obj_node_v2 = obj_node_v2->next;
                continue; //* not break;
            }
            
            //   printf("\tv2:%d\n",obj_v2->id);
            
            if( obj_v1 == obj_v2){
                obj_set_v = alloc_obj_set();
                add_obj_set_entry(obj_node_v->obj_v, obj_set_v);
                add_obj_set_entry(obj_v1, obj_set_v);
                c_can = MCC(obj_set_v);
                release_obj_set(obj_set_v);
            }else{
                obj_set_v = alloc_obj_set();
                add_obj_set_entry(obj_node_v->obj_v, obj_set_v);
                add_obj_set_entry(obj_v1, obj_set_v);
                add_obj_set_entry(obj_v2, obj_set_v);
                c_can = MCC(obj_set_v);
                release_obj_set(obj_set_v);
            }
            // printf("c_can->radius:%f \t c_cur->radius:%f\n", c_can->radius, c_cur->radius);
            
            if ( c_can->radius < c_cur->radius){
                
                S = alloc_obj_set();
                psi_v = alloc_psi();
                copy_k_list(psi_v->k_head, q->psi_v->k_head);
                psi_v->key_n = q->psi_v->key_n;
                
                psi_exclusion(psi_v, obj_node_v->obj_v->fea);
                
                //for each object in O
                obj_node_temp = olist->head->next;
                //3
                while (obj_node_temp != NULL){
                    
                    if (is_enclosed(obj_node_temp->obj_v->MBR, c_can)){
                        int key_n_prev = psi_v->key_n;
                        psi_exclusion(psi_v, obj_node_temp->obj_v->fea);
                        if (key_n_prev != psi_v->key_n)//if some keyword is removed
                            add_obj_set_entry(obj_node_temp->obj_v, S);
                        
                    }
                    if(psi_v->key_n==0){
                        //all keyword are covered
                        release_disk(c_cur);
                        c_cur = c_can;
                        c_can = NULL;
                        
                        release_obj_set(S_cur);
                        add_obj_set_entry(obj_node_v->obj_v, S);
                        S_cur = S;
                        S = NULL;
                        
                        
                        break;
                    }
                    
                    obj_node_temp=obj_node_temp->next;
                } //end while 3
                release_psi(psi_v);
                if (S != NULL)
                    release_obj_set(S);
            }
            if (c_can!=NULL)
                release_disk(c_can);
            obj_node_v2 = obj_node_v2->next;
            
        }//end while 2
        
        release_loc( loc_v1);
        obj_node_v1 = obj_node_v1->next;
        
    }//end while 1
J:
    release_loc(loc_v);
    release_obj_set(olist);
    
    return;
}
//===============================================================
//===============================================================
//===============================================================
void findAppOSKEC(obj_set_t* O_t, obj_node_t* o, B_KEY_TYPE& c_cur_diam, query_t* q, B_KEY_TYPE eps, B_KEY_TYPE LB, obj_set_t*& cur_group, B_KEY_TYPE& cost_curGroup){
    
    B_KEY_TYPE UB, diam, alpha;
    bool updatedResult;
    disk_t* disk_v;
    obj_node_t* o2;
    
    UB = c_cur_diam;                     //cir dia
    
    
    //---------------------------------------------------
    //perform range query on each obj in region_u with UB radius
    
    //----
    disk_v = alloc_disk(2);
    disk_v->radius = UB;
    disk_v->center->coord[0] = o->obj_v->MBR[0].min;
    disk_v->center->coord[1] = o->obj_v->MBR[1].min;
    o->range = range_query(disk_v, q);
    //compute dist for each o' in range
    o2 = o->range->head->next;
    while (o2 != NULL){
        o2->dist = calc_dist_obj(o->obj_v, o2->obj_v);
        o2 = o2->next;
    }
    release_disk(disk_v);
    
    //--------------------------------------------------
    
    int cnt=1;
    updatedResult = circleScan (O_t, o->obj_v, UB, q, 2, cur_group,  cost_curGroup);
    if (!updatedResult) return  ;
    
    alpha = LB * eps;
    while ( (UB-LB) > alpha ){
        //      printf("UB:%f\t LB:%f\n",UB,LB);
        diam = (UB+LB)/(double)2.0;
        
        updatedResult = circleScan (O_t, o->obj_v, diam, q, 2, cur_group,  cost_curGroup);
        cnt++;
        if (updatedResult){
            // printf("oskec not null\n");
            UB = diam;
            c_cur_diam = diam;
        }else{//no result with diam found
            LB = diam;
        }
        
    }
    //  printf("#of circleScan called:%d\n",cnt);
    return ;
}




//===============================================================
//===============================================================
//===============================================================
//SKECa+

obj_set_t* SKECaplus( query_t* q, B_KEY_TYPE eps, bool callByExact, obj_set_t*& region_u, disk_t* disk_v)
{
    //obj_set_t* region_u;
    obj_node_t *o;
    disk_t* c_cur;
    obj_set_t* G_gkg, * G_skec = NULL;
    B_KEY_TYPE UB, LB, diam, alpha, cost_curGroup;
    bool foundResult, updatedResult;
    
    G_gkg = GKG(q, region_u, disk_v);
    if ( comp_diameter(G_gkg) <= dist_thr)
    {
        return G_gkg;
    }
    
    G_skec = G_gkg;
    cost_curGroup = comp_diameter(G_gkg);
    
    c_cur = MCC(G_gkg);
    
    UB = 2 * c_cur->radius;                     //cir dia
    LB = cost_curGroup / (double) 2.0;  //S dia / 2
    
    release_disk(c_cur);
    
    
    //retreive all relevant objects
    //region_u = range_query(q);
    
    //--------------------------------------------------
    alpha = LB * eps;
    while ( (UB-LB) > alpha )
    {
        diam = (UB+LB) / (double)2.0;
        foundResult = false;
        
        //  printf("LB:%f \t UB:%f \t diam:%f \t alpha:%f\n",UB,LB,diam,alpha);
        
        //---
        //even if we find a feasible set in here, the diameter > dist_thr
        if(diam > dist_thr)
        {
            release_obj_set(G_skec);
            return NULL;
        }
        //---
        
        o = region_u->head->next;
        while ( o!= NULL)
        {
            if ( diam < o->maxInvalidRange)
            {
                o=o->next;
                continue;
            }
            
            updatedResult = circleScan(region_u, o->obj_v, diam, q, 2, G_skec, cost_curGroup);
            if (updatedResult)
            {
                UB = diam;
                foundResult = true;
                //---
                if(diam <= dist_thr)
                    return G_skec;
                //---
                break;
            }else //no result with diam found
            {
                if( diam > o->maxInvalidRange)
                    o->maxInvalidRange = diam;
            }
            
            o=o->next;
        }//end inner while
        if (!foundResult)
            LB = diam;
    }//end outer while
    
    /*
     //release resources
     if (!callByExact){
     o = region_u->head->next;
     while ( o!= NULL){
     release_obj_set(o->range);
     o = o->next;
     }
     release_obj_set( region_u );
     }
     */
    
    return G_skec;
}


//s_tag=1 = exact
//s_tag=2 = approx = SKECa+
bool circleScan(obj_set_t* region_u, obj_t* o, B_KEY_TYPE diam, query_t* q, int s_tag, obj_set_t*& curGroup, B_KEY_TYPE& cost_curGroup){
    
    int top;
    obj_set_t * O_t, * list,* G = NULL;
    obj_node_t* obj_node_v;
    // disk_t* disk_v1 = NULL;
    B_KEY_TYPE angle_in, angle_out;
    bst_t* Tab;
    obj_t * obj_v, * obj_v1;
    b_heap_t* R_heap = NULL;
    B_KEY_TYPE minAngle=INFINITY;
    bool flagSetUpdated = false;
    psi_t* psi_v;
    k_node_t* k_head;
    
    Tab = const_tab( q->psi_v);
    add_Tab_key(Tab, o->fea);
    
    list = alloc_obj_set();
    G = alloc_obj_set();
    //add_obj_set_entry(o->obj_v, G);
    
    O_t = region_u;
    //---------------------------------------
    
    //psi_v is for checking whether all keyword can be covered by the objects in the sweeping area
    k_head = key_exclusion( q->psi_v->k_head, o->fea);
    psi_v = const_psi( k_head);
    
    //for each object in O_t, we have 2 circles with r = diam/2
    //only consider objects with d(obj_v,o) <= diam
    obj_node_v = O_t->head->next;
    while (obj_node_v != NULL){
        
        obj_v = obj_node_v->obj_v;
        
        if( obj_v == o || calc_dist_obj(obj_v, o)> diam)
        {
            obj_node_v = obj_node_v->next;
            continue;
        }
        
        //--
        if(psi_v->key_n > 0)
            psi_exclusion(psi_v, obj_v->fea);
        //--
        
        getAngle2(o, obj_v, diam/2, angle_in, angle_out);
        
        add_obj_set_entry_extended2(obj_v, list, angle_out, false);
        if (angle_out > angle_in)
        {
            add_obj_set_entry(obj_v, G);
            add_Tab_key(Tab, obj_v->fea);
            //we store the disk center with min angle (use if tab is full)
            if (angle_in < minAngle)
            {
                minAngle = angle_in;
                obj_v1 = obj_v;
            }
        }
        else
        {
            add_obj_set_entry_extended2(obj_v, list, angle_in , true);
        }
        
        obj_node_v = obj_node_v->next;
    }///end while
    
    if(psi_v->key_n > 0){ //obj set is NOT feasible
        // printf("psi:%d\n", psi_v->key_n);
        goto J;
    }
    
    if (is_full(Tab->root)){
        if ( s_tag==2 ){
            flagSetUpdated = true;
            add_obj_set_entry(o, G);
            release_obj_set(curGroup);
            curGroup = G;
            G = NULL;
            goto J;
        }else //exact
        {
            search(q, o, obj_v1 , G, curGroup, cost_curGroup);
            if(cost_curGroup <= dist_thr)
            {
                goto J;
            }
        }
    }
    
    //-------------------------------------------
    //Sort the objects in R by their angle.
    R_heap = heap_sort_obj_set( list);
    
    //for each tuple in list
    while( true)
    {
        top = b_h_get_top( R_heap);
        
        if( top == 0)
            goto J;
        else
            obj_v1 = R_heap->obj_arr[ top].obj_v;
        
        if (  R_heap->obj_arr[ top].type == true)//in
        {
            add_Tab_key(Tab, obj_v1->fea);
            //in each possible circle, we enforce o, o1 in the selectedSet to avoid duplicate enum
            //thus o1 is added to G after search
            if (is_full(Tab->root)){
                if ( s_tag == 1 ){//exact
                    search(q, o, obj_v1 , G, curGroup, cost_curGroup);
                    //---
                    if(cost_curGroup <= dist_thr)
                    {
                        goto J;
                    }
                    //---
                }else {//appro
                    add_obj_set_entry(obj_v1, G);
                    break;
                }
            }
            add_obj_set_entry(obj_v1, G);
            
        }else //out
        {
            remove_obj_set_entry(G, obj_v1);
            remove_Tab_key(Tab, obj_v1->fea);
        }
        
    }//end while
    
    //Tab is full
    if (s_tag==2){
        flagSetUpdated = true;
        add_obj_set_entry(o, G);
        release_obj_set(curGroup);
        curGroup = G;
        G = NULL;
    }
    
J:
    
    if( R_heap!= NULL)
        release_b_heap(R_heap);
    
    release_obj_set(G);
    bst_release(Tab);
    release_obj_set(list);
    release_psi(psi_v);
    
    return flagSetUpdated;
    
}

void getAngle2( obj_t* o, obj_t* o_j, B_KEY_TYPE r, B_KEY_TYPE& o_in, B_KEY_TYPE& o_out)
{
    
    B_KEY_TYPE PI = 3.14159265;
    B_KEY_TYPE x1,y1,x2,y2,q;
    B_KEY_TYPE a,b;
    
    x1=o->MBR[0].min;
    y1=o->MBR[1].min;
    x2=o_j->MBR[0].min;
    y2=o_j->MBR[1].min;
    
    q = sqrt(pow(x1-x2,(double)2.0) + pow(y1-y2,(double)2.0))/2;
    b = acos(q/r);
    a = atan2(y2-y1, x2-x1);
    
    // printf("q:%f\t r:%f\t b:%f\t a:%f\n",q,r,b,a);
    
    o_in =  a + b;
    o_out = a - b;
    
    if (o_in < 0 ) o_in+=2*PI;
    if (o_out < 0 ) o_out+=2*PI;
    
    //printf("(%f,%f),",x2,y2);
    //  printf("diam:%f \t o_in:%f \t o_out:%f \n",r,o_in,o_out);
    
}


void getAngle( obj_node_t* o, obj_node_t* o_j, B_KEY_TYPE r, B_KEY_TYPE& o_in, B_KEY_TYPE& o_out){
    
    //c1 c2 are the circle
    //x1 y1 is o
    //x2 y2 is o_j
    //x3 y3 is mid pt of o, o_j
    // q is d(o,o_j)
    B_KEY_TYPE PI = 3.14159265;
    B_KEY_TYPE c1x,c1y,c1a,c2x,c2y,c2a,x1,y1,x2,y2,x3,y3,q,p,temp;
    
    x1=o->obj_v->MBR[0].min;
    y1=o->obj_v->MBR[1].min;
    x2=o_j->obj_v->MBR[0].min;
    y2=o_j->obj_v->MBR[1].min;
    
    x3 = (x1+x2)/2;
    y3 = (y1+y2)/2;
    q = sqrt(pow(x1-x2,(double)2.0) + pow(y1-y2,(double)2.0));
    p = sqrt(pow( r,2)- pow(q/2,2))/q;
    //-----
    c1x = x3 + p * (y1-y2);
    c1y = y3 + p * (x2-x1);
    temp =  atan2(c1y-y1, c1x-x1);
    if (temp <0)
        temp += 2* PI;
    c1a = temp ;
    
    c2x = x3 - p * (y1-y2);
    c2y = y3 - p * (x2-x1);
    temp = atan2(c2y-y1, c2x-x1);
    if (temp < 0)
        temp += 2* PI;
    
    c2a = temp ;
    //------
    /*
     printf("x1:%0.3lf\t y1:%0.3lf \t\n",x1,y1);
     printf("x2:%0.3lf\t y2:%0.3lf \t\n",x2,y2);
     printf("x3:%0.3lf\t y3:%0.3lf \t\n",x3,y3);
     
     printf("c1x:%0.3lf\tc1y:%0.3lf \t c1a:%f\n",c1x,c1y, c1a);
     printf("c2x:%0.3lf\tc2y:%0.3lf \t c2a:%f\n",c2x,c2y, c2a);
     */
    //usually angle smaller = out
    if ( c1a < c2a){
        o_out = c1a;
        o_in = c2a;
        
    }else{
        o_out = c2a;
        o_in = c1a;
    }
    
    //special case: if they pass 0
    //swap in and out
    if( o_in - o_out > PI){
        B_KEY_TYPE temp = o_in;
        o_in = o_out;
        o_out = temp;
    }
    //    printf("x_in:%0.3lf\ty_in:%0.3lf\n",x_in,y_in);
    //    printf("x_o:%0.3lf\ty_o:%0.3lf\n",x_out,y_out);
    
    //check
    /*
     loc_t* loc_v1 = alloc_loc(2);
     loc_v1->coord[0] = c1x;
     loc_v1->coord[1] = c1y;
     loc_t* loc_v2 = alloc_loc(2);
     loc_v2->coord[0] = c2x;
     loc_v2->coord[1] = c2y;
     
     loc_t* loc_v = get_obj_loc(o->obj_v);
     loc_t* loc_j = get_obj_loc(o_j->obj_v);
     
     B_KEY_TYPE d1 = calc_dist_loc(loc_v1,loc_v );
     B_KEY_TYPE d2 = calc_dist_loc(loc_v1,loc_j );
     
     B_KEY_TYPE d3 = calc_dist_loc(loc_v2,loc_v );
     B_KEY_TYPE d4 = calc_dist_loc(loc_v2,loc_j );
     
     printf("d1:%f \t d2:%f \t r:%f\n",d1,d2, r);
     printf("d3:%f \t d4:%f \t r:%f\n",d3,d4, r);
     printf("a1:%f \t a2:%f \n",c1a,c2a);
     */
}

/*
 *	Construct the IF on a set of objects @obj_set_v for the keywords in @psi_v.
 *
 *	1. The IF structure is indexed by a binary search tree.
 *	2. No ordering is imposed in IF.
 *  3. this is not inverted list, just store keyword and freq
 */
bst_t* const_tab( psi_t* psi_v)
{
    int i;
    bst_t* IF_v;
    k_node_t* k_node_v;
    bst_node_t* bst_node_v;
    
    //IF_v = alloc_IF( psi_v->key_n);
    IF_v = bst_ini( );
    
    k_node_v = psi_v->k_head->next;
    ///for each keyword in psi_v
    for( i=0; i<psi_v->key_n; i++)
    {
        //IF_v->entry_v[ i].key = k_node_v->key;
        bst_node_v = ( bst_node_t*)malloc( sizeof( bst_node_t));
        memset( bst_node_v, 0, sizeof( bst_node_t));
        
        /*s*/
        stat_v.memory_v += sizeof( bst_node_t);
        /*s*/
        //--
        
        bst_node_v->freq = 0;
        //--
        bst_node_v->key = k_node_v->key;
        bst_node_v->p_list_obj = NULL;
        
        bst_insert( IF_v, bst_node_v);
        
        k_node_v = k_node_v->next;
    }
    
    /*s*/
    if( stat_v.memory_v > stat_v.memory_max)
        stat_v.memory_max = stat_v.memory_v;
    /*s*/
    
    return IF_v;
}



bool is_full(bst_node_t* bst_node_v){
    
    if (bst_node_v == NULL)
        return true;
    
    if (bst_node_v->freq > 0)
        return is_full(bst_node_v->left) && is_full(bst_node_v->right);
    else
        return false;
    
    
}

/*
 *	Add an object entry @obj_v to @obj_set_v.
 */
void add_obj_set_entry_extended2( obj_t* obj_v, obj_set_t* obj_set_v, B_KEY_TYPE angle, bool type)
{
    obj_node_t* obj_node_v;
    
    obj_node_v = ( obj_node_t*)malloc( sizeof( obj_node_t));
    memset( obj_node_v, 0, sizeof( obj_node_t));
    
    obj_node_v->obj_v = obj_v;
    //--
    obj_node_v->angle = angle;
    obj_node_v->type  = type;
    //--
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
 *	Sort the objects in @obj_set_v by their angle.
 *
 *	Method: heap-sort.
 *	we multiply -1 because min heap is used, and we want larger first
 */
b_heap_t* heap_sort_obj_set( obj_set_t* obj_set_v)
{
    int cur;
    b_heap_t* b_h;
    obj_node_t* obj_node_v;
    
    b_h = alloc_b_heap( obj_set_v->obj_n + 1);
    
    cur = 1;
    obj_node_v = obj_set_v->head->next;
    while( obj_node_v != NULL)
    {
        b_h->obj_arr[ cur].key = -1 * obj_node_v->angle;
        b_h->obj_arr[ cur].obj_v = obj_node_v->obj_v;
        b_h->obj_arr[ cur].type = obj_node_v->type;
        
        b_h_insert( b_h, cur);
        
        cur ++;
        
        obj_node_v = obj_node_v->next;
    }
    
    return b_h;
}


void add_Tab_key(bst_t* Tab, k_node_t* k_head){
    
    k_node_t* k_node_v ;
    bst_node_t* bst_node_v;
    
    k_node_v = k_head->next;
    while ( k_node_v !=NULL){
        bst_node_v =  bst_search(Tab, k_node_v->key);
        if (bst_node_v!=NULL) //if the keyword is one of the query keyword
            bst_node_v->freq++;
        
        k_node_v=k_node_v->next;
    }
    
}


void add_Tab_key(bst_t* Tab, FEA_TYPE fea)
{
    
    //    k_node_t* k_node_v ;
    bst_node_t* bst_node_v;
    
    //    k_node_v = k_head->next;
    //    while ( k_node_v !=NULL){
    bst_node_v =  bst_search(Tab, fea);
    if (bst_node_v!= NULL) //if the keyword is one of the query keyword
        bst_node_v->freq++;
    
    //       k_node_v=k_node_v->next;
    //   }
    
}

void remove_Tab_key(bst_t* Tab, k_node_t* k_head){
    k_node_t* k_node_v ;
    bst_node_t* bst_node_v;
    
    k_node_v = k_head->next;
    while ( k_node_v != NULL){
        bst_node_v =  bst_search(Tab, k_node_v->key);
        if (bst_node_v != NULL) //if the keyword is one of the query keyword
            bst_node_v->freq--;
        
        k_node_v=k_node_v->next;
    }
    
}

void remove_Tab_key(bst_t* Tab, FEA_TYPE fea)
{
    //   k_node_t* k_node_v ;
    bst_node_t* bst_node_v;
    
    //    k_node_v = k_head->next;
    //    while ( k_node_v != NULL){
    bst_node_v =  bst_search(Tab, fea);
    if (bst_node_v != NULL) //if the keyword is one of the query keyword
        bst_node_v->freq--;
    
    //        k_node_v=k_node_v->next;
    //    }
}

void print_Tab(bst_t* Tab)
{
    bst_node_t* bst_node_v;
    
    bst_node_v = bst_get_min(Tab->root);
    while ( bst_node_v !=NULL){
        printf("key:%.0f\t freq:%d \n", bst_node_v->key, bst_node_v->freq);
        bst_node_v= bst_successor(bst_node_v);
    }
}

/*
 *	Check wheher the objects in @obj_set_v contain an obj obj_v.
 */
bool is_contained_obj_set( obj_set_t* obj_set_v, obj_t* obj_v)
{
    obj_node_t* obj_node_iter;
    
    obj_node_iter = obj_set_v->head->next;
    while( obj_node_iter != NULL)
    {
        if( obj_node_iter->obj_v==obj_v)
            return true;
        obj_node_iter = obj_node_iter->next;
    }
    return false;
}


//==================================================

void search( query_t* q, obj_t* o, obj_t* o1, obj_set_t* G, obj_set_t*& curGroup, B_KEY_TYPE& cost_curGroup){
    
    B_KEY_TYPE cost_selectedSet;
    obj_set_t* selectedSet;
    psi_t* psi_v;
    
    selectedSet = alloc_obj_set();
    add_obj_set_entry(o, selectedSet);
    if(o1!=NULL)
        add_obj_set_entry(o1, selectedSet);
    
    cost_selectedSet = comp_cost(cost_tag, selectedSet, q);
    psi_v = uncover_keyword(q, selectedSet);
    
    search_sub(q, selectedSet, G, 0, cost_selectedSet, psi_v, curGroup, cost_curGroup);
    
    release_psi(psi_v);
    release_obj_set(selectedSet);
}



void search_sub( query_t* q, obj_set_t* selectedSet, obj_set_t* candidateSet, int maxId, B_KEY_TYPE& cost_selectedSet, psi_t* psi_uncover_selectedSet, obj_set_t*& curGroup, B_KEY_TYPE& cost_curGroup){
    //  printf("\t stat_v.memory_v:%f\n",stat_v.memory_v/(1024*1024));
    
    obj_set_t* nextSet;
    psi_t * psi_v, * psi_next;
    obj_node_t* o_c, * o_n, * o_s;
    B_KEY_TYPE selectedDiam;
    
    //    printf("psi_uncover:%d\n",psi_uncover_selectedSet->key_n);
    
    if (psi_uncover_selectedSet->key_n==0){
        //    if (is_covered_obj_set(selectedSet, q)){
        if(cost_selectedSet <= cost_curGroup){
            release_obj_set(curGroup);
            curGroup = copy_obj_set( selectedSet );
            cost_curGroup = cost_selectedSet;
        }
        return;
    }
    if (cost_selectedSet > cost_curGroup)
        return;
    
    nextSet = alloc_obj_set();
    
    //    psi_v = uncover_keyword(q, selectedSet);
    psi_v = alloc_psi();
    copy_k_list(psi_v->k_head, psi_uncover_selectedSet->k_head);
    psi_v->key_n = psi_uncover_selectedSet->key_n;
    
    //for each obj in cand set
    o_c = candidateSet->head->next;
    while( o_c != NULL) {
        
        if(o_c->obj_v->id < maxId){
            o_c=o_c->next;
            continue;
        }
        //---1---
        selectedDiam = cost_selectedSet;
        o_s = selectedSet->head->next;
        while( o_s != NULL){
            selectedDiam = fmax(selectedDiam, calc_dist_obj(o_s->obj_v, o_c->obj_v));
            o_s = o_s->next;
        }
        if ( selectedDiam > cost_curGroup){
            o_c = o_c->next;
            continue;
        }
        //---2----
        if(!is_relevant_obj(o_c->obj_v, psi_v)){
            //        if (number_intersection(psi_v->k_head, o_c->obj_v->k_head)==0){
            o_c=o_c->next;
            continue;
        }
        //---------
        
        add_obj_set_entry(o_c->obj_v, nextSet);
        nextSet->head->next->dist = selectedDiam;
        
        //** remove keywords appear in o_c
        psi_exclusion(psi_v, o_c->obj_v);
        o_c = o_c->next;
        
    }//end while
    
    //-----3-----
    
    if (psi_v->key_n != 0){
        release_psi(psi_v);
        goto E;
    }
    release_psi(psi_v);
    //-------------
    
    // newcandSet = copy_obj_set(nextSet);
    
    o_n = nextSet->head->next;
    while (o_n != NULL){
        obj_t* obj_v = o_n->obj_v;
        
        add_obj_set_entry(obj_v, selectedSet);
        //--
        psi_next = alloc_psi();
        copy_k_list(psi_next->k_head, psi_uncover_selectedSet->k_head);
        psi_next->key_n = psi_uncover_selectedSet->key_n;
        psi_exclusion(psi_next, obj_v->fea);
        //--
        if (selectedSet->obj_n <= q->psi_v->key_n)
        {
            //not necessary to remove o_n from nextSet
            //infinite loop may occur !
            B_KEY_TYPE cost_newselSet = fmax(cost_selectedSet, o_n->dist);
            search_sub(q, selectedSet, nextSet, obj_v->id, cost_newselSet, psi_next, curGroup, cost_curGroup);
            //---
            if(cost_curGroup <= dist_thr)
            {
                release_psi(psi_next);
                release_obj_set(nextSet);
                return;
            }
            //---
        }
        remove_obj_set_entry(selectedSet);
        release_psi(psi_next);
        
        o_n=o_n->next;
    }
    
E:
    release_obj_set(nextSet);
    
    return;
}

/*
 *	Check whether the keywords in the query @q are covered by a set of objs in @obj_set_v.
 *   only check obj with dist  < d
 */
bool is_covered_obj_set( obj_set_t* obj_set_v, query_t* q, B_KEY_TYPE d)
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
            // printf("dist:%f\n",obj_node_iter->dist);
            if( obj_node_iter->dist <= d && has_key_obj( obj_node_iter->obj_v, key))
                break;
            
            obj_node_iter = obj_node_iter->next;
        }
        
        if( obj_node_iter == NULL)
            return false;
        
        k_node_iter = k_node_iter->next;
    }
    
    return true;
}

//======================================================================
//======================================================================
//======================================================================

/*
 *	 range query.
 *
 *	DFS: recursive implementation.
 * modified from range_query and range_query_sub, this version does not need the disk
 * all objs are retrived
 */
obj_set_t* range_query3()
{
    obj_set_t* obj_set_v;
    
    obj_set_v = alloc_obj_set( );
    retrieve_sub_tree3(  IRTree_v.root, obj_set_v);
    
    return obj_set_v;
}

/*
 *	Retrieve all the objects located at the sub-tree rooted at @node_v.
 *	The retrieved objects are stored in obj_set_v.
 */
void retrieve_sub_tree3( node_t* node_v, obj_set_t* &obj_set_v)
{
    int i;
    //	BIT_TYPE p_list;
    
    if( node_v->level == 0)
    {
        //node_v is a leaf-node.
        //Retrieve all its objects.
        for( i=0; i<node_v->num; i++)
        {
            //if( is_relevant_obj( ( obj_t*)( node_v->child[ i]), q))
            add_obj_set_entry( ( obj_t*)( node_v->child[ i]), obj_set_v);
        }
    }
    else
    {
        //node_v is an inner-node.
        //Invoke the function recursively.
        //	p_list = is_relevant_node( node_v, q);
        for( i=0; i<node_v->num; i++)
        {
            //	if( get_k_bit( p_list, i))
            retrieve_sub_tree3( ( node_t*)( node_v->child[ i]), obj_set_v);
        }
    }
}


/*
 *	The implementation of gen_query_set is problematic.
 *
 */
/*
 query_t** gen_query_set3( int query_n, int key_n, range* MBR, data_t* data_v, int low, int high, B_KEY_TYPE maxDia)
 {
	int i, j, k, cnt, k_max, k_num, rand_tmp, valid_k_sta, valid_k_end, valid_k_num;
	float p_low, p_high;
	int* freq_v, *rand_v;
	loc_t* loc_v;
	psi_t* psi_v;
	query_t** q_set;
	key_freq_pair_t* key_freq_pair_v;
 obj_set_t* obj_set_v;
 
	if( low > high || low < 0 || low > 100 || high < 0 || high > 100)
 return NULL;
 
 
 
	//Collect the "keyword range" information.
	k_max = get_max_key( data_v);
	
	
	//Collect the "frequency" information.
	freq_v = ( int*)malloc( ( k_max + 1) * sizeof( int));
	memset( freq_v, 0, ( k_max + 1) * sizeof( int));
 
	collect_key_freq( data_v, freq_v, k_max + 1);
 printf( "Collected keyword freq\n");
	
	//Restore the frequency information.
	//1. Collect the number of keys.
	k_num = get_key_num( freq_v, k_max + 1);
 
	if( k_num < key_n)
	{
 printf( "Insufficient keywords!\n");
 exit( 0);
	}
	printf( "Sufficient keywords!\n");
	
	//2. Restore.
	key_freq_pair_v = ( key_freq_pair_t*)malloc( k_num * sizeof( key_freq_pair_t));
	memset( key_freq_pair_v, 0, k_num * sizeof( key_freq_pair_t));
 
	cnt = 0;
	for( i=0; i<=k_max; i++)
	{
 if( freq_v[ i] != 0)
 {
 key_freq_pair_v[ cnt].key = i;
 key_freq_pair_v[ cnt++].freq = freq_v[ i];
 }
	}
 
	//3. Sort the keywords id based on the freq.
	qsort( key_freq_pair_v, k_num, sizeof( key_freq_pair_t), compare_key_freq);
	
	//Locate the keywords within the percentile range.
	p_low = ( 100 - high) / 100.0;
	p_high = ( 100 - low) / 100.0;
	
	valid_k_sta = 0;
	while( ( valid_k_sta + 1) / ( float)k_num < p_low && valid_k_sta < k_num)
 valid_k_sta ++;
 
	if( valid_k_sta == k_num)
	{
 printf( "Gen_query2 bug.\n");
 return NULL;
	}
 
 
	valid_k_end = valid_k_sta;
	while( ( valid_k_end + 1) / ( float)k_num < p_high && valid_k_end < k_num)
 valid_k_end ++;
 
	valid_k_num = valid_k_end - valid_k_sta + 1;
	if( valid_k_num < key_n)
	{
 printf( "Insufficient valid keywords!\n");
 exit( 0);
	}
 
	//Generate the queries.
	q_set = ( query_t**)malloc( sizeof( query_t*) * query_n);
	memset( q_set, 0, sizeof( query_t*) * query_n);
 printf( "Gen_query2 no bug.\n");
 
 i=0;
 while (i <query_n)
 //	for( i=0; i<query_n; i++)
	{
 
 ///printf( "Gen_query %d: ",i);
 
 ///allocate memory
 query_t* q = alloc_query();
 
 
 //q_set[ i] = alloc_query( );
 psi_v = alloc_psi( );
 
 rand_v = ( int*)malloc( key_n * sizeof( int));
 memset( rand_v, 0, sizeof( int));
 
 /// generate key_n query keywords and store in rand_v
 for( k=0; k<key_n; k++)
 {
 /// printf( "%d ",k);
 
 
 rand_tmp = rand_i( 0, valid_k_num-1);
 while( is_old( rand_v, k, rand_tmp))
 rand_tmp = rand_i( 0, valid_k_num-1);
 
 rand_v[ k] = rand_tmp;
 
 add_psi_entry( psi_v, key_freq_pair_v[ valid_k_sta + rand_v[ k]].key);
 }
 // printf( "\n");
 
 q->psi_v = psi_v;
 //q_set[ i]->psi_v = psi_v;
 
 
 //query location.
 loc_v = alloc_loc( data_v->dim);
 for( j=0; j<data_v->dim; j++)
 loc_v->coord[ j] = rand_f( MBR[ j].min, MBR[ j].max);
 
 q->loc_v = loc_v;
 //q_set[ i]->loc_v = loc_v;
 
 
 //--
 obj_set_t* obj_set_v =  Gen(q, 1, 6, NULL);
 B_KEY_TYPE cost = comp_diameter(obj_set_v);
 release_obj_set(obj_set_v);
 
 if (cost > maxDia){//not OK
 //printf("cost:%f \t maxDia:%f\n",cost, maxDia);
 
 release_query(q);
 continue;
 }
 
 q_set[i] = q;
 //  printf("i:%d\n",i);
 i++;
 //--
 free( rand_v);
	}
 
	free( freq_v);
	free( key_freq_pair_v);
 
	return q_set;
 }
 */


/*
 *	The implementation of the "Cao-Appro1" algrithm.
 *
 *	return NULL if no feasible set is possible,
 *	i.e., the keywords in @q could not be covered,
 *	otherwise, return the approximate solution.
 */
obj_set_t* Cao_Appro1( query_t* q)
{
    obj_set_t* S;
    obj_t* obj_v;
    k_node_t* k_node_v;
    
    //printf( "Cao-Appro1:\n");
    
    S = alloc_obj_set( );
    
    k_node_v = q->psi_v->k_head->next;
    while( k_node_v != NULL)
    {
        obj_v = const_NN_key( q->loc_v, k_node_v->key, NULL);
        
        if( obj_v == NULL)
        {
            printf( "No solution exists!\n");
            release_obj_set( S);
            return NULL;
        }
        
        add_obj_set_entry( obj_v, S);
        
        k_node_v = k_node_v->next;
    }
    
    return S;
}



/*
 *	The implementation of the "Cao-Appro1" algrithm.
 *
 *	return NULL if no feasible set is possible,
 *	i.e., the keywords in @q could not be covered,
 *	otherwise, return the approximate solution.
 *
 *  this version only search the objects in the disk @disk_v (but not the whole space)
 */
obj_set_t* Cao_Appro1( query_t* q, disk_t* disk_v)
{
    obj_set_t* S;
    obj_t* obj_v;
    k_node_t* k_node_v;
    
    //printf( "Cao-Appro1:\n");
    
    S = alloc_obj_set( );
    
    k_node_v = q->psi_v->k_head->next;
    while( k_node_v != NULL)
    {
        obj_v = const_NN_key( q->loc_v, k_node_v->key, disk_v);
        
        if( obj_v == NULL)
        {
            printf( "No solution exists!\n");
            release_obj_set( S);
            return NULL;
        }
        
        add_obj_set_entry( obj_v, S);
        
        k_node_v = k_node_v->next;
    }
    return S;
    
}



/*
 *	Remove the specified entry from the list of objects @obj_set_v.
 *  Modified from remove_obj_set_entry
 */
void remove_obj_set_entry( obj_set_t* obj_set_v, obj_t* obj_v )
{
    obj_node_t* obj_node_prev,* obj_node_temp;
    
    obj_node_prev = obj_set_v->head;
    // printf("finding object to remove\n");
    
    while(obj_node_prev->next != NULL){
        if(obj_node_prev->next->obj_v==obj_v){
            obj_node_temp=obj_node_prev->next;
            obj_node_prev->next = obj_node_prev->next->next;
            obj_set_v->obj_n --;
            free(obj_node_temp);
            //  printf("an object removed\n");
            /*s*/
            stat_v.memory_v -= sizeof( obj_node_t);
            /*s*/
            
            return ;
        }
        
        obj_node_prev=obj_node_prev->next;
    }
    
}



//remove the keywords in the object in the object set obj_set_v
psi_t*  uncover_keyword(query_t* q, obj_set_t* obj_set_v){
    
    psi_t* psi_v;
    obj_node_t* obj_node_v;
    k_node_t* k_node_v,* k_node_v2;
    bool flag;
    
    psi_v = alloc_psi();
    
    //for each query keyword
    k_node_v = q->psi_v->k_head->next;
    while (k_node_v!= NULL){
        
        flag = false;
        
        //for each object
        obj_node_v = obj_set_v->head->next;
        while (obj_node_v != NULL) {
            
            //for each keyword in the object
            //            k_node_v2 = obj_node_v->obj_v->k_head->next;
            //            while (k_node_v2 != NULL){
            //
            //                if (k_node_v->key == k_node_v2->key){
            //                    flag = true;
            //                    break;
            //                }
            //                k_node_v2 = k_node_v2->next;
            //
            //            }
            //            if (flag)
            if (obj_node_v->obj_v->fea == k_node_v->key)
            {
                flag = true;
                break;
            }
            obj_node_v=obj_node_v->next;
        }
        //if no object contain this keyword
        if (!flag){
            add_psi_entry(psi_v, k_node_v->key);
        }
        
        k_node_v=k_node_v->next;
    }
    
    return psi_v;
}



/*
 *	Check whether an object @obj_v is "relevant" to the query @q.
 *	That is, whether @obj_v contains a keyword in the query @q.
 *  modified from is_relevant_obj( obj_v, q)
 */
int is_relevant_obj( obj_t* obj_v, psi_t* psi_v)
{
    int cnt=0;
    k_node_t* k_node_v;
    
    k_node_v = psi_v->k_head->next;
    while( k_node_v != NULL)
    {
        if( has_key_obj( obj_v, k_node_v->key))
            cnt++;
        
        k_node_v = k_node_v->next;
    }
    
    return cnt;
}




/*
 *	Exclude the keywords that occur in @obj_v from @psi_v1.
 *
 *	
 */
void psi_exclusion( psi_t* psi_v1, obj_t* obj_v)
{
    
    psi_exclusion(psi_v1, obj_v->fea);
    return;
}
