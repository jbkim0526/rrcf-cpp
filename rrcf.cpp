#include "rrcf.h"
#include <stdio.h>
#include <assert.h>

Node::Node(int t){
    type = t;
    if(t == LEAF) n = 1;
}

Leaf::Leaf(int* point, int size, int index,int depth){
    s = size;
    x = point;
    i = index;
    d = depth;
    b = point;
}
Branch::Branch(int _q, float _p, int** _b){
    q = _q;
    p = _p;
    b = _b;
}

RCTree::RCTree(int n){
    root = NULL;
    ndim = n;
    
}

RCTree::~RCTree(){
	_remove_all(root);
}

RCF::RCF(int num, int size,int v_size){
    num_trees = num;
    tree_size = size;
    vect_size = v_size;

    rct = (RCTree**) malloc(sizeof(RCTree*)*num_trees);
    
    for(int i = 0 ; i< num_trees ; i++){
        rct[i] = new RCTree(v_size);
    }
}

int RCTree::size(){
    return leaves.size();
}

Node* RCTree::forget_point(int index){

    Node* node = leaves[index];
    
    if(node->n>1){
        _update_leaf_count_upwards(node,-1);
        leaves.erase(index);
        return node;
    }

    //weird cases
    // if node is root
    if(node == root){
        root = NULL;
        ndim = 0;
        leaves.erase(index);
        return node;
    }

    // else find parent
    Node* parent = node->u;
    // find sibling
    Node* sibling;
    if(node == parent->l) sibling = parent->r;
    else sibling = parent->l;

    // parenet is root
    if(parent == root){
        //delete parent 
        free(parent->branch->b[0]);
        free(parent->branch->b[1]);
        free(parent->branch->b);
        delete parent->branch; 
        delete parent; 

        //set sibling as new root
        sibling->u = NULL;
        root = sibling;
        if(sibling->type == LEAF) sibling->leaf->d = 0;
        else map_leaves_decrement_depth(sibling);
        leaves.erase(index);
        return node;
    }

    //find grandparent
    Node* grandparent = parent->u;

    //set parent of sibling to grand
    sibling->u = grandparent;
    if(parent == grandparent->l) grandparent->l = sibling;
    else grandparent->r = sibling;

    //update depths

    map_leaves_decrement_depth(sibling);
    _update_leaf_count_upwards(grandparent,-1);
    int* point = node->leaf->x;
    
    _relax_bbox_upwards(grandparent,point);

    free( parent->branch->b[0]);
    free(parent->branch->b[1]);
    free(parent->branch->b);
    delete parent->branch;
    delete parent;

    leaves.erase(index);
    return node;

}

Leaf* RCTree::insert_point( int* point , int index ){
    // if tree is empty insert new.
    if(root == NULL){
        Leaf* leaf = new Leaf(point,ndim,index,0);
        Node* node = new Node(LEAF);
        node->leaf = leaf;
        root = node;
        leaves[index] = node;
        //leaves.insert(make_pair(index,leaf));
        return leaf;
    }
    // Point must be same dimension as existing points in tree.
    //assert(pdim == ndim); 
    //points must have unique key
    assert(leaves.find(index) == leaves.end());
   // if points are duplicate increment
    Node* duplicate = find_duplicate(point); // we know duplicate node is LEAF

    if(duplicate){
        _update_leaf_count_upwards(duplicate,1);
        leaves[index] = duplicate;
        //leaves.insert(make_pair(index,duplicate->leaf));
	free(point);
        return duplicate->leaf;
    }
    // if tree has points and not duplicate continue main algorithm.
    Node* node = root;
    Node* parent;
    int maxdepth = _get_maxdepth(); 
    int depth = 0;
    Node* new_node = NULL;
    Node* new_b_node = NULL;
    char side;
    for( int i = 0 ; i < maxdepth + 1 ; i++){

        // When node is LEAF
        if(node->type == LEAF){  
            Leaf* leaf = node->leaf;
            parent = node->u;

            int* bbox = leaf->b;
            tuple<int,float> tup = _insert_point_cut(point,bbox);
            int cut_dimension = get<0>(tup);
            float cut = get<1>(tup);
            if(cut <= bbox[cut_dimension]){
                Leaf* new_leaf = new Leaf(point,ndim,index,depth);
                new_node = new Node(LEAF);
                new_node->leaf = new_leaf;

                Branch* new_branch = new Branch(cut_dimension,cut);
                new_b_node = new Node(BRANCH);
                new_b_node->branch = new_branch;
                new_b_node->l = new_node;
                new_b_node->r = node;
                new_b_node->n =new_node->n + node->n;
                break;

            }
            else if(cut > bbox[cut_dimension]){
                Leaf* new_leaf = new Leaf(point,ndim,index,depth);
                new_node = new Node(LEAF);
                new_node->leaf = new_leaf;

                Branch* new_branch = new Branch(cut_dimension,cut);
                new_b_node = new Node(BRANCH);
                new_b_node->branch = new_branch;
                new_b_node->l = node;
                new_b_node->r = new_node;
                new_b_node->n =new_node->n + node->n;
                break;
            }
        }

        // When node is BRANCH
        else{
            Branch* branch = node->branch;
            parent = node->u;
            int** bbox = branch->b;

            tuple<int,float> tup = _insert_point_cut(point,bbox);
            int cut_dimension = get<0>(tup);
            float cut = get<1>(tup);
            if(cut <= bbox[0][cut_dimension]){
                Leaf* new_leaf = new Leaf(point,ndim, index,depth);
                new_node = new Node(LEAF);
                new_node->leaf = new_leaf;

                Branch* new_branch = new Branch(cut_dimension,cut);
                new_b_node = new Node(BRANCH);
                new_b_node->branch = new_branch;
                new_b_node->l = new_node;
                new_b_node->r = node;
                new_b_node->n = new_node->n + node->n;
                break;

            }
            else if(cut > bbox[1][cut_dimension]){
                Leaf* new_leaf = new Leaf(point,ndim, index,depth);
                new_node = new Node(LEAF);
                new_node->leaf = new_leaf; 

                Branch* new_branch = new Branch(cut_dimension,cut);
                new_b_node = new Node(BRANCH);
                new_b_node->branch = new_branch;
                new_b_node->l = node;
                new_b_node->r = new_node;
                new_b_node->n = new_node->n + node->n;
                break;
            }
            else{
                depth += 1;
                if(point[branch->q] <= branch->p){
                    parent = node;
                    node = node->l;
                    side = 'l';
                }
                else{
                    parent = node;
                    node = node->r;
                    side = 'r';
                }
            }

        }
    
    }

    //assert The cut was found -> new leaf and old branch 
    assert(new_b_node != NULL); 

    node->u = new_b_node;
    new_node->u = new_b_node;
    new_b_node->u = parent;

    // if parent is not null(new_b_node is not root) so there exists above node
    if(parent != NULL){
        if(side =='r') parent->r = new_b_node;
        else parent->l = new_b_node;
    }
    else root = new_b_node;


    int inc = 1;
    //since new branch is added increment the below leaf's depth
    map_leaves_increment_depth(new_b_node); 

    //since new leaf is added increment the above leaf count
    _update_leaf_count_upwards(parent,1); //why parent? -> branch's count was already calculated
    _tighten_bbox_upwards(new_b_node);
 
    leaves[index] = new_node;
    //leaves.insert(make_pair(index,new_node->leaf));
    
    return new_node->leaf;
}

float RCTree::codisp(int index){

    //assert that index is inside map
    assert(leaves.find(index) != leaves.end()); 
    Node* node = leaves[index];
    //if leaf is root -> there is only one leaf
    if(node == root) return 0;

    float result = 0;
    int d = node->leaf->d;

    for( int i = 0; i< d; i++){
        Node* parent = node->u;
        Node* sibling;
        if(parent == NULL) break;
        
        if(node == parent->l) sibling = parent->r;
        else sibling = parent->l;

        int num_deleted = node->n;
        int displacement = sibling->n;

        float temp_result = displacement/num_deleted;
        if(temp_result > result) result = temp_result;

        node = parent;
    }

    return result;
}

void RCTree::map_leaves_increment_depth(Node* node){
    if(node->type ==BRANCH){
        if(node->l) map_leaves_increment_depth(node->l);
        if(node->r) map_leaves_increment_depth(node->r);
    }
    else _increment_depth(node);
}

void RCTree::map_leaves_decrement_depth(Node* node){
    if(node->type ==BRANCH){
        if(node->l) map_leaves_decrement_depth(node->l);
        if(node->r) map_leaves_decrement_depth(node->r);
    }
    else _decrement_depth(node);
}


void RCTree::_tighten_bbox_upwards(Node* node){ // node is Branch type
    
    int** bbox = _lr_branch_bbox(node);

    //change to new bbox.
    if(node->branch->b != NULL){
        free( node->branch->b[0]);
        free( node->branch->b[1]);
        free( node->branch->b);
    }
    node->branch->b = bbox;
    node = node->u;

    while(node){
        int flag = 0;
        int** temp_bbox= node->branch->b;

        for(int i = 0 ; i < ndim; i++){
            if( bbox[0][i] < temp_bbox[0][i]){
                temp_bbox[0][i] = bbox[0][i];
                flag = 1;
            }
            if(temp_bbox[1][i] < bbox[1][i]){
                temp_bbox[1][i] = bbox[1][i];
                flag = 1;
            }
        }
        if(flag == 1) break; 
        node = node->u;
    }

}

//allocate memory and create bbox of child nodes
int** RCTree::_lr_branch_bbox(Node* node){ // node is Branch type
    int** bbox;
    bbox = (int**) malloc(sizeof(int*)*2);
    bbox[0] = (int*) malloc(sizeof(int)*ndim);
    bbox[1] = (int*) malloc(sizeof(int)*ndim);

    Node* r_child = node->r;
    Node* l_child = node->l;
    int l_child_type = l_child->type;
    int r_child_type = r_child->type;

    if(l_child_type  == LEAF && r_child_type == LEAF){

        int* l_v = l_child->leaf->b;
        int* r_v = r_child->leaf->b;

        for(int i = 0 ; i < ndim ; i++){
            if(l_v[i] > r_v[i]){
                bbox[0][i] = r_v[i];
                bbox[1][i] = l_v[i];
            }
            else{
                bbox[0][i] = l_v[i];
                bbox[1][i] = r_v[i];
            }
        }
    }


    if(l_child_type == LEAF && r_child_type == BRANCH){

        int* l_v = l_child->leaf->b;
        int** r_v = r_child->branch->b;

        for(int i = 0 ; i < ndim ; i++){
            if(l_v[i] < r_v[0][i]){
                bbox[0][i] = l_v[i];
            }
            else{
                bbox[0][i] = r_v[0][i];
            }

            if(l_v[i] > r_v[1][i]){
                bbox[1][i] = l_v[i];
            }
            else{
                bbox[1][i] = r_v[1][i];
            }
        }
    }

    if(l_child_type == BRANCH && r_child_type == LEAF){

        int** l_v = l_child->branch->b;
        int* r_v = r_child->leaf->b;

        for(int i = 0 ; i < ndim ; i++){
            if(l_v[0][i] < r_v[i]){
                bbox[0][i] = l_v[0][i];
            }
            else{
                bbox[0][i] = r_v[i];
            }

            if(l_v[1][i] > r_v[i]){
                bbox[1][i] = l_v[1][i];
            }
            else{
                bbox[1][i] = r_v[i];
            }
        }
    }

    if(l_child_type == BRANCH && r_child_type == BRANCH){

        int** l_v = l_child->branch->b;
        int** r_v = r_child->branch->b;

        for(int i = 0 ; i < ndim ; i++){
            if(l_v[0][i] < r_v[0][i]){
                bbox[0][i] = l_v[0][i];
            }
            else{
                bbox[0][i] = r_v[0][i];
            }

            if(l_v[1][i] > r_v[1][i]){
                bbox[1][i] = l_v[1][i];
            }
            else{
                bbox[1][i] = r_v[1][i];
            }
        }
    }
    return bbox;
}


Node* RCTree::find_duplicate(int* point){

    Node* nearest = _query(point,root); // we know it is LEAF type
    
    Leaf* l = nearest->leaf;

    for(int i = 0; i< ndim; i++){
        if(point[i] != l->x[i]) return NULL;
    }

    return nearest;
}


Node* RCTree::_query(int* point,  Node* node){

    if(node->type == LEAF)
        return node;
    else{
        Branch* b = node->branch;
        if(point[b->q] <= b->p ) //b->q = the dimension of cut. b->p = value of the cut.
            return _query(point,node->l);
        else
            return _query(point,node->r);
    }
}

void RCTree::_update_leaf_count_upwards(Node* node,int inc){
    while( node ){
        if(node->type == LEAF){
            node->n += inc;
            node = node->u;
        }
        else{
            node->n += inc;
            node = node->u;
        }
    } 
}

// called when deleting the point.
// resize all the bbox above the deleted point 
// -> check if point is defined as boundary.
void RCTree::_relax_bbox_upwards(Node* node, int* point){
    while(node){
        int flag = 0;
        int** node_bbox = node->branch->b;
        int** bbox = _lr_branch_bbox(node); // this node is always branch
        for(int i = 0; i < ndim ; i++){
            if( node_bbox[0][i] == point[i] || node_bbox[1][i] == point[i]){ //if there is same element
                flag = 1;
                break;
            } 
        }
        if(flag == 0){ 
            free(bbox[0]) ;
            free(bbox[1]) ;
            free(bbox);
            break;  // there are no element in point that is boundary. 
        }
        else{
            for(int i = 0; i < ndim ; i++){
                node_bbox[0][i] = bbox[0][i];
                node_bbox[1][i] = bbox[1][i];
            }
            free(bbox[0]) ;
            free(bbox[1]) ;
            free(bbox);
            node = node->u;
        }
    }

}

int RCTree::_get_maxdepth(){
    map<int,Node*>::iterator it;
    int maxdepth = 0;
    for(it = leaves.begin(); it != leaves.end(); it++){
        int curr_depth = it->second->leaf->d;
        if(maxdepth < curr_depth) maxdepth = curr_depth;
    }
    return maxdepth;
}

float RandomFloat(float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b;
    float r = random * diff;
    return r;

}

tuple <int,float> RCTree::_insert_point_cut(int* point, int* bbox){
    int cut_dimension;
    float cut;
    int* bbox_hat[2];
    bbox_hat[0] = (int*) malloc(sizeof(int)*ndim);
    bbox_hat[1] = (int*) malloc(sizeof(int)*ndim);      // bbox_hat[1] is upper bound
    int* b_span = (int*) malloc(sizeof(int)*ndim); 

    for( int i = 0; i < ndim ; i++ ){
        int p = point[i];
        int b = bbox[i];
        if( p > b){
            bbox_hat[0][i] = b;
            bbox_hat[1][i] = p;
            b_span[i] = p-b;
        }
        else{
            bbox_hat[1][i] = b;
            bbox_hat[0][i] = p;
            b_span[i] = b-p;
        }
    }

    int b_range = 0;
    int* span_sum;
    span_sum = (int*) malloc(sizeof(int)*ndim); 

    for( int i = 0 ; i < ndim ; i++){
        b_range += b_span[i];
        span_sum[i] = b_range;
    }

    // // THIS MAY BE VERY SLOW - check later
    float r = RandomFloat(b_range);

    for(int i = 0 ; i < ndim; i++){
        if(span_sum[i] >= r){
            cut_dimension = i;
            break;
        }
    }
    cut = bbox_hat[0][cut_dimension] + span_sum[cut_dimension] - r;

    free(bbox_hat[0]);
    free(bbox_hat[1]);
    free(b_span);
    free(span_sum);

    return make_tuple(cut_dimension,cut);

}

tuple <int,float> RCTree::_insert_point_cut(int* point, int** bbox){
    int cut_dimension;
    float cut;

    int* bbox_hat[2];
    bbox_hat[0] = (int*) malloc(sizeof(int)*ndim);  // bbox_hat[0] is lower bound
    bbox_hat[1] = (int*) malloc(sizeof(int)*ndim);  // bbox_hat[1] is upper bound
    
    int* b_span = (int*) malloc(sizeof(int)*ndim);

    for( int i = 0; i < ndim ; i++ ){
        int p = point[i];
        int b_l = bbox[0][i];
        int b_h = bbox[1][i];
        if( p > b_h){
            bbox_hat[1][i] = p;
        }
        else{
            bbox_hat[1][i] = b_h;
        }

        if( p < b_l){
            bbox_hat[0][i] = p;
        }
        else{
            bbox_hat[0][i] = b_l;
        }
        b_span[i] = bbox_hat[1][i] - bbox_hat[0][i] ;

    }

    int b_range = 0;
    int* span_sum = (int*) malloc(sizeof(int)*ndim);
     
    for( int i = 0 ; i < ndim ; i++){
        b_range += b_span[i];
        span_sum[i] = b_range;
    }
    float r = RandomFloat(b_range);


    for(int i = 0 ; i < ndim; i++){
        if(span_sum[i] >= r){
            cut_dimension = i;
            break;
        }
    }

    cut = bbox_hat[0][cut_dimension] + span_sum[cut_dimension] - r;

    free(bbox_hat[0]);
    free( bbox_hat[1]);
    free( b_span);
    free( span_sum);
    return make_tuple(cut_dimension,cut);


}

void RCTree::_decrement_depth(Node* n, int dec){ //default dec = -1
    n->leaf->d += dec;
}

void RCTree::_increment_depth(Node* n, int inc){ //default inc = 1
    n->leaf->d += inc;
}
void RCTree::print_tree(){
    _print_tree(root);
}

void RCTree::_print_tree(Node* node){
    if(node->type == LEAF){
        cout<<"LEAF - index : " << node->leaf->i << " count : "<< node->n << " depth : " <<node->leaf->d << "\n";
        return;
    } 
    if(node->type == BRANCH){
        cout <<"BRANCH - cut dimension : " <<node->branch->q <<" | cut value " << node->branch->p << " \n";
    }
    _print_tree(node->l);
    _print_tree(node->r);
}

void RCTree::_remove_all(Node* node){


    if(node == NULL) return;
    _remove_all(node->l);
    _remove_all(node->r);
    
    if(node->type == LEAF){
	leaves.erase(node->leaf->i);
	free(node->leaf->x);
	delete node->leaf;
	delete node;
    }
    else{

	free( node->branch->b[0]);
	free( node->branch->b[1]);
	free( node->branch->b);
	delete node->branch;
	delete node;	

    }





}
