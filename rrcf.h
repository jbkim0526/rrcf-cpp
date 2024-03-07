#ifndef RRCF_H
#define RRCF_H

#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <numeric>
#include <random>



using namespace std;
#define LEAF 1024
#define BRANCH 2024

class Node;

class Branch{
    /*
    q: Dimension of cut
    p: Value of cut
    l: Pointer to left child
    r: Pointer to right child
    u: Pointer to parent
    n: Number of leaves under branch
    b: Bounding box of points under branch (2 x d)
    */
    public:
        int q;
        float p;
        int** b;
    private:
    public:
        Branch(int q, float p,int** b = NULL);
        

};

class Leaf{
    /*
        i: Index of leaf (user-specified)
        d: Depth of leaf
        u: Pointer to parent
        x: Original point (1 x d) 
        n: Number of points in leaf (1 if no duplicates)
        b: Bounding box of point (1 x d)
    */
    public:
        int i ;
        int s ;
        int d ;
        int* x; // in frequency model the value of array is integer.
        int* b; 

    private:
    public:
        Leaf(int* point, int size, int index, int depth);
};

class Node{
    public:
        int type;
        int n;
        Node* u = NULL;
        Node* l = NULL;
        Node* r = NULL;
        Leaf* leaf = NULL;
        Branch* branch = NULL;

        Node(int type);
};


class RCTree{
    public:
        Node* root;
        int ndim;
        map< int, Node* > leaves;

    private:

    public:
        RCTree(int ndim);
	~RCTree();
        int size();
        Leaf* insert_point(int* point , int index ); // Inserts a point into the tree, creating a Leaf. returns inserted Leaf.
        Node* forget_point(int index); // Delete leaf from tree. returns deleted Leaf.
        float codisp(int index);
        Node* find_duplicate(int* point);
        void map_leaves_increment_depth(Node* node);
        void map_leaves_decrement_depth(Node* node);
        void print_tree();
        //disp;
    private:
        Node* _query(int* point, Node* node);
        void _update_leaf_count_upwards(Node* node,int inc = 1);
        int _get_maxdepth();
        tuple<int,float> _insert_point_cut(int* point, int* bbox);
        tuple<int,float> _insert_point_cut(int* point, int** bbox);
        void _increment_depth(Node* node, int inc = 1);
        void _decrement_depth(Node* node, int dec = -1);
        void _relax_bbox_upwards(Node* node, int* point);
        void _tighten_bbox_upwards(Node* node);
        int** _lr_branch_bbox(Node* node);
        void _print_tree(Node* node);
	void _remove_all(Node* node);
        

};

class RCF{

    public:
        int num_trees;
        int tree_size;
        int vect_size;
        RCTree** rct;
    private:
    public:
        RCF(int num_trees,int tree_size,int vect_size);

};



#endif
