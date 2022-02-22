#include <iostream>
#include <algorithm>
#include <numeric>

using namespace std;

// int main (int numpoints, int numtrials, int dimension) {

// }

//Prim's pseudocode

/*

    def Prim ()
        // generate some starting vertex s
        
        //the graph (first vertex we need to explore in the graph)
        v_to_explore = [s : 0]

        //empty MST
        MST = {}

        //calculate euclidian distance from s to the following vertices (last 2 cases) 
        //or in first case just randomly generate
        dist_from_s = {s:0,  other vertices:∞}
        //parent list
        prev = {s:null}

        //so while we still have nodes to explore
        while v_to_explore != []
            //delete either the first vertex or the minimum edge weight available
            v = deletemin(v_to_explore)
            
            //any vertex connected to current vertex (not yet in MST)
            for (v,w) in E where w not in MST
                //if current calculated distance to w from s is > 
                if dist_from_s[w] > length(v,w)
                    prev[w] = v
                    dist_from_s[w] = length(v,w)
                    insert (w, dist[w], v_to_explore)

*/



struct nodes {
    int key; //give each vertex a number
    long vertex; //in case where its one dimension
    long* edges; //always n - 1 elements 
    bool visited; //whether or not node has been added to the tree
}
//each vertex of the graph will be a struct (it will contain an array of all edges from that vertex)

nodes* head = nullptr; 
nodes* tail = nullptr; 

long *point = new long(n)
long **vertices = new long*(n)
long A[numpoints][numpoints]

int random (int dimension) {
    for (int i = 0; i < numpoints; i ++) {
        
        for (int j = 0; j < dimension; j ++) {
            point[j] = rand() 1; //rand can't return duplicates
        } 
        vertices[i][j] = point[i]  
    }
}

void one_d () {
    long graph_1 = [numpoints]
    for (int i = 0; i < numpoints; i ++) {
        long v = rand() 1; 
        graph_1[i] = v;
    }
}
//prims understand edges on fly and why can you do that
//prims you might not need a heap, in fact you really need one distance per node without using a heap 
//usually its n^2, takes around 10 hours to run totalish 

//when d = 0 easy, d > 2 different case
//at any point in time don't expect to store every single edge weight (that'll be n^2 so it will take waaay too long) 
//if prim, you can generate edge weights on the fly for the d = 0 case and make sure to justify if you do this 
//maybe compute them on the fly, only compute the edges that you need (generate on the fly all the edges from this vertex)

//if kruskal, you'll have to keep all the edges, consider pruning 
//in an mst you're not gonna use an expensive edge, generate/compute all the edges but only keep 
//use a hard bound for this--intuitively it should depend on n and d for d = 2,3,4
//intuitively in higher d itll be larger bc points are wider, when n is larger theres more density so itll be smaller
//in d = 0, n you draw the graph, look at longest edge on mst and plot that against n and d, how far can we set the cutoff 
//whenever you compute and edge and its too long then you can throw it out
//either you make your upper bound like really leneiet (like keeping twice as many edges)
//orr you set it to be pretty tight, and then 20/100 times it'll break, what is correct way to update from here 

//option #1: 2-D array where we calculate each edge weight and then rep. each vertex (random generation function will take a LONG time)

//option #2: priority queue (implemented using heap)

//option #3 trying to sample from the distribution of edge weights rather than compute them 


//Attempt !! 

int main (int numpoints, int numtrials, int dimension) {
    
    vector<pair<long, long> > adj[numpoints];
    int i, j;
    for (int i = 0; i < numpoints; i ++) {
        long weight = rand() 1;
        adj[i].push_back(make_pair(j, weight));
        adj[j].push_back(make_pair(i, weight));
    }

//every time you only generate the weights for the edges that havent been visited yet 
//basic case!! every time you generate n - i weights and then pick the shortest and add to a list

void one_dim () {
    avg = 0;
    for (k = 0; k < numtrials; k++) {
    //initialize array for final edges of MST (N - 1 edges)
    long MST_0 = [numpoints - 1];
    //iterate through # of vertices
    for (i = 0; i < numpoints; i++) {
        //initalize array for consideration edges from current vertex (aka N - i)
        long considerations = [numpoints - i];
        for (j = 0; j < numpoints - i; j++) {
            //add elements to array for consideration
            considerations[j] = rand() 1;
        }
        //add minimum element from the considerations to the final MST edges
        MST_0[i] = min_element(begin(considerations), end(considerations));
    }
    long sum = 0
    avg += accumulate(MST_0, MST_0+(numpoints-1), sum);
    }
    printf("Average MST_0 weights: %lu\n", avg / numtrials);
}

void two_dim () {
    avg = 0;
    //conduct multiple trials
    for (k = 0; k < numtrials; k++) {
    //initialize array for final edges of MST (N - 1 edges)
    long MST_2 = [numpoints - 1];
    long vertices = [numpoints];
    //create list of vertices (sampling each one from unit square)
    for (i = 0; i < numpoints; i++) {
        vertices[i] = (rand() 1, rand() 1);
    }
    int index = rand() numpoints;
    //iterate through # of vertices
    for (i = 0; i < numpoints; i++) {
        if (i != index) {
        //initalize array for consideration edges from current vertex (aka N - i)
        long considerations = [numpoints - i];
        for (j = 0; j < numpoints - i; j++) {
            //add elements to array for consideration
            considerations[j] = (rand() 1, rand() 1, rand() 1);
        }
        //add minimum element from the considerations to the final MST edges
        MST_2[i] = min_element(begin(considerations), end(considerations));
            
        }
    }
    long sum = 0
    avg += accumulate(MST_2, MST_2+(numpoints-1), sum);
    }
    printf("Average MST_2 weights: %lu", avg / numtrials);
}

//for two dim, think every time we're deletemin the first elemnt of the array 
//for the dimensions we have an array of tuples (each tuple representing a vertex)
//iterate through array, compute weights from 1st elt to all other elts
//then once you get min, add that to a different array and delete 1st elt from the array (or just like iterate to the next elt of the array) 
//keep going until you have nothing left / entire array traversed
//another nested for-loop situation !!! 


//within one run, you compute min edge weight, and 

/*
    for i, i < n, i++
        
        
        
        generate { vertex(range from n-i) - edge (from 0-1) } pairs
        

*/
}

