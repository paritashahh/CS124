#include <iostream>
#include <utility>
#include <vector>
#include <array>
#include <cstdlib>
#include <ctime>
#include <stdlib.h>
#include <cmath>
#include <tuple>
#include <stdexcept>
#include <algorithm>
#include <cstdio>
#include <tuple>

using namespace std;

class IndexPrioQueue {

    //pm maps the vertex to the node
    //im maps the node to the vertex
    //value maps the vertex to the weight 
    int SZ;
    vector<int> im;
    vector<int> pm;
    vector<tuple<int, int, float> > values;

    void swap(int i, int j) {
        int temp = im[i];
        im[i] = im[j]; 
        im[j] = temp;
        pm[im[i]] = i; 
        pm[im[j]] = j;
    }

    void shift_up(int i) {
        while(i > 1 && get<2>(values[im[i/2]]) > get<2>(values[im[i]])) {
            swap(i, i / 2);
            i = i / 2;
        }
    }

    void shift_down(int i) {
        int j;
        int child = 2 * i;
        while (child <= SZ) {
            j = child;
            if(j < SZ && get<2>(values[im[j]]) > get<2>(values[im[j+1]])) {
                j++;
            }
            if(values[im[i]] <= values[im[j]]) {
                break;
            }
            swap(i, j);
            i = j;
        }
    }

public:
    //constructor
    IndexPrioQueue() {
        SZ = 0;
        fill(pm.begin(), pm.end(), -1);
    }

    // check whether or not heap is empty
    bool is_empty() {
        return SZ == 0;
    }

    //check whether or not there is a node associated with that key
    bool contains(int ki) {
        return pm[ki] != -1;
    }

    // // size of heap
    // int size() {
    //     return SZ;
    // }

    //insert new pair, making sure to satify heap invariant
    void insert(int i, int ki, float value) {
        SZ++;
        pm[ki] = SZ;
        im[SZ] = ki;
        values[ki] = make_tuple(i, ki, value);
        shift_up(SZ);
    }

    //remove top element from priority queue
    tuple<int, int, float> dequeue() {
        int min = im[1];
        tuple<int, int, float> min_weight = values[min];
        swap(1, SZ--);
        shift_down(1);
        pm[min] = -1;
        im[SZ+1] = -1;
        return min_weight;
    }

    //decrease key value (ensuring heap invariant is satisfied)
    void decreaseKey(int s, int ki, float value) {
        //can do euclidian distnace on the fly here
        if (value < get<2>(values[ki])) {
            values[ki] = make_tuple(s, ki, value);
            shift_up(pm[ki]);
        }
    }
};

//redeed the random sample generator
//generate random sample weight 
float rand_sample() {
	return ((float) rand() / (float) RAND_MAX);
}

//this will calculate euclidian distances
//v1 and v2 will be rows of the graph matrix
float eucl_dist(vector<float> v1, vector<float> v2, int dim) {
    float sum = 0;
    for (int i = 0; i < dim; i++) {
        float a = v2[i] - v1[i];
        sum += a * a;
    }
    return sqrt(sum);
}

//represent graph as an array in which the rows = numpoints 
// columns = elements of euclidian dis calculation
vector<vector<float> > graph_md(int dim, int numpoints) {
    vector<vector<float> > graph(numpoints, vector<float>(dim));
    for (int i = 0; i < numpoints; i++) {
        for (int j = 1; j < dim; j ++) {
            graph[i][j] = rand_sample();
        }
    }
    return graph; 
}

//weights = values (adjacency matrix representation)
vector<vector<float> > graph_od(int numpoints) {
    vector<vector<float> > graph(numpoints, vector<float>(numpoints));
    for (int i = 0; i < numpoints; i++) {
        for (int j = 0; j < numpoints; j++) {
            if (i == j) {
                graph[i][j] = -1;
            }
            int val = rand_sample();
            graph[i][j] = val;
            graph[j][i] = val;
        }
    }
    return graph;
}

//create adjacency matrix for the multidimensional graph
vector<vector<float> > graph_mod(vector<vector<float> > graph, int numpoints, int dim) {
    vector<vector<float> > adj(numpoints, vector<float>(numpoints));
    for (int i = 0; i < numpoints; i++) {
        for (int j = 0; j < numpoints; j ++) {
            //if both vertices same (index i = index j) then no edge between them
            if (i == j) {
                adj[i][j] = -1;
            }
            else {
                //compute euclidian distance between points at row i and j
                int val = eucl_dist(graph[i], graph[j], dim);
                //could include pruning here
                adj[i][j] = val;
                adj[j][i] = val;
            }
        }
    }
    return adj;
}

//put vertex on priority queue
//if there already exists an incoming edge to that vertex, check if current val is < 
//if current val is <, update, else do nothing
//if there doesn't exist an incoming edge to that vertex
//add edge to graph
void relax(int s, int numpoints, bool visited[], IndexPrioQueue prioqueue, vector<vector<float> > graph) {
        visited[s] = true;
        vector<float> outgoing_edges = graph[s];
        for (int i = 0; i < (numpoints - 1); i ++) {
            int destNode = i;
            if (visited[destNode]) {
            continue;
            }
            if (!prioqueue.contains(destNode)) {
                prioqueue.insert(s, destNode, outgoing_edges[i]);
            }
            else {
                prioqueue.decreaseKey(s, destNode, outgoing_edges[i]);
            }
        }
}

//iterate through graph
//add first vertex to priority queue (add all of it's OUTgoing edges too)
//now, pick the smallest edge (aka top of prioqeue) and set that as destination node
//add edge to mst_cost
//now repeat process, making new "s" the destination node
vector<tuple<int, int, float> > eager_Prims(int numpoints, vector<vector<float> > graph) {
    IndexPrioQueue prioqueue;
    int m = numpoints - 1;
    int edge_count = 0;
    //int mst_cost = 0;
    //create array of length numpoints
    //indicate whether or not vertex has been visited
    bool visited[numpoints];
    for (int i = 0; i < numpoints; i++) {
        visited[i] = false;
    }
    //for final mst, vector of tuples
    //node from, node to, edge weight
    vector<tuple<int, int, float> > mstEdges(m);
    int s = 0;
        relax(s, numpoints, visited, prioqueue, graph);

        while (!prioqueue.is_empty() and edge_count != m) {
            tuple<int, int, float> deq = prioqueue.dequeue();
            int srcNode = get<0>(deq);
            int destNode = get<1>(deq);
            float edge = get<2>(deq);

            mstEdges.push_back(deq); //change to push_back
            //mst_cost += edge;

            relax(destNode, numpoints, visited, prioqueue, graph);
        }
        return mstEdges;
}

//function to sum up costs of MST
float sum_cost (vector<tuple<int, int, float> > edges) {
    float sum = 0;
    for (int i = 0; i != edges.size(); i++) {
        tuple<int, int, float> edg = edges[i];
        sum+= get<2>(edg);
    }
    return sum;
}

//function to get average from all trials
float trial_avg (float vals[], int numtrials) {
    float sum = 0;
    for (int i = 0; i < numtrials; i++) {
        sum+= vals[i];
    }
    return sum / numtrials;
}

//function to print adjacency matrix
void print_graph (vector<vector<float> > graph) {
    for (int i = 0; i < graph.size(); i++) {
        for (int j = 0; j < graph[i].size(); j++) {
            cout << graph[i][j];
        }
    }
}

float run_trials (int numpoints, int numtrials, int dimension) {
    //keep track of total cost of mst during each trial 
    float cost[numtrials];
    for (int i = 0; i < numtrials; i++) {
        vector<vector<float> > adj_mat;
        if (dimension == 0) {
            adj_mat = graph_od(numpoints);
        }
        else {
            //generate vertices of graph
            vector<vector<float> > vertices = graph_md(dimension, numpoints);
            //create adjacency matrix and calculate weights
            adj_mat = graph_mod(vertices, numpoints, dimension);
        }
        //use prims alg to find mst
        vector<tuple<int, int, float> > mst = eager_Prims(numpoints, adj_mat);
        //calculate cost of mst
        float sum = sum_cost(mst);
        //depending on dimemsion, add to total cost of mst
    }
    return trial_avg(cost, numtrials);
}

int main() {

    int numpoints;
    int numtrials;
    int dimensions; 
    run_trials(numpoints, numtrials, dimensions);
    return 0;

    //cost should be 0.3 with edges going from 1-2 and 2-3
    vector<vector<float> > test_1 
        {{-1, 0.1, 0.3},
        {0.1, -1, 0.2},
        {0.3, 0.2, -1}};

    //cost should be 0.9 with edges going from 1-2, 2-3, and 3-4
    vector<vector<float> > test_2 
        {{-1, 0.2, 0.5, 0.1},
        {0.2, -1, 0.3, 0.6},
        {0.5, 0.3, -1, 0.4},
        {0.1, 0.6, 0.4, -1}};

    //cost should be 1 with edges going from 1-2, 2-3, 3-4, 4-5
    vector<vector<float> > test_3
        {{-1, 0.1, 0.6, 0.7, 0.5},
        {0.1, -1, 0.2, 0.85, 0.8},
        {0.6, 0.2, -1, 0.3, 0.9},
        {0.7, 0.85, 0.3, -1, 0.4},
        {0.5, 0.8, 0.9, 0.4, -1}};
}
