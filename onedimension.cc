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

using namespace std;

class IndexPrioQueue {

    //heap is im, index is pm, keys is value
    //pm maps the vertex to the node
    //im maps the node to the verte
    //value maps the vertex to the weight 
    int MAX, sz, *values, *im, *pm;

    //compare values, return true if value at index i < value at index j
    //compare values for different keys
    bool comp(int i, int j) {
        return (values[im[i]] < values[im[j]]);
    }

    //swap the position of two values in the heap
    void swap(int i, int j) {
        pm[im[j]] = i;
        pm[im[i]] = j;
        float tmp = im[i];
        im[i] = im[j];
        im[j] = tmp;
    }

    //move up in heap
    void shift_up(int i) {
        for (int j = parent(i); i > 0 && comp(i, j); j = parent(i)) {
            swap(i, j);
            i = j;
        }
    }

    //move down in heap
    void shift_down(int i) {
        while (true) {
            int smallest = left_kid(i);
            if (right_kid(i) < sz && comp(right_kid(i), left_kid(i))) {
                smallest = right_kid(i);
            }
            if (left_kid(i) >= sz || comp(i, smallest)) {
                break;
            }
            swap(smallest, i);
            i = smallest;
        }
    }

public:
    //index priority queue initialization
    //heap can contain at max n - 1 elements 
    //constructor
    IndexPrioQueue(int MAX) {
        this->MAX = MAX;
        sz = 0;
        values = new int[MAX + 1];
        im = new int[MAX + 1];
        pm = new int[MAX + 1];
        for (int i = 0; i < (MAX + 1); i++) {
            pm[i] = -1;
        }
    }

    //destructor
    ~IndexPrioQueue() {
        delete [] values;
        delete [] im;
        delete [] pm;
    } 

    //check whether or not there is a node associated with that key
    bool contains(int i) {
        if (pm[i] == -1) {
            return false;
        }
        else {
            return true;
        }
    }

    //return size of heap
    // int size() {
    //     sz;
    // }

    //get parent of node
    int parent(int i) {
        return ((i - 1) / 2);
    }

    //get left child of node
    int left_kid(int i) {
        return ((2 * i) + 1);
    }

    //get right child of node
    int right_kid(int i) {
        return ((2 * i) + 2);
    }

    //check whether or not the heap is empty
    bool is_empty() {
        return sz == 0;
    }

    //insert new pair, making sure to satify heap invariant
    void insert(int ki, float value) {
        values[ki] = value;
        pm[ki] = sz;
        im[sz] = ki;
        shift_up(sz);
        sz++;
    }

    //delete pair, making sure to satisfy heap invariant
    void remove(int ki) {
        int i = pm[ki];
        swap(i, sz);
        shift_down(i);
        shift_up(i);
        values[ki] = -1;
        pm[ki] = -1;
        im[sz] = -1;
    }

    //increase key value
    void increaseKey(int ki, float value) {
        values[ki] = value;
        shift_down(pm[ki]);
    }

    //decrease key value
    void decreaseKey(int ki, float value) {
        values[ki] = value;
        shift_up(pm[ki]);
    }

    //update value, then make sure heap invariant is still satisfied
    void update(int ki, float value) {
        int i = pm[ki];
        values[ki] = value;
        shift_down(i);
        shift_up(i);
    }

    pair<int, float> dequeue() {
        int min = im[1];
        int min_weight = values[min];
        swap(1, sz);
        shift_down(1);
        im[min] = -1;
        pm[sz] = -1;
        return make_pair(min, min_weight);
    }
};



float rand_sample() {
    /* initialize random seed: */
    srand (time(NULL));
	return ((float) rand() / (float) RAND_MAX);
}

//store the randomly computed dimensions in an array
//this will calculate euclidian distances
//v1 and v2 will be rows of the graph matrix
float eucl_dist(float v1[], float v2[], int dim) {
    float sum = 0;
    for (int i = 0; i < dim; i++) {
        float a = v2[i] - v1[i];
        sum += a * a;
    }
    return sqrt(sum);
}

int numpoints= 128;
float two_dvertices[128][2];
float three_dvertices[128][3];
float four_dvertices[128][4];
float one_dweights[128][128];

//represent graph as an array in which the rows = numpoints columns = elements of euclidian dis calculation
void graph_2d(int dim, int numpoints) {
    for (int i = 0; i < numpoints; i++) {
        for (int j = 1; j < dim; j ++) {
            two_dvertices[i][j] = rand_sample();
        }
    } 
}

void graph_3d(int dim, int numpoints) {
    for (int i = 0; i < numpoints; i++) {
        for (int j = 1; j < dim; j ++) {
            three_dvertices[i][j] = rand_sample();
        }
    } 
}

void graph_4d(int dim, int numpoints) {
    for (int i = 0; i < numpoints; i++) {
        for (int j = 1; j < dim; j ++) {
            four_dvertices[i][j] = rand_sample();
        }
    } 
}

//represent 2-4 dim graph as an adjacency matrix
// vector <pair<int, float> > adj_mat(vector <pair<int, float> > adj[], float v1[], float v2[], int dim) {
//     for (int i = 0; i < numpoints; i++) {
//         for (int j = 0; j < numpoints; j++) {
//             if (j != i) {
//                 float dist = eucl_dist(vertices[i], vertices[j], dim);
//                 adj[i].push_back(make_pair(j, dist));
//                 adj[j].push_back(make_pair(i, dist));
//             }
//         }
//     }
//     return *adj;
// }

//weights = values (adjacency matrix representation)
void graph_od(int numpoints) {
    float one_dweights[numpoints][numpoints];
    for (int i = 0; i < numpoints; i++) {
        for (int j = 1; j < numpoints; j ++) {
            if (i == j) {
                one_dweights[i][j] = 0;
            }
            int val = rand_sample();
            one_dweights[i][j] = val;
            one_dweights[j][i] = val;
        }
    }
}

bool visited_list(int numpoints) {
    bool visited[numpoints];
     for (int i = 0; i < numpoints; i++) {
    visited[i] = false;
    }
    return visited;
}

bool mst_list(int numpoints) {
    bool mst_Edges[numpoints - 1];
    for (int i = 0; i < m; i++) {
    mst_Edges[i] = false;
    }
    return mst_Edges;
}

IndexPrioQueue prioqueue(numpoints);
int m = numpoints - 1;
int edge_count = 0;
int mst_cost = 0;
bool visited = visited_list(numpoints);
bool mstEdges = mst_list(numpoints);
int s = 0;

void relax(int s) {
        visited[s] = true;
        float* outgoing_edges = one_dweights[s];
        for (int i = 0; i < (numpoints - 1); i ++) {
            int destNode = i;
            if (visited[destNode]) {
            continue;
            }
            if (!prioqueue.contains(destNode)) {
                prioqueue.insert(destNode, outgoing_edges[i]);
            }
            else {
                prioqueue.decreaseKey(destNode, outgoing_edges[i]);
            }
        }
}

float eager_Prims(int numpoints) {
        int s = 0;
        relax(s);

        while (prioqueue.is_empty() and edge_count != m) {
            pair<int, float> deq = prioqueue.dequeue();
            int destNode = deq.first;
            float edge = deq.second;

            mstEdges[edge_count++] = edge;
            mst_cost += edge;

            relax(destNode);
        }
        return mst_cost;
    }
