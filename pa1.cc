#include <iostream>
#include <utility>
#include <vector>
#include <array>
#include <cstdlib>
#include <ctime>
#include <stdlib.h>
#include <cmath>

using namespace std;

class indexed_priority_queue {
    public:
    float ki;
    float value;
    float values[];
    int pm[];
    float im[];
    int sz;

    bool less(int i, int j) {
        return values[im[i]] < values[im[j]];
    }
    ;
    void swim(int i){
        for (int j = (i - 1)/2; i > 0 && less(i, j); j = (i - 1) / 2) {
            swap(i, j);
            i = j;
        }
    };

    void sink(int i) {
        while (true) {
            int left = 2 * i + 1;
            int right = 2 * i + 2;
            int smallest = left;
            if (right < sz && less(right, left)) {
                smallest = right;
            }
            if (left >= sz || less(i, smallest)) {
                break;
            }
            swap(smallest, i);
            i = smallest;
        }
    };

    void update(int ki, float value) {
        int i = pm[ki];
        values[ki] = value;
        sink(i);
        swim(i);
    };

    void increaseKey(int ki, float value) {
        if (less(values[ki], value)) {
            values[ki] = value;
            sink(pm[ki]);
        }
    };

    void decreaseKey(int ki, float value) {
        if (less(value, values[ki])) {
            values[ki] = value;
            swim(pm[ki]);
        }
    };

    void swap(int i, int j) {
        pm[im[j]] = i;
        pm[im[i]] = j;
        float tmp = im[i];
        im[i] = im[j];
        im[j] = tmp;
    }

    void insert(int ki, float value) {
        values[ki] = value;
        //sz is current size of heap
        pm[ki] = sz;
        im[sz] = ki;
        swim(sz);
        sz = sz + 1;
    }

    void remove(int ki) {
        int i = pm[ki];
        swap(i, sz);
        sink(i);
        swim(i);
        values[ki] = -1;
        pm[ki] = -1;
        im[sz] = -1;
    }
};


// void genEdge(vector <pair<int, int> > adj[], int u, int v, int wt) {
//     adj[u].push_back(make_pair(v, wt));
//     adj[v].push_back(make_pair(u, wt));
// }

struct nodes {
    int from; 
    int to;
    float weight;
} ;

vector<nodes> node;

float rand_sample() {
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

numpoints= 5;
float vertices[numpoints][2];

//represent graph as an array in which the rows = numpoints columns = elements of euclidian dis calculation
float *graph_2d(int dim, int numpoints, float vertices[][2]) {
    for (int i = 0; i < numpoints; i++) {
        for (int j = 1; j < dim; j ++) {
            vertices[i][j] = rand_sample();
        }
    } 
    return reinterpret_cast<float *>(vertices);
}


//represent 2-4 dim graph as an adjacency matrix
void adj_mat(vector <pair<int, float> > adj[], float v1[], float v2[], int dim) {
    for (int i = 0; i < numpoints; i++) {
        for (int j = 0; j < numpoints; j++) {
            if (j != i) {
                float dist = eucl_dist(vertices[i], vertices[j], dim)
                adj[i].push_back(make_pair(j, wt));
                adj[j].push_back(make_pair(i, wt));
            }
        }
    }
}


//modiy one d case to generate edges on the fly 
//represent graph as an array where rows and columns = vertices
//weights = values (adjacency matrix representation)
float *graph_od(int numpoints) {
    /* initialize random seed: */
    srand (time(NULL));
    float vertices[numpoints][numpoints];
    for (int i = 0; i < numpoints; i++) {
        for (int j = 1; j < numpoints; j ++) {
            if (i == j) {
                vertices[i][j] = 0;
            }
            int val = rand_sample();
            vertices[i][j] = val;
            vertices[j][i] = val;
        }
    }
    return reinterpret_cast<float *>(vertices); 
}


int numpoints = 5; //number of nodes in graph
ipq = ;
float G_O = graph_od(numpoints); //graph

    //generate list indicating that all n vertices have not been visited
    bool gen_visit(int numpoints) {
        bool visited[numpoints];
        for (int i = 0; i < numpoints; i++) {
            visited[i] = false;
        }
        return visited;
    }

bool visited = gen_visit(numpoints);

    float Prims() {

        int m = numpoints - 1; //final number of edges in MST
        int edge_num = 0;
        int mst_count = 0;
        bool mstEdges[m];
        float graph[numpoints][numpoints];
        vector<nodes> node;
        int dim = 2;
        int s = 0; //starting vertex
        for (int i = 0; i < m; i++) {
            mstEdges[i] = false;
        }
        
        //in first case, vertex is 0, iterate through 
        for (int i = 0; i < numpoints; i++) {
            if (i != 0) {
                node.push_back(nodes());
                node[i].from = s;
                node[i].to = i;
                node[i].weight = eucl_dist(graph[i], graph[s], dim);
            }
            if (i == numpoints - 1) {
                i = -1;
            }
        }
    }

    pair eager_Prims(s = 0) {
        int m = n - 1;
        int edgeCount, mstCount = 0, 0;
        mstEdges[m];
        for (int i = 0; i < m; i++) {
            mstEdges[i] = false;
        }
        
        relax(s);

        while (!ipq.isEmpty() and edgeCount != m) {
            destNodeIndex, edge = ipq.dequeue();

            mstEdges[edgeCount++] = edge;
            mstCost += edge.cost;

            relax(destNodeIndex);
        }
        if (edgeCount != m) {
            return (false, false);
        }
        return (mstCost, mstEdges);
    }

    void relax(int s) {
        visited[s] = true;
        edges = graph[s];
        for (int i = 0; i < edges; i ++) {
            destNodeIndex = edgeto

            if visited[destNodeIndex]{
                continue;
            }

            if !ipq.contains(destNodeIndex) {
                ipq.insert(destNodeIndex, edge);
            }
            else {
                ipq.decreaseKey(destNodeIndex, edge);
            }
        }
    }


    //so overall, 
    //we use an indexed priority queue implementation of prim's alg
    //this reduces runtime significantly
    //bascially think about MST
    //you'll always end up with exactly 1 incoming edge for each vertex
    //except for the node that you started out with
    //so we iterate through using an indexed priority queue structure
    //for this explanation, vertices refers to points in our original graph while nodes 
    //refers to vertices in our MST
    //so we map a row of our matrix (this is an array of 2-4 elts) (this is generated by our graph functions)
    //to some key aka like the first node (0.1, 0.2) is node #1
    //then we map that key to some value aka the calculated edge weight
    //im thinking rn that we can store this as like a struct (so maybe we convert our 
    //array information into a struct so we can easily refer to individual elements?? idk)
    

    //then we have an im (inverse map) array where, knowing the node in the MST, you can determine what the key is for the original grpah
    //then we gave a pm (position map) where knowing a particular key, we can determine the index of the node in the MST
    //ultimately,we traverse through each time
    //find first vertex, calc all edges, pick shortest edge and add to ipq
    //here, the "name" is the index of the node that the edge is coming into (we can use the array here)
    //this maps to (where its coming from, this current node, cost) 
    //we then go to this node 
    //add all of node 2's edges to the ipq
    //once you add to prioqueue, we keep going through and checking that third elt, the cost, 
    //if an edge going to the same elt (aka having same "name") has lower cost, we change its value 