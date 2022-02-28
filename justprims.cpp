#include <iostream>
#include <utility>
#include <vector>
#include <array>
#include <cstdlib>
#include <ctime>
#include <stdlib.h>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cstdio>
#include <tuple>
#include <fstream>
#include <sys/time.h>
#include <chrono>
#include <ctime>
#include <stdio.h>    
#include <math.h>  
#include <cmath>
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::system_clock;

using namespace std;

long getTime() {
    return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}


class MinHeapPQ {

    vector<int> expl_vertices;

public:
    vector<float> dist;

    MinHeapPQ(int numpoints) {
        expl_vertices.reserve(numpoints);
        dist.reserve(numpoints);
    }

    int parent(int i) {
        return (i - 1) / 2;
    }

    int lc(int i) {
        return (2 * i) + 1;
    }

    int rc(int i) {
        return (2 * i) + 2;
    }

    void Swap(int i, int j) {
        int temp = expl_vertices[i];
        expl_vertices[i] = expl_vertices[j];
        expl_vertices[j] = temp;
    }

    int Size() {
        return expl_vertices.size();
    }

    bool Empty() {
        return Size() == 0;
    }

    void BubbleDown(int i) {
        int size = Size();
        int left = lc(i);
        int right = rc(i);
        int smallest = i;
        if (left < size && dist[i] > dist[left]) {
            smallest = left;
        }
        if (right < size && (dist[smallest] > dist[right])) {
            smallest = right;
        }
        if (smallest != i) {
            Swap(i, smallest);
            BubbleDown(smallest);
        }
    }

    void BubbleUp(int i) {
        if (Size() == 2) {
            if (dist[expl_vertices[0]] > dist[expl_vertices[1]]) {
                Swap(0, 1);
            }
        }
       // printf("current index %d and parent dist %f and parent index %d\n", i, dist[parent(i)], parent(i));
        if (i && dist[parent(i)] > dist[i]) {
            Swap(parent(i), i);
            BubbleUp(parent(i));
        }
    }

    void Insert(int v) {
        int size = Size();
        expl_vertices.push_back(v);
        BubbleUp(size);
    }

    bool Contains(int v) {
        return (find(expl_vertices.begin(), expl_vertices.end(), v) != expl_vertices.end());
    }

    void DeleteMin() {
        int size = Size();
        if (size == 0) {
            return;
        }
        expl_vertices[0] = expl_vertices.back();
        expl_vertices.pop_back();
        BubbleDown(0);
    }

    int GetMin() {
        return expl_vertices[0];
    }

    void MinHeapify() {
        int size = Size();
        for (int i = size - 1; i >= 0; --i) {
            BubbleDown(i);
        }
    }
};

/////////////////////////////////////////////////////////////////////////////////////////

bool threshold_check (int numpoints, int dim, float val) {
    float threshold = pow((1.0/numpoints),(1.0/dim)); 
    float margin = 3.0;
    return val > (margin * threshold);
}
//generate random sample 
float rand_sample() {
	return ((float) rand() / (float) RAND_MAX);
}

//this will calculate euclidian distances given 2 points represented as a list of floats
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
vector<vector<float> > graph(int dim, int numpoints) {
    vector<vector<float> > graph(numpoints, vector<float>(dim));
    for (int i = 0; i < numpoints; i++) {
        for (int j = 0; j < dim; j++) {
            graph[i][j] = rand_sample();
        }
    }
    return graph; 
}

//function to print graph (or any 2-D array)
void print_graph (vector<vector<float> > graph) {
    ofstream myfile;
    myfile.open("graph_generation.txt");
    for (int i = 0; i < graph.size(); i++) {
        for (int j = 0; j < graph[i].size(); j++) {
            myfile << graph[i][j] << " ";
        }
        myfile << endl;
    }
    myfile.close();
}

int vertices;
float max_flt = INFINITY;
int min_vertex;

void trial () {
    MinHeapPQ prioqueue(5);
    prioqueue.Insert(0);
    prioqueue.Insert(1);
    prioqueue.Insert(2);
    prioqueue.Insert(3);
    prioqueue.Insert(4);
    prioqueue.dist[0] = 0;
    prioqueue.dist[1] = 1;
    prioqueue.dist[2] = 2;
    prioqueue.dist[3] = 3;
    prioqueue.dist[4] = 4;
     while (!prioqueue.Empty()) {
        printf("del node val %d\n", prioqueue.GetMin());
        //printf("del node index %d\n", prioqueue.GetMin());
        prioqueue.DeleteMin();
        //printf("del node new val %f\n", prioqueue.dist[0]);
    }
}


float prims (vector<vector<float> > graph, int numpoints, int dimension) {
    //this is my priority queue here
    MinHeapPQ prioqueue(numpoints);
    int prev[numpoints];
    bool visited[numpoints];
    int edgeCount = 0;
    float mstCost = 0;

    for (int i = 0; i < numpoints; i++) {
        prioqueue.dist[i] = max_flt;
        visited[i] = false;
        prev[i] = -1;
    }

    int v = 0;
    visited[0] = true;

    while(edgeCount != numpoints - 1) {
        if (prioqueue.Size() != 0) {
           // printf("size is %d\n", prioqueue.Size());
            v = prioqueue.GetMin();
            // printf("min vertex is %d\n", v);
            prioqueue.DeleteMin();
            // printf("now min vertex is %d\n", prioqueue.GetMin());
            //printf("mst edge from vertex %d to vertex %d\n", v, prev[v]);
            visited[v] = true;
            mstCost += prioqueue.dist[v];
            //printf("added edge to %d with value %f to the MST\n", v, prioqueue.dist[v]);
            edgeCount++;
        }
       for (int i = 1; i < numpoints; i++) {
            float edge_dist = eucl_dist(graph[v], graph[i], dimension);
            ///might need to reconsider WHEN we are pruning here!!!
            if (i == v || visited[i] || threshold_check(numpoints, dimension, edge_dist)) {
               continue;
            }
            //printf("calculated distance to vertex %d: from vertex %d is %f\n", i, v, edge_dist);
            if (!prioqueue.Contains(i))  {
                prioqueue.dist[i] = edge_dist;
                prev[i] = v;
                prioqueue.Insert(i);
               // printf("added vertex %d with edge weight %f\n", i, prioqueue.dist[i]);
               // printf("this is min vertex %d with weight %f in prioqueue rn\n", prioqueue.GetMin(), prioqueue.dist[prioqueue.GetMin()]);
            }
            else if (prioqueue.dist[i] > edge_dist) {
                prioqueue.dist[i] = edge_dist;
                prev[i] = v;
                prioqueue.MinHeapify();
           }
       }
    }
    // for (int i = 0; i < prioqueue.dist.size(); i++) {
    //     printf("value at index %d is %f\n", i, prioqueue.dist[i]);
    // }
    assert (edgeCount == numpoints - 1);
    return mstCost;
}

float run_trials (int numpoints, int dimension, int numtrials) {
    float trial_avg = 0;
    long graph_avg = 0;
    long prims_avg = 0;
    for (int i = 0; i < numtrials; i++) {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        srand((unsigned) tv.tv_usec);

        long gtime1 = getTime();
        vector<vector<float> > gener_graph = graph(dimension, numpoints);
        long gtime2 = getTime();

        long ptime1 = getTime();
        float mstcost = prims(gener_graph, numpoints, dimension);
        long ptime2 = getTime();

        trial_avg += mstcost;
        graph_avg += (gtime2-gtime1);
        prims_avg += (ptime2-ptime1);
    }
    // ofstream myfile;
    // myfile.open("graph_generation.txt");
    // myfile << "graph_avg" << (graph_avg / numtrials) << endl;
    // myfile << "prims_avg" << (graph_avg / numtrials) << endl;
    // myfile.close();
    return (trial_avg / numtrials);
}

int main () {
    printf("2D Trials in order from smallest numpoints to largest");
    printf("%f\n", run_trials(128, 2, 5));
    printf("%f\n", run_trials(256, 2, 5));
    printf("%f\n", run_trials(512, 2, 5));
    printf("%f\n", run_trials(1024, 2, 5));
    printf("%f\n", run_trials(2048, 2, 5));
    printf("%f\n", run_trials(4096, 2, 5));
    printf("%f\n", run_trials(8192, 2, 5));
    printf("%f\n", run_trials(16384, 2, 5));
    printf("%f\n", run_trials(32768, 2, 5));
    printf("%f\n", run_trials(65536, 2, 5));
    printf("%f\n", run_trials(131072, 2, 5));
    printf("%f\n", run_trials(262144, 2, 5));
    // vector<vector<float> > testing{{0.1, 0.2, 0.3},
    // {0.25, 0.35, 0.36}, 
    // {0.44, 0.15, 0.22}, 
    // {0.14, 0.31, 0.21}
    // };

    // vector<vector<float> > testing2{{0.2, 0.3},
    // {0.1, 0.2}, 
    // {0.45, 0.68}};

    // vector<vector<float> > testing3{{0.1, 0.2, 0.3},
    // {0.5, 0.6, 0.7}, 
    // {0.9, 0.8, 0.7}};

    // vector<vector<float> > testing4{{0.2}, {0.1}, {0.8}, {0.9}};

    // float od = prims(graph(1, 300000), 300000, 1);
    // printf("od mst cost %f\n", od);
    // float oad = prims(testing, 4, 3);
    // printf("oad mst cost %f\n", oad);
    
    // float biad = prims(testing2, 3, 2);
    // printf("biad mst cost %f\n", biad);

    // float triad = prims(testing3, 3, 3);
    // printf("triad mst cost %f\n", triad);
    return 0;
}