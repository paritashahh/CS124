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
#include <fstream>
#include <sys/time.h>
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::system_clock;

using namespace std;

long getTime(){
    return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

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
    IndexPrioQueue(int numpoints) {
        SZ = 0;
        im.reserve(numpoints);
        pm.reserve(numpoints);
        values.reserve(numpoints);
        for (int i = 0; i < numpoints; i++) {
            im.push_back(-1);
            pm.push_back(-1);
            values.push_back(make_tuple(-1, -1, -1));
        }
    }

    // check whether or not heap is empty
    bool is_empty() {
        return SZ == 0;
    }

    //check whether or not there is a node associated with that key
    bool contains(int ki) {
        return pm[ki] != -1;
    }

    // size of heap
    int size() {
        return SZ;
    }

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
        //can do euclidian distance on the fly here
        if (value < get<2>(values[ki]) && value != -1) {
            values[ki] = make_tuple(s, ki, value);
            shift_up(pm[ki]);
        }
    }
};

//reseed the random sample generator
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

bool threshold_check (int numpoints, int dim, float val) {
    float threshold = 0;
    if (numpoints <= 1024 && dim == 0) {
        threshold = 0.1;
    }  
    if (numpoints > 1024 && numpoints < 8192 && dim == 0) {
        threshold = 0.02;
    }
     if (numpoints >= 8192 && numpoints < 32768 && dim == 0) {
        threshold = 0.0015;
    }
    if (numpoints >= 32768 && numpoints < 131072 && dim == 0) {
        threshold = 0.0005;
    }
    if (numpoints >= 131072 && numpoints < 262144 && dim == 0) {
        threshold = 0.00011;
    }
    if (numpoints >= 262144 && dim == 0) {
        threshold = 0.00008;
    }
    if (numpoints == 128 && dim == 2) {
        threshold = 0.3;
    }   
    if (numpoints == 128 && dim == 3) {
        threshold = 0.4;
    }
    if (numpoints == 128 && dim == 4) {
        threshold = 0.5;
    }
    if (numpoints == 256 && dim == 2) {
        threshold = 0.2;
    }
    if (numpoints == 256 && dim == 3) {
        threshold = 0.3;
    }
    if (numpoints == 256 && dim == 4) {
        threshold = 0.4;
    }
    if (numpoints == 512 && dim == 2) {
        threshold = 0.1;
    }
    if (numpoints == 512 && dim == 3) {
        threshold = 0.25;
    }
    if (numpoints == 512 && dim == 4) {
        threshold = 0.3;
    }
    if (numpoints == 1024 && dim == 2) {
        threshold = 0.07;
    }
    if (numpoints == 1024 && dim == 3) {
        threshold = 0.2;
    }
    if (numpoints == 1024 && dim == 4) {
        threshold = 0.3;
    }
    if (numpoints == 2048 && dim == 2) {
        threshold = 0.05;
    }
    if (numpoints >= 2048 && dim == 3) {
        threshold = 0.2;
    }
    if (numpoints >= 2048 && dim == 4) {
        threshold = 0.25;
    }
    if (numpoints >= 4096 && dim == 2) {
        threshold = 0.05;
    }
    return val > threshold;
}


//represent graph as an array in which the rows = numpoints 
// columns = elements of euclidian dis calculation
vector<vector<float> > graph_md(int dim, int numpoints) {
    vector<vector<float> > graph(numpoints, vector<float>(dim));
    for (int i = 0; i < numpoints; i++) {
        for (int j = 0; j < dim; j ++) {
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
            float val = rand_sample();
            if (threshold_check(numpoints, 0, val)) {
                graph[i][j] = -1;
            }
            else {
                graph[i][j] = val;
                graph[j][i] = val;
            }
        }
    }
    return graph;
}

void max_weight(vector<vector<float> > graph) {
    ofstream myfile;
    myfile.open("maxweight.txt");
    float max_weight = 0;
    for (int i = 0; i < graph.size(); i++) {
        for (int j = 0; j < graph.size(); j++) {
            if (graph[i][j] > 0) {
                max_weight = graph[i][j];
            }
        }
    }
    myfile << max_weight << endl;
    myfile.close();
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
                float val = eucl_dist(graph[i], graph[j], dim);
                //could include pruning here
                adj[i][j] = val;
                adj[j][i] = val;
            }
        }
    }
    return adj;
}

//function to print adjacency matrix
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

//put vertex on priority queue
//if there already exists an incoming edge to that vertex, check if current val is < 
//if current val is <, update, else do nothing
//if there doesn't exist an incoming edge to that vertex
//add edge to graph
IndexPrioQueue relax(int s, int numpoints, bool visited[], IndexPrioQueue prioqueue, vector<vector<float> > graph, int dim) {
        visited[s] = true;
        vector<float> outgoing_edges = graph[s];
        float val = 0;
        for (int i = 0; i < numpoints; i++) {
            int destNode = i;
            //if we've already visited the destination node, skip this iteration of loop
            float val = 0;
            float threshold = 0;
            if (dim != 0) {
                val = eucl_dist(graph[s], graph[i], dim);
            }
            else {
                val = outgoing_edges[i];
            }
            if (visited[destNode] || s == i || val == -1) {
            continue;
            }
            if (!prioqueue.contains(destNode)) {
                prioqueue.insert(s, destNode, val);
            }
            else {
                 prioqueue.decreaseKey(s, destNode, val);
            }
        }
        return prioqueue;
}

//iterate through graph
//add first vertex to priority queue (add all of it's OUTgoing edges too)
//now, pick the smallest edge (aka top of prioqeue) and set that as destination node
//add edge to mst_cost
//now repeat process, making new "s" the destination node
float eager_Prims(int numpoints, vector<vector<float> > graph, int dim) {
    IndexPrioQueue prioqueue(numpoints);
    int m = numpoints - 1;
    int edge_count = 0;
    float mst_cost = 0;
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
        prioqueue = relax(s, numpoints, visited, prioqueue, graph, dim);

        while (!prioqueue.is_empty() and edge_count != m) {
            tuple<int, int, float> deq = prioqueue.dequeue();
            int srcNode = get<0>(deq);
            int destNode = get<1>(deq);
            float edge = get<2>(deq);
            //printf("edge from %d to %d with weight %f\n", srcNode, destNode, edge);

            mstEdges.push_back(deq);
            mst_cost += edge;

        prioqueue = relax(destNode, numpoints, visited, prioqueue, graph, dim);
        }
        //return mstEdges;
        return mst_cost;
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

float run_trials (int numpoints, int numtrials, int dimension) {
    //keep track of total cost of mst during each trial 
    float cost[numtrials];
    for (int i = 0; i < numtrials; i++) {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        srand((unsigned) tv.tv_usec);
        vector<vector<float> > adj_mat;
        if (dimension == 0) {
            adj_mat = graph_od(numpoints);
            max_weight(adj_mat);
        }
        else {
            //generate vertices of graph
            adj_mat = graph_md(dimension, numpoints);
            //create adjacency matrix and calculate weights
            // adj_mat = graph_mod(vertices, numpoints, dimension);
        }
        //use prims alg to find mst
        ofstream myfile;
        myfile.open("final times.txt");
        long t1 = getTime();
        float mst = eager_Prims(numpoints, adj_mat, dimension);
        long t2 = getTime();
        long prims_time = t2 - t1;
        myfile << "dimension " << dimension << "prim's time ";
        myfile << prims_time << endl;
        cost[i] = mst;
    }
    return trial_avg(cost, numtrials);
}

int main() {

    ofstream myfile;
    myfile.open("final data.txt");
    vector<vector<float> > testing
    {{0.1, 0.2, 0.3},
    {0.25, 0.35, 0.36}, 
    {0.44, 0.15, 0.22}, 
    {0.14, 0.31, 0.21}
    };

 float cool = eager_Prims(4, testing, 3);
    printf("mst cost %f", cool);
    // myfile << "numpoints = 128 dimension = 0 averge" << run_trials(128, 5, 0) << endl;
    // myfile << "numpoints = 128 dimension = 2 averge" << run_trials(128, 5, 2) << endl;
    // myfile << "numpoints = 128 dimension = 3 averge" << run_trials(128, 5, 3) << endl;
    // myfile << "numpoints = 128 dimension = 4 averge" << run_trials(128, 5, 4) << endl;
    //myfile << "numpoints = 256 dimension = 0 averge" << run_trials(256, 1, 0) << endl;
    // myfile << "numpoints = 256 dimension = 2 averge" << run_trials(256, 5, 2) << endl;
    // myfile << "numpoints = 256 dimension = 3 averge" << run_trials(256, 5, 3) << endl;
    // myfile << "numpoints = 256 dimension = 4 averge" << run_trials(256, 5, 4) << endl;
    // myfile << "numpoints = 512 dimension = 0 averge" << run_trials(512, 5, 0) << endl;
    // myfile << "numpoints = 512 dimension = 2 averge" << run_trials(512, 5, 2) << endl;
    // myfile << "numpoints = 512 dimension = 3 averge" << run_trials(512, 5, 3) << endl;
    // myfile << "numpoints = 512 dimension = 4 averge" << run_trials(512, 5, 4) << endl;
    // myfile << "numpoints = 1024 dimension = 0 averge" << run_trials(1024, 5, 0) << endl;
    // myfile << "numpoints = 1024 dimension = 2 averge" << run_trials(1024, 5, 2) << endl;
    // myfile << "numpoints = 1024 dimension = 3 averge" << run_trials(1024, 5, 3) << endl;
    // myfile << "numpoints = 1024 dimension = 4 averge" << run_trials(1024, 5, 4) << endl;
    // myfile << "numpoints = 2048 dimension = 0 averge" << run_trials(2048, 5, 0) << endl;
    // myfile << "numpoints = 2048 dimension = 2 averge" << run_trials(2048, 5, 2) << endl;
    // myfile << "numpoints = 2048 dimension = 3 averge" << run_trials(2048, 5, 3) << endl;
    // myfile << "numpoints = 2048 dimension = 4 averge" << run_trials(2048, 5, 4) << endl;
    // myfile << "numpoints = 4096 dimension = 0 averge" << run_trials(4096, 5, 0) << endl;
    // myfile << "numpoints = 4096 dimension = 2 averge" << run_trials(4096, 5, 2) << endl;
    // myfile << "numpoints = 4096 dimension = 3 averge" << run_trials(4096, 5, 3) << endl;
    // myfile << "numpoints = 4096 dimension = 4 averge" << run_trials(4096, 5, 4) << endl;
    // myfile << "numpoints = 8192 dimension = 0 averge" << run_trials(8192, 5, 0) << endl;
    // myfile << "numpoints = 8192 dimension = 2 averge" << run_trials(8192, 5, 2) << endl;
    // myfile << "numpoints = 8192 dimension = 3 averge" << run_trials(8192, 5, 3) << endl;
    // myfile << "numpoints = 8192 dimension = 4 averge" << run_trials(8192, 5, 4) << endl;
    // myfile << "numpoints = 16384 dimension = 0 averge" << run_trials(16384, 5, 0) << endl;
    // myfile << "numpoints = 16384 dimension = 2 averge" << run_trials(16384, 5, 2) << endl;
    // myfile << "numpoints = 16384 dimension = 3 averge" << run_trials(16384, 5, 3) << endl;
    // myfile << "numpoints = 16384 dimension = 4 averge" << run_trials(16384, 5, 4) << endl;
    //myfile << "numpoints = 32768 dimension = 0 averge" << run_trials(32768, 5, 0) << endl;
    // myfile << "numpoints = 32768 dimension = 2 averge" << run_trials(32768, 5, 2) << endl;
    // myfile << "numpoints = 32768 dimension = 3 averge" << run_trials(32768, 5, 3) << endl;
    // myfile << "numpoints = 32768 dimension = 4 averge" << run_trials(32768, 5, 4) << endl;
    // myfile << "numpoints = 65536 dimension = 0 averge" << run_trials(65536, 5, 0) << endl;
    // myfile << "numpoints = 65536 dimension = 2 averge" << run_trials(65536, 5, 2) << endl;
    // myfile << "numpoints = 65536 dimension = 3 averge" << run_trials(65536, 5, 3) << endl;
    // myfile << "numpoints = 65536 dimension = 4 averge" << run_trials(65536, 5, 4) << endl;
    // myfile << "numpoints = 131072 dimension = 0 averge" << run_trials(131072, 5, 0) << endl;
    // myfile << "numpoints = 131072 dimension = 2 averge" << run_trials(131072, 5, 2) << endl;
    // myfile << "numpoints = 131072 dimension = 3 averge" << run_trials(131072, 5, 3) << endl;
    // myfile << "numpoints = 131072 dimension = 4 averge" << run_trials(131072, 5, 4) << endl;
    //myfile << "numpoints = 262144 dimension = 0 averge" << run_trials(262144, 1, 0) << endl;
    // myfile << "numpoints = 262144 dimension = 2 averge" << run_trials(262144, 5, 2) << endl;
    // myfile << "numpoints = 262144 dimension = 3 averge" << run_trials(262144, 5, 3) << endl;
    // myfile << "numpo ints = 262144 dimension = 4 averge" << run_trials(262144, 5, 4) << endl;
    myfile.close();
    return 0;
    // ofstream myfile;
    // myfile.open("data_generation.txt");
    // int numpoints = 100000;
    // int numtrials = 5;
    // int dimensions = 0; 
    // //graph_od(numpoints);
    // vector<vector<float> > g1 = graph_md(2, numpoints);
    // //graph_mod(g1, numpoints, 2);
    // vector<vector<float> > g2 = graph_md(3, numpoints);
    // //graph_mod(g2, numpoints, 3);
    // vector<vector<float> > g3 = graph_md(4, numpoints);
    // //graph_mod(g3, numpoints, 4);
    // myfile << "numpoints = 5, numtrials = 3, dimensions = 0 case\n";
    // myfile << "avg cost was" << endl;
    // float avg = run_trials(numpoints, numtrials, dimensions);
    // myfile << avg << endl;
    // dimensions = 2; 
    // myfile << "numpoints = 5, numtrials = 3, dimensions = 2 case" << endl;
    // avg = run_trials(numpoints, numtrials, dimensions);
    // myfile << "average weight was " << avg << endl;
    // dimensions = 3; 
    // myfile << "numpoints = 5, numtrials = 3, dimensions = 3 case" << endl;
    // avg = run_trials(numpoints, numtrials, dimensions);
    // myfile << "average weight was " << avg << endl;
    // dimensions = 4; 
    // myfile << "numpoints = 5, numtrials = 3, dimensions = 4 case" << endl;
    // avg = run_trials(numpoints, numtrials, dimensions);
    // myfile << "average weight was " << avg << endl;
    // // myfile << "0D avg " << run_trials(5, 5, 0) << endl;
    // // myfile << "2D avg " << run_trials(5, 5, 0) << endl;
    // // myfile << "3D avg " << run_trials(5, 5, 0) << endl;
    // // myfile << "4D avg " << run_trials(5, 5, 0) << endl;
    // myfile.close();
    // return 0;

    // //cost should be 0.3 with edges going from 1-2 and 2-3
    // vector<vector<float> > test_1 
    //     {{-1, 0.1, 0.3},
    //     {0.1, -1, 0.2},
    //     {0.3, 0.2, -1}};

    // //cost should be 0.15 with edges going from 0-1 and 1-2
    // vector<vector<float> > test_4 
    //     {{-1, 0.1, 0.05},
    //     {0.1, -1, 0.2},
    //     {0.05, 0.2, -1}};

    // //cost should be 0.6 with edges going from 1-2, 2-3, and 3-4
    // vector<vector<float> > test_2 
    //     {{-1, 0.2, 0.5, 0.1},
    //     {0.2, -1, 0.3, 0.6},
    //     {0.5, 0.3, -1, 0.4},
    //     {0.1, 0.6, 0.4, -1}};

    // //cost should be 1 with edges going from 1-2, 2-3, 3-4, 4-5
    // vector<vector<float> > test_3
    //     {{-1, 0.1, 0.6, 0.7, 0.5},
    //     {0.1, -1, 0.2, 0.85, 0.8},
    //     {0.6, 0.2, -1, 0.3, 0.9},
    //     {0.7, 0.85, 0.3, -1, 0.4},
    //     {0.5, 0.8, 0.9, 0.4, -1}};

    // vector<vector<float> > rand_1 = graph_od(3);

    // vector<vector<float> > rand_2 = graph_md(2, 3);

    // vector<vector<float> > rand_3 = graph_md(3, 3);

    // vector<vector<float> > rand_4 = graph_md(4, 3);

    // cout << "test_1";
    // print_graph(test_1);
    // cout << endl;
    // cout << "test_2";
    // print_graph(test_2);
    // cout << endl;
    // cout << "test_3";
    // print_graph(test_3);
    // cout << endl;
    // cout << "testing onedimensional";
    // srand ( time(NULL) );
    // cout << endl;
    // print_graph(rand_1);
    // cout << endl;
    // cout << "hey hey lets do more dimensiosn!";
    // cout << endl;
    // print_graph(rand_2);
    // cout << endl;
    // cout << "a third!!";
    // cout << endl;
    // print_graph(rand_3);
    // cout << endl;
    // cout << "a fourth!!";
    // cout << endl;
    // print_graph(rand_4);
    // cout << endl;
    // cout <<"now lets get some adjacency matrices!!";
    // cout <<endl << endl;
    // cout <<"graph 2d";
    // cout << endl;
    // print_graph(graph_mod(rand_2, 3, 2));
    // cout << endl;
    // cout <<"graph 3d";
    // cout << endl;
    // print_graph(graph_mod(rand_3, 3, 3));
    // cout << endl;
    // cout <<"graph 4d";
    // cout << endl;
    // print_graph(graph_mod(rand_4, 3, 4));

    // IndexPrioQueue prioqueue(5);
    // printf("%d", prioqueue.size());
    // prioqueue.insert(0, 1, 0.1);
    // prioqueue.insert(1, 2, 0.2);
    // prioqueue.insert(2, 3, 0.3);
    // prioqueue.insert(3, 4, 0.05);
    // prioqueue.decreaseKey(3, 2, 0.03);
    // printf("new size aye%d", prioqueue.size());
    // printf("now we r dequeueing!!!\n");
    // tuple<int, int, float> deq = prioqueue.dequeue();
    // int srcNode = get<0>(deq);
    // int destNode = get<1>(deq);
    // float edge = get<2>(deq);
    // printf("%d\n", srcNode);
    // printf("%d\n", destNode);
    // printf("%f\n", edge);

//     vector<vector<float> > testing
//     {{0.1, 0.2},
//     {0.3, 0.4}, 
//     {0.5, 0.6}, 
//     {0.7, 0.8}
//     };

//     float cool = eager_Prims(4, testing, 2);
//     printf("mst cost %f", cool);
}
