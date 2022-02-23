// int main (N, m, a, b, wt) {
//     int N, m
//     vector<pair<int, int> > adj[N]

//     int a, b, wt;
//     for (int i = 0; i < m; i++) {
//         adj[a].push_back(make_pair(b, wt));
//         adj[b].push_back(make_pair(a, wt));
//     }
//     int parent[N]
//     int key[N]
//     bool mstSet[N]

//     for (int i = 0; i < N; i++){
//         key[i] = INT_MAX, mstSet[i] = false, parent[i] = -1
//     }

//     key[0] = 0;
//     parent[0] = -1;

//     for (int count = 0; count < N - 1; count ++) {
//         int mini = INT_MAX, u;
        
//         for (int v = 0; v < N; v++) {
//             if (mstSet[v] == false && key[v] = mini) {
//                 mini = key[v], u = v;
//             }
//         }

//         mstSet[u] = true;

//         for (auti it : adj[u]) {
//             int v = it.first;
//             int weight = it.second;
//             if (mstSet[v] == false && weight < key[v]) {
//                 parent[v] = u, key[v] = weight;
//             }
//         }
//     }
// }

#include <iostream>
#include <vector>
#include <array>
using namespace std;


int main () {
    array<pair<long, long> 2D_Graph[100]
    for (int i = 0); i < 100; i++) {
        2D_Graph[i].push_back(make_pair(rand () 1, rand () 1))
    }
}

void graph_gen () {
    for (k = 0; k < numtrials; k++) {
    //initialize array for vertices
    vector<pair<long, long> > 2D_Graph[numpoints]
    for (int i = 0; i < numpoints; i++) {
        2D_Graph[i].push_back(make_pair(rand () 1, ran () 1));
    }
    //iterate through # of vertices
    for (i = 0; i < numpoints; i++) {
        //initalize array for consideration edges from current vertex (aka N - i)
        long considerations = [numpoints - i];
    }
    }
}


