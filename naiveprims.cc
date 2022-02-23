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

using namespace std;

int main() {
    int numpoints = 128;
	vector<pair<int,float> > adj[numpoints]; 

	int a,b,wt;
	for(int i = 0; i<m ; i++){
		cin >> a >> b >> wt;
		adj[a].push_back(make_pair(b,wt));
		adj[b].push_back(make_pair(a,wt));
	}	
	
	int prev[N]; 
      
    int key[N]; 
      
    bool mstSet[N]; 
  
    for (int i = 0; i < N; i++) 
        key[i] = INT_MAX, mstSet[i] = false; 
  

    key[0] = 0; 
    prev[0] = -1; 
    int ansWeight = 0;
    for (int count = 0; count < N - 1; count++)
    { 
        
        int mini = INT_MAX, u; 
  
        for (int v = 0; v < N; v++) 
            if (mstSet[v] == false && key[v] < mini) 
                mini = key[v], u = v; 
                
        mstSet[u] = true; 
        
        for (auto it : adj[u]) {
            int v = it.first;
            int weight = it.second;
            if (mstSet[v] == false && weight < key[v]) 
                prev[v] = u, key[v] = weight; 
        }
            
    } 
    
    for (int i = 1; i < N; i++) 
        cout << prev[i] << " - " << i <<" \n"; 
	return 0;
}


//generate random sample weight 
float rand_sample() {
	return ((float) rand() / (float) RAND_MAX);
}

//weights = values (adjacency matrix representation)
void graph_od(int numpoints, float** graph) {
    for (int i = 0; i < numpoints; i++) {
        for (int j = 1; j < numpoints; j ++) {
            if (i == j) {
                graph[i][j] = 0;
            }
            int val = rand_sample();
            graph[i][j] = val;
            graph[j][i] = val;
        }
    }
}