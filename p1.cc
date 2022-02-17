int main (int numpoints, int numtrials, int dimension) {

}

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

//option #1: 2-D array where we calculate each edge weight and then rep. each vertex (random generation function will take a LONG time)

//option #2: priority queue (implemented using heap)

//option #3 trying to sample from the distribution of edge weights rather than compute them 