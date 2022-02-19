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
