/*
 *	CSR Representation of graph
 *
 *	Maintains two arrays, edgeArray and vertexArray
 *	edgeArray[|E|] contains the dest vertex for each edge, sorted by source
 *	vertexArray[|V|] contains offsets into edgeArray. Each i in vertexArray 
 *	is the index into edgeArray of the first outgoing edge of vertex i
 */

typedef struct graph {
	int numVertices;					// The number of vertices in the graph
	int numEdges;						// The number of edges in the graph

	double avgOutDegree;
	int maxOutDegree;
	int maxDegreeVertex;
	int *edgeArray;
	int *vertexArray;
	int *edgeLabels;
} Graph;

/**
 *	Initialises a graph of |V| = n and |E| = m, with the given source and dest vertices
 */
Graph * initGraph(int n, int m, int *srcVerts, int *destVerts);

/**
 *	Destroys the given graph
 */
void destroyGraph(Graph *graph);

/**
 *	Returns the out degree of v in G
 */
int outDegree(Graph *graph, int v);

/**
 *	Returns the offset of v in the edgeArray
 */
int edgeOffset(Graph *graph, int v);

/**
 *	Returns edgeArray, pointing to the first adjacent vertex of v
 */
int * adjacentVertices(Graph *graph, int v);

