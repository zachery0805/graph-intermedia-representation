#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graph.h"


/**
 *	Initialises a graph of |V| = n and |E| = m, with the given source and dest vertices
 */
Graph * initGraph(int n, int m, int *srcVerts, int *destVerts) {
	Graph *graph = (Graph *) malloc(sizeof(Graph));

	graph->numVertices = n;
	graph->numEdges = m;
	graph->avgOutDegree = (double) m / (double) n;

	// Initialise arrays for CSR
	graph->avgOutDegree = 0;
	graph->edgeArray = (int *) malloc(sizeof(int) * m);
	graph->edgeLabels = (int *) malloc(sizeof(int) * m);
	graph->vertexArray = (int *) malloc(sizeof(int) * (n + 1));
	graph->maxOutDegree = 0;
	graph->maxDegreeVertex = 0;
	
	// Vertex 0 always starts at 0 in edgeArray
	graph->vertexArray[0] = 0;

	// Work out the outDegrees of each vertex to find offsets for vertexArray
	int *outDegrees = (int *) malloc(sizeof(int) * n);
	for(int i = 0; i < n; i++) {
		outDegrees[i] = 0;
	}
	for(int i = 0; i < m; i++) {
		outDegrees[srcVerts[i]]++;
	}

	// Calculate vertexArray values
	for(int i = 0; i < n; i++) {
		graph->vertexArray[i + 1] = graph->vertexArray[i] + outDegrees[i];

		// Check if this is the largest out degree so far
		if(outDegrees[i] > graph->maxOutDegree) {
			graph->maxOutDegree = outDegrees[i];
			graph->maxDegreeVertex = i;
		}
	}

	// Calculate edgeArray values
	for(int i = 0; i < m; i++) {
		graph->edgeArray[i] = destVerts[i];
	}


	free(outDegrees);
	return graph;
}

/**
 *	Destroys the given graph by freeing edge and vertex arrays
 */
void destroyGraph(Graph *graph) {
	free(graph->edgeArray);
	free(graph->vertexArray);
	free(graph->edgeLabels);
	free(graph);
}

/**
 *	Returns the out degree of v in G
 */
int outDegree(Graph *graph, int v) {
	return graph->vertexArray[v+1] - graph->vertexArray[v];
}

/**
 *	Returns the offset of v in the edgeArray
 */
int edgeOffset(Graph *graph, int v) {
	return graph->vertexArray[v];
}

/**
 *	Returns edgeArray, pointing to the first adjacent vertex of v
 */
int * adjacentVertices(Graph *graph, int v) {
	return &graph->edgeArray[graph->vertexArray[v]];
}