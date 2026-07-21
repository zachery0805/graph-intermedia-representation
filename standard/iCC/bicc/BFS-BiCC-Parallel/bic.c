#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graph.h"
#include <sys/time.h>
// #include "omp.h"

const int MAX_LEVELS = 1000000;

int numBiconnectedComponents = 0;
int maxBiconnectedComponent = 0;

double timer() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
}


/**
 *  Reads the graph from the given file, and fills the arrays of source and dest vertices
 *  Each edge in the graph is represented by (srcVerts[i], destVerts[i])
 */
Graph * parseGraphFile(char* filename) {
    int n, m;
    int *srcVerts, *destVerts;
    FILE *file = fopen(filename, "r");
    Graph *graph;
    if(file != NULL) {
        char line[1000];
        int lineNum = 0, vertIndex = 0;
        while(fgets(line, sizeof line, file) != NULL) {
            // Strip the trailing new line character
            size_t ln = strlen(line) - 1;
            if(line[ln] == '\n') {
                line[ln] = '\0';
            }

            if(lineNum == 0) {
                // First line, get n and m values
                char *tok = strtok(line, " ");
                n = strtol(tok, NULL, 10);

                tok = strtok(NULL, " ");
                m = strtol(tok, NULL, 10);

                // Initialise srcVerts and destVerts to number of edges
                srcVerts = (int *) malloc(sizeof(int) * m);
                destVerts = (int *) malloc(sizeof(int) * m);

            } else {
                // Edge definition, put into src and dest arrays
                char *tok = strtok(line, " ");
                int src = strtol(tok, NULL, 10);
                srcVerts[vertIndex] = src;

                tok = strtok(NULL, " ");
                int dest = strtol(tok, NULL, 10);
                destVerts[vertIndex] = dest;
                vertIndex++;
            }
            lineNum++;
        }

        // Generate the graph with the edges from the file
        graph = initGraph(n, m, srcVerts, destVerts);
        free(srcVerts);
        free(destVerts);
    }

    fclose(file);
    return graph;
}

/**
 *  Parallel BFS to find P, L, and LQ values of the given graph
 */
void bfs(Graph *graph, int v, int *levels, int *P, int *L, int **LQ, int *LQCounts, int *visited) {
    int numThreads = omp_get_max_threads();

    // List of vertices to visit - appended to during the search
    int *toVisit = malloc(sizeof(int) * graph->numVertices);
    int back = 0;
    int *toVisitThreads = malloc(sizeof(int) * (numThreads * graph->numVertices));
    int *threadOffsets = malloc(sizeof(int) * numThreads);

    P[v] = v;
    L[v] = 0;
    LQ[0] = malloc(sizeof(int));
    LQ[0][0] = v;
    LQCounts[0] = 1;

    toVisit[0] = v;
    back = 1;

    // Direction optimising BFS vars
    double alpha = 15.0;
    double beta = 25.0;
    int useBottomUp = 0;
    int nf = 0;
    int localnf;

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int *toVisitThread = &toVisitThreads[tid * graph->numVertices];
        int threadOffset = 0;
        int level = 1;
        int prevLevel;
        int vert, adjEnd, u;
        int * adjVertices;

        while(back) {
            threadOffset = 0;

            if(!useBottomUp) {
                /* Using top down approach */
                #pragma omp for schedule(static) reduction(+:localnf)
                for(int i = 0; i < back; i++) {
                    vert = toVisit[i];
                    adjVertices = adjacentVertices(graph, vert);
                    adjEnd = outDegree(graph, vert);
                    // Go through each vertex adjacent to v
                    for(int j = 0; j < adjEnd; j++) {
                        u = adjVertices[j];
                        if(L[u] < 0) {
                            L[u] = level;
                            P[u] = vert;
                            // Add the adjacent vertex to the threads buffer
                            toVisitThread[threadOffset++] = u;
                            localnf++;
                        }
                    }
                }
            } else {
                /* Using bottom up approach */
                prevLevel = level - 1;

                #pragma omp for schedule(static) reduction(+:localnf)
                for(int i = 0; i < graph->numVertices; i++) {
                    vert = i;
                    if(L[vert] < 0) {
                        adjVertices = adjacentVertices(graph, vert);
                        adjEnd = outDegree(graph, vert);
                        for(int j = 0; j < adjEnd; j++) {
                            u = adjVertices[j];
                            if(L[u] == prevLevel) {
                                L[vert] = level;
                                P[vert] = u;
                                toVisitThread[threadOffset++] = vert;
                                localnf++;
                                break;
                            }
                        }
                    }
                }
            }

            threadOffsets[tid] = threadOffset;

            #pragma omp barrier

            #pragma omp single
            {
                nf += localnf;

                // Determine whether to switch BFS approaches
                if(useBottomUp) {
                    double mf = (double) localnf * graph->avgOutDegree;
                    double mu = (double) (graph->numVertices - nf) * graph->avgOutDegree;

                    // Switch if heuristic is true, and if we haven't switched yet
                    if(mf > (mu / alpha) && mu > 0) {
                        useBottomUp = 1;
                    }
                } else {
                    if(nf < ((double) graph->numVertices / beta)) {
                        useBottomUp = 0;
                    }
                }

                int oldBack = back;
                back = 0;
                // Copy each threads buffer of vertices into main array
                for(int i = 0; i < numThreads; i++) {
                    if(threadOffsets[i]) {
                        int offset = i * graph->numVertices;
                        for(int j = 0; j < threadOffsets[i]; j++) {
                            toVisit[back++] = toVisitThreads[offset + j];
                        }
                    }
                }

                // Store this level in LQ
                LQ[level] = (int *) malloc(sizeof(int) * back);
                for(int i = 0; i < back; i++) {
                    LQ[level][i] = toVisit[i];
                }
                LQCounts[level] = back;
                (*levels)++;
            }
            level++;

            #pragma omp barrier
        }
    }
    
    free(toVisit);
    free(toVisitThreads);
    free(threadOffsets);
}

void findArticulationPoints(Graph *graph, int *P, int *L, int *Par, int *Low, int **LQ, int *LQCounts, int *Art, int levels) {
    #pragma omp parallel 
    {
        // Each thread maintains it's own visited array
        int *visited = (int *) malloc(sizeof(int) * graph->numVertices);
        for(int i = 0; i < graph->numVertices; i++) {
            visited[i] = 0;
        }

        // Tracks the list of unique vertices encountered
        int *stack = (int *) malloc(sizeof(int) * graph->numVertices);
        int *queue = (int *) malloc(sizeof(int) * graph->numVertices);
        int *nextQueue = (int *) malloc(sizeof(int) * graph->numVertices);
        int maxVert = 0;
        int stackBack;
        int back;
        int nextBack;
        int curLow;
        int isCurrArt;

        for(int l = levels - 1; l > 1; l--) {
            int *currQueue = LQ[l];
            int levelSize = LQCounts[l];

            #pragma omp for schedule(guided)
            for(int v = 0; v < levelSize; v++) {
                int vert = currQueue[v];
                if(Low[vert] > -1) {
                    continue;
                }

                int vertLevel = l;
                int vertParent = P[vert];

                queue[0] = vert;
                back = 1;
                stack[0] = vert;
                stackBack = 1;
                visited[vert] = 1;
                isCurrArt = 0;
                curLow = graph->numVertices;

                while(back) {
                    nextBack = 0;
                    for(int j = 0; j < back; j++) {
                        int newVert = queue[j];
                        int isNewArt = Art[newVert];

                        int *newOuts = adjacentVertices(graph, newVert);
                        int newOutDegree = outDegree(graph, newVert);
                        for(int n = 0; n < newOutDegree; n++) {
                            int newOut = newOuts[n];

                            if(visited[newOut] || newOut == vertParent) {
                                continue;
                            } else if(isNewArt && Low[newOut] > -1) {
                                continue;
                            } else if(L[newOut] < vertLevel) {
                                goto next_out;
                            } else {
                                visited[newOut] = 1;
                                nextQueue[nextBack++] = newOut;
                                stack[stackBack++] = newOut;
                                if(newOut < curLow) 
                                    curLow = newOut;
                            }
                        }
                    }

                    int *tmp = nextQueue;
                    nextQueue = queue;
                    queue = tmp;
                    back = nextBack;
                }

                Art[vertParent] = 1;
                isCurrArt = 1;

                next_out:
                if(isCurrArt) {
                    for(int j = 0; j < stackBack; j++) {
                        int sVert = stack[j];
                        visited[sVert] = 0;
                        Par[sVert] = vert;
                        Low[sVert] = curLow;  
                    }
                } else {
                    for(int j = 0; j < stackBack; j++) {
                        int sVert = stack[j];
                        visited[sVert] = 0;        
                    }
                }
            }
        }

        free(visited);
        free(stack);
        free(queue);
        free(nextQueue);

    }
}

/**
 *  Finds the articulation points of the given graph by inspecting the queues in LQ in reverse order, and 
 *  determining Par and Low values
 */
void findArticulationPoints2(Graph *graph, int *P, int *L, int *Par, int *Low, int **LQ, int *LQCounts, int *Art, int levels) {
    #pragma omp parallel 
    {
        // Each thread maintains it's own visited array
        int *visited = (int *) malloc(sizeof(int) * graph->numVertices);
        for(int i = 0; i < graph->numVertices; i++) {
            visited[i] = 0;
        }

        // Tracks the list of unique vertices encountered
        int *Vu = (int *) malloc(sizeof(int) * graph->numVertices);
        int *queue = (int *) malloc(sizeof(int) * graph->numVertices);
        int *nextQueue = (int *) malloc(sizeof(int) * graph->numVertices);
        int *queueWhole = (int *) malloc(sizeof(int) * graph->numVertices);

        int maxVert = 0;
        int stackEnd;
        int queueEnd, queueEndWhole;
        int nextQueueEnd;
        int isCurrArt;

        // Lowest vertex identifier encountered
        int vidLow;

        // Inspect LQ in reverse order (LQm..1)
        for(int level = levels - 1; level >= 1; level--) {
            int levelSize = LQCounts[level];

            #pragma omp for schedule(guided)
            for(int uIndex = 0; uIndex < levelSize; uIndex++) {
                int u = LQ[level][uIndex];
                
                if(Low[u] == -1) {
                    // Get the parent of u
                    int v = P[u];

                    queue[0] = u;
                    queueEnd = 1;
                    Vu[0] = u;
                    stackEnd = 1;
                    visited[u] = 1;
                    isCurrArt = 0;
                    vidLow = u;

                    int queueIndex = 0;
                    while(queueIndex < queueEnd) {
                        int x = queue[queueIndex++];

                        int isXArt = Art[x];

                        int *xAdj = adjacentVertices(graph, x);
                        int xAdjEnd = outDegree(graph, x);
                        for(int wIndex = 0; wIndex < xAdjEnd; wIndex++) {
                            int w = xAdj[wIndex];

                            if(visited[w] || w == v) {
                                // This thread has already visited, or it's an edge back to parent, ignore!
                                continue;
                            } else if(isXArt && Low[w] > -1) {
                                // Already know it's an articulation point, ignore!
                                continue;
                            } else if(L[w] < L[u]) {
                                // Ref. lines 11 & 12, Alg. 8 (BFS-LV)
                                // Jump out of this BFS and set everything on the stack to unvisited
                                goto nextOut;
                            } else {
                                // Can keep searching!
                                // Ref. lines 14--18, Alg. 8 (BFS-LV)

                                // Add to the next queue so we process it after everything in the current queue
                                queue[queueEnd++] = w;
                                Vu[stackEnd++] = w;
                                visited[w] = 1;
                                if(w < vidLow) {
                                    vidLow = w;
                                }
                            }
                        }
                    }

                    // Ref. line 14 - 21, Alg. 8 (BFS-LV)
                    Art[v] = 1;
                    numBiconnectedComponents++;
                    isCurrArt = 1;
                    // For each of the unique vertices encountered
                    for(int j = 0; j < stackEnd; j++) {
                        int w = Vu[j];
                        Low[w] = vidLow;
                        Par[w] = v;
                        visited[w] = 0;
                    }

                    nextOut:            // Label to get out of nested loop above
                    if(!isCurrArt) {
                        // We'll get here from the jump above when L[w] < L[u] - set everything on stack to unvisited
                        for(int j = 0; j < stackEnd; j++) {
                            int s = Vu[j];
                            visited[s] = 0;
                        }
                    }
                }
            }
        }

        free(visited);
        free(Vu);
        free(queue);
        free(nextQueue);
    }
}

void findLowValues(Graph *graph, int root, int *P, int *L, int *Par, int *Low, int **LQ, int *LQCounts, int *Art) {
    int numThreads = omp_get_max_threads();

    // List of vertices to visit - appended to during the search
    int *toVisit = malloc(sizeof(int) * graph->numVertices);
    int back = 0;
    int *toVisitThreads = malloc(sizeof(int) * (numThreads * graph->numVertices));
    int *threadOffsets = malloc(sizeof(int) * numThreads);
    int *stack = (int *) malloc(sizeof(int) * (graph->numVertices * 2));
    int *lows = (int *) malloc(sizeof(int) * numThreads);
    int stackBack = 0;

    int *visited = (int *) malloc(sizeof(int) * graph->numVertices);
    for(int i = 0; i < graph->numVertices; i++) {
        visited[i] = 0;
    }

    int levelSize = LQCounts[1];
    visited[root] = 1;

    // Direction optimising BFS vars
    double alpha = 14.0;
    double beta = 24.0;
    int useBottomUp = 0;
    int nf = 0;
    int localnf;

    for(int l = 0; l < levelSize; l++) {
        int globalLow = graph->numVertices;
        int v = LQ[1][l];
        if(Low[v] == -1) {
            toVisit[0] = v;
            back = 1;
            Par[v] = 0;
            stack[0] = v;
            stackBack = 1;

            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                int *toVisitThread = &toVisitThreads[tid * graph->numVertices];
                int threadOffset = 0;
                int prevLevel;
                int vert, adjEnd, u;
                int *adjVertices;
                int threadLow = graph->numVertices;

                while(back) {
                    threadOffset = 0;

                    if(!useBottomUp) {
                        #pragma omp for schedule(static) reduction(+:localnf)
                        for(int i = 0; i < back; i++) {
                            vert = toVisit[i];
                            adjEnd = outDegree(graph, vert);
                            adjVertices = adjacentVertices(graph, vert);
                            for(int j = 0; j < adjEnd; j++) {
                                u = adjVertices[j];

                                if(!visited[u] && Low[u] < 0) {
                                    visited[u] = 1;
                                    Par[u] = 0;
                                    toVisitThread[threadOffset++] = u;
                                    if(u < threadLow) {
                                        threadLow = u;
                                    }
                                }
                            }
                        }
                    } else {
                        #pragma omp for schedule(static) reduction(+:localnf)
                        for(int i = 0; i < graph->numVertices; i++) {
                            vert = i;
                            if(!visited[vert] && Low[vert] < 0) {
                                adjEnd = outDegree(graph, vert);
                                adjVertices = adjacentVertices(graph, vert);
                                for(int j = 0; j < adjEnd; j++) {
                                    u = adjVertices[j];
                                    if(visited[j]) {
                                        visited[vert] = 1;
                                        Par[vert] = 0;
                                        toVisitThread[threadOffset++] = vert;
                                        if(vert < threadLow) {
                                            threadLow = vert;
                                        }
                                        localnf++;
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    threadOffsets[tid] = threadOffset;
                    lows[tid] = threadLow;

                    #pragma omp barrier

                    #pragma omp single
                    {
                        nf += localnf;

                        if(useBottomUp) {
                            double mf = (double) localnf * graph->avgOutDegree;
                            double mu = (double) (graph->numVertices - nf) * graph->avgOutDegree;

                            // Switch if heuristic is true, and if we haven't switched yet
                            if(mf > (mu / alpha) && mu > 0) {
                                useBottomUp = 1;
                            }
                        } else {
                            if(nf < ((double) graph->numVertices / beta)) {
                                useBottomUp = 0;
                            }
                        }

                        back = 0;
                        for(int i = 0; i < numThreads; i++) {
                            if(lows[i] < globalLow) {
                                globalLow = lows[i];
                            }

                            if(threadOffsets[i]) {
                                int offset = i * graph->numVertices;
                                for(int j = 0; j < threadOffsets[i]; j++) {
                                    toVisit[back++] = toVisitThreads[offset + j];
                                }
                            }
                        }

                        for(int i = 0; i < back; i++) {
                            stack[i + stackBack] = toVisit[i];
                        }
                        stackBack += back;
                    }
                }

                #pragma omp for
                for(int i = 0; i < stackBack; i++) {
                    int v = stack[i];
                    Low[v] = globalLow;
                    visited[v] = 0;
                }
            }
        }
    }

    // Deal with root
    int end = outDegree(graph, root);
    int *adjVertices = adjacentVertices(graph, root);
    for(int i = 0; i < end; i++) {
        int currLow = Low[adjVertices[i]];
        for(int j = i + 1; j < end; j++) {
            if(Low[adjVertices[j]] == currLow) {
                Low[root] = currLow;
                break;
            }
        }
        if(Low[root] >= 0) {
            break;
        }
    }

    if(Low[root] < 0) {
        Art[root] = 1;
        Low[root] = root;
    }

    free(visited);
    free(toVisit);
    free(toVisitThreads);
    free(stack);
    free(threadOffsets);
    free(lows);
}

/**
 *  Labels the edges of the given graph using Par and Low values.
 *  Edge labels are stored in the graph itself as values for each edge in the CSR representation
 */
void labelEdges(Graph *graph, int *Par, int *Low) {
    #pragma omp parallel for schedule(static)
    for(int u = 0; u < graph->numVertices; u++) {
        int uLow = Low[u];
        int uPar = Par[u];

        int start = edgeOffset(graph, u);
        int end = edgeOffset(graph, u + 1);
        for(int i = start; i < end; i++) {
            int v = graph->edgeArray[i];
            if(Low[v] == uLow || uPar == v) {
                graph->edgeLabels[i] = uLow;
            } else {
                graph->edgeLabels[i] = Low[v];
            }
        }
    }
}

/**
 *  Discovers the biconnected components in the given graph
 */
void bicc(Graph *graph) {
    double elt, elt2;
    elt = timer();
    elt2 = timer();

    // Get memory for the needed arrays
    int numVertices = graph->numVertices;
    int arrSize = sizeof(int) * numVertices;
    int *P = (int *) malloc(arrSize);
    int *L = (int *) malloc(arrSize);
    int *Art = (int *) malloc(arrSize);
    int *Low = (int *) malloc(arrSize);
    int *Par = (int *) malloc(arrSize);
    int *visited = (int *) malloc(arrSize);
    int **LQ = (int **) malloc(sizeof(int *) * MAX_LEVELS); 
    int *LQCounts = (int *) malloc(sizeof(int) * MAX_LEVELS); 
    int levels = 0;

    // Initialise the arrays
    #pragma omp parallel for
    for(int i = 0; i < numVertices; i++) {
        Art[i] = 0;
        Low[i] = -1;
        Par[i] = -1;
        L[i] = -1;
        P[i] = -1;
        visited[i] = 0;
    }

    elt = timer() - elt;
    printf("init time: %9.6lf\n", elt); 
    elt = timer();

    // Pick the vertex with the highest out degree as root
    int root = graph->maxDegreeVertex;
    
    // Initial BFS to get P, L and LQ
    bfs(graph, root, &levels, P, L, LQ, LQCounts, visited);

    elt = timer() - elt;
    printf("BFS time: %9.6lf\n", elt);
    elt = timer();

    // BFS to discover articulation points
    findArticulationPoints(graph, P, L, Par, Low, LQ, LQCounts, Art, levels);

    elt = timer() - elt;
    printf("articulation time: %9.6lf\n", elt);
    elt = timer();

    findLowValues(graph, root, P, L, Par, Low, LQ, LQCounts, Art);

    elt = timer() - elt;
    printf("low values time: %9.6lf\n", elt);
    elt = timer();

    // Label the edges
    labelEdges(graph, Par, Low);

    elt = timer() - elt;
    printf("label time: %9.6lf\n", elt);

    // for(int u = 0; u < numVertices; u++) {
    //     int start = edgeOffset(graph, u);
    //     int end = edgeOffset(graph, u + 1);

    //     for(int i = start; i < end; i++) {
    //         int v = graph->edgeArray[i];
    //         int label = graph->edgeLabels[i];
    //         printf("(%d, %d) = %d\n", u, v, label);
    //     }
    // }

    elt2 = timer() - elt2;
    printf("BCC_BFS %d total: %9.6lf\n", omp_get_max_threads(), elt2);

    printf("num bcc: %d\n", numBiconnectedComponents);
    // printf("max bcc: %d\n", maxBiconnectedComponent);

    free(P);
    free(L);
    free(Par);
    free(Low);
    free(Art);
    for(int i = 0; i < levels; i++) {
        free(LQ[i]);
    }
    free(LQ);
    free(LQCounts);
}

int main(int argc, char *argv[]) {
    // Check that all params are given
    if(argc != 2) {                             
        fprintf(stderr, "Usage: ./bic graphFile\n");
        return 1;
    }

    // Parse the input file and generate graph
    char* filename = argv[1];
    Graph *graph = parseGraphFile(filename);

    // Find biconnected components
    bicc(graph);

    destroyGraph(graph);
}

