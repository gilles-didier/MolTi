/*
    'MolTi' and 'molti-console' detects communities from multiplex networks / 'bonf' computes q-values of annotations enrichment of communities / 'test' simulates random multiplex to test community detection approaches
    Copyright (C) 2015  Gilles DIDIER

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#ifndef GraphF
#define GraphF

#include <stdlib.h>
#include <stdio.h>
#include "Utils.h"
#include "Partition.h"


#define INC_SIZE 100
#define MAX_CURRENT 1000000
#define UNKNOWN -1

#ifdef __cplusplus
extern "C" {
#endif


typedef enum DISPLAY_G {
    display_none=0,
    display_name,
    display_index,
    display_both,
    display_time_none,
    display_time_name,
    display_time_index,
    display_time_both
} TypeDisplayG;

typedef float TypeEdgeG;
typedef int TypeSizeG;


typedef struct GRAPH {
    char **name;
    TypeEdgeG **edge;
    TypeSizeG sizeGraph;
} TypeGraph;

typedef struct MULTI_GRAPH {
    char **name;
    TypeEdgeG ***edge;
    TypeSizeG sizeGraph;
    int sizeTable, **present;
} TypeMultiGraph;

TypeMultiGraph *toMultiGraph(TypeGraph **graph, int sizeTable);
void fprintGraphTable(FILE *f, TypeGraph *g);
void fprintGraph(FILE *f, TypeGraph *Graph);
TypeGraph *cpyGraph(TypeGraph *graph);
/*read graph as a matrix*/
TypeGraph *readMatrixGraph(FILE *f);
/*read graph as list of edges*/
TypeGraph *readGraph(FILE *f);
void freeGraph(TypeGraph *graph);
/*free multi graph in standard format*/
void freeMultiGraph(TypeMultiGraph *g);
/*free multi graph in standard format*/
void freeGraph(TypeGraph *g);
TypeGraph *newGraph();
/*print ith graph in standard format*/
void fprintMultiGraph(FILE *f, int i, TypeMultiGraph *g);
/*print ith graph in standard format*/
void fprintMultiGraphDebug(FILE *f, TypeMultiGraph *g);
TypeGraph *sumMultiGraph(TypeMultiGraph *multi);
TypeGraph *interMultiGraph(TypeMultiGraph *multi);
TypeGraph *unionMultiGraph(TypeMultiGraph *multi);
TypeMultiGraph *sumMultiGraphMulti(TypeMultiGraph *multi);
TypeMultiGraph *interMultiGraphMulti(TypeMultiGraph *multi);
TypeMultiGraph *unionMultiGraphMulti(TypeMultiGraph *multi);

/*print ith graph in Octave matrix format*/	
void fprintMultiGraphOctaveMatrix(FILE *f, int i, TypeMultiGraph *g, char *id);
void fprintGraphInci(FILE *f, TypeGraph *g);
void fprintMultiGraphTable(FILE *f, int t, TypeMultiGraph *g);
int countEdge(TypeMultiGraph *g);
TypeMultiGraph *cpyMultiGraph(TypeMultiGraph *h);
double getModularity(TypeGraph *g, TypePartition *part);
double getModularityMulti(TypeMultiGraph *g, TypePartition *part);
double getERModularity(TypeGraph *g, TypePartition *part);
double getERModularityMulti(TypeMultiGraph *g, TypePartition *part);
double getERFullModularityMulti(TypeMultiGraph *g);
double getFullModularityMulti(TypeMultiGraph *g);
double getERFullModularityMultiOne(int t, TypeMultiGraph *g);
void fixEdgeMultiGraph(TypeMultiGraph *g);
void fixEdgeGraph(TypeGraph *g);
double getTestMulti(TypeMultiGraph *g);
/*reorder graph vertices wrt their names*/
void orderGraphName(TypeGraph *g);
double getGraphDensity(TypeGraph *g);

#ifdef __cplusplus
}
#endif

#endif
