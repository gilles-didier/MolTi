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




#ifndef LouvainF
#define LouvainF

#include <stdlib.h>
#include <stdio.h>
#include "Utils.h"
#include "Graph.h"
#include "Partition.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum TYPE_PARTITION_METHOD {
    LouvainType=0
} TypePartitionMethod;

typedef struct LOUVAIN_PARAM {
    double  m, *k;
} TypeLouvainParam;

typedef struct LOUVAIN_PARAM_MULTI {
    int size, sizeBuf;
    double alpha;
    TypeLouvainParam *table;
} TypeLouvainParamMulti;

typedef struct EDGE_ELEMENT {
    int neighbour, next;
} TypeElementEdge;

typedef struct PROTO_CLASS {
    int size, first, next, *prec;
} TypeProtoClass;

typedef struct ELEMENT_CLASS {
    int first, atom, next, size, *prec;
} TypeElementClass;

typedef double (*TypeVariation)();

typedef void (*TypeUpdateParam)();

typedef struct FUNCTION_STATE {
    double (*variation)();
    void (*init)(), (*update)(), (*transfert)();
} TypeFunctionState;

typedef struct COMMUNITY_STATE {
    int sizeElement, sizeProto, first, trash, *elementInit, sizeInit, nedge, nedgeTot;
    TypeProtoClass *community;
    TypeElementClass *element;
    TypeElementEdge *edgeList;
    TypeMultiGraph *graph, *gsave;
} TypeCommunityState;

typedef struct LOUVAIN_GLOBAL {
    TypeCommunityState *state;
    void *param;
    double (*variation)(int, int, TypeCommunityState*, void*);
    void (*init)(TypeMultiGraph*, void*, TypeCommunityState**, void**), (*update)(TypeCommunityState**, void**), (*transfert)(int, int, TypeCommunityState*, void*), (*freeParam)(void*);
} TypeLouvainGlobal;

/****************************************************************************/
/* General functions on state*/
/****************************************************************************/

void freeCommunityState(TypeCommunityState* state);
TypePartition communityStateToPartition(TypeCommunityState *state);
TypeMultiGraph *getUpdatedGraph(TypeCommunityState *state, int *index);
TypeCommunityState *initCommunityState(TypeMultiGraph *graph);
TypeCommunityState *getUpdatedCommunityState(TypeCommunityState *state);
void printCommunityState(TypeCommunityState *state);
void printCommunityStateBis(TypeCommunityState *state);
int countCommunity(TypeCommunityState *state);
/*test if e is the single element of protoatom i*/
int isSingle(int e, TypeCommunityState *state);
/*transfer element e from protoatom i to protoatom j  !! e must be in protoatom i !!*/
void transfertStandard(int e, int j, TypeCommunityState *state, void *param);

/****************************************************************************/
/* Louvain functions */
/****************************************************************************/

double variationLouvain(int e, int j, TypeCommunityState *state, void *param);
void computeParamLouvain(TypeCommunityState *state, TypeLouvainParamMulti *param);
void *getParamLouvain(TypeCommunityState *state, double alpha);
void freeParamLouvain(void *param);
void updateLouvain(TypeCommunityState **state, void **param);
void initLouvain(TypeMultiGraph *graph, void *info, TypeCommunityState **stateP, void **paramP);

/****************************************************************************/
/* Global functions */
/****************************************************************************/

int iterateGlobal(TypeLouvainGlobal *global);
TypePartition getPartition(TypeMultiGraph *graph, TypePartitionMethod type, void *info);

#ifdef __cplusplus
}
#endif

#endif
