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




#ifndef RandomGraphF
#define RandomGraphF

#include <stdlib.h>
#include <stdio.h>
#include "Graph.h"
#include "ERMG.h"


/* return  a random type of event wrt the rates*/
int getClass(double *prob, int size);
/* return  a random type of event wrt the rates*/
int getEdge(double prob);
/*simulate a random graph following model*/
TypeGraph *getRandomGraph(TypeERMG *model, TypeSizeG size, TypePartition *part);
/*simulate a random multi graph following a unique model*/
TypeMultiGraph *getRandomMultiGraphUni(TypeERMG *model, int sizeTable, TypeSizeG sizeGraph, TypePartition *part);
TypeMultiGraph *getRandomMultiGraphUniSpecial(TypeERMG *model, int sizeTable, TypePartition *part);
/*simulate a random multi graph following model*/
TypeMultiGraph *getRandomMultiGraph(TypeMultiERMG *model, TypeSizeG sizeGraph, TypePartition *part);
/*simulate a random multi graph following a unique model*/
TypeMultiGraph *getRandomMultiGraphER(double p, int sizeTable, TypeSizeG sizeGraph);
/*create a weighted multi graph following the edge probabilities of model*/
TypeMultiGraph *getMultiGraphWeighted(TypeMultiERMG *model, TypeSizeG sizeGraph, TypePartition *part);
int *getRandomClass(TypeMultiERMG *model, int sizeGraph);
int *getBalancedClass(int sizeAtom, int sizeGraph);
void setRandomPartition(TypeMultiERMG *model, int sizeGraph, TypePartition *part);
void setBalancedPartition(int sizeAtom, int sizeGraph, TypePartition *part);
int *getPresent(int sizeGraph, double prob);
/*simulate a random multi graph following model*/
TypeMultiGraph *getRandomMultiGraphMissing(TypeMultiERMG *model, TypeSizeG sizeGraph, TypePartition *part, double ppres);
#endif
