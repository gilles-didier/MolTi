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




#include <assert.h>
#include <math.h>
#include "EdgesComposition.h"
#include "Utils.h"




int **getEdgesNumbers(TypePartition *part, TypeMultiGraph *graph) {
	int **eff, i, j, c, t;
	TypePartitionList *list;
	
	list = getPartitionList(part);
	if(part->sizeItem != graph->sizeGraph) {
		printf("Partition %d and Graph %d\n", (int) part->sizeItem, (int) graph->sizeGraph);
		exit(1);
	}
	eff = (int**) malloc(part->sizeAtom*sizeof(int*));
	for(c=0; c<part->sizeAtom; c++) {
		eff[c] = (int*) malloc(graph->sizeTable*sizeof(int));
		for(t=0;t<graph->sizeTable; t++)
			eff[c][t] = 0;
	}
	for(c=0; c<part->sizeAtom; c++) {
		int i,j, t;
		for(i=list->start[c]; list->next[i]>=0; i=list->next[i]) {
			for(j=list->next[i]; j>=0; j=list->next[j]) {
				for(t=0; t<graph->sizeTable; t++) {
					eff[c][t] += graph->edge[t][i][j];
				}
			}
		}
	}
	freePartitionList(list);
	return eff;
}
