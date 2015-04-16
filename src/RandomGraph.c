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




#include "RandomGraph.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAX_NAME_SIZE 3000
#define EPSILON 0.001
#define MAX_ITER 50





/*****************************************************************************************************************/
/*Simulation*/
/*****************************************************************************************************************/

/* return  a random type of event wrt the rates*/
int getClass(double *prob, int sizeGraph){
	int i;
	double uni = UNIF_RAND, cumul=prob[0];
	for(i=0; i<(sizeGraph-1) && uni>cumul; i++)
		cumul += prob[i+1];
	return i;
}




/* return  a random type of event wrt the rates*/
int getEdge(double prob) {
	return (UNIF_RAND<prob);
}

int *getPresent(int sizeGraph, double prob) {
	int *present, i;
	present = (int*) malloc(sizeGraph*sizeof(int));
	for(i=0; i<sizeGraph; i++)
		present[i] = (UNIF_RAND<=prob);
	return present;
}

/*simulate a random multi graph following a unique model*/
TypeMultiGraph *getRandomMultiGraphER(double p, int sizeTable, TypeSizeG sizeGraph) {
	TypeMultiGraph *g;
	TypeSizeG n, m;
	int t;
	g = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
	g->sizeGraph = sizeGraph;
	g->sizeTable = sizeTable;
	g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
	for(n=0; n<g->sizeGraph; n++) {
		char tmp[100];
		sprintf(tmp, "%d", n+1);
		g->name[n] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
		strcpy(g->name[n], tmp);
	}
	g->edge = (TypeEdgeG***) malloc(g->sizeTable*sizeof(TypeEdgeG**));
	g->present = (int**) malloc(g->sizeTable*sizeof(int*));
	for(t=0; t<g->sizeTable; t++) {
		g->present[t] = (int*) malloc(g->sizeGraph*sizeof(int));
		g->edge[t] = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
		for(n=0; n<g->sizeGraph; n++) {
			g->present[t][n] = 1;
			g->edge[t][n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
			g->edge[t][n][n] = 0;
			for(m=0; m<n; m++) {
				g->edge[t][n][m] = getEdge(p);
				g->edge[t][m][n] = g->edge[t][n][m];
			}
		}
	}
	return g;
}

/*simulate a random multi graph following a unique model*/
TypeMultiGraph *getRandomMultiGraphUniSpecial(TypeERMG *model, int sizeTable, TypePartition *part) {
	TypeMultiGraph *g;
	TypeSizeG n, m;
	char tmp[100];
	int *class, t, bad = 0, *test, c;
	class = part->atom;

	g = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
	g->sizeGraph = part->sizeItem;
	g->sizeTable = sizeTable;
	g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
	for(n=0; n<g->sizeGraph; n++) {
		char tmp[100];
		sprintf(tmp, "%d", n+1);
		g->name[n] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
		strcpy(g->name[n], tmp);
	}
	g->edge = (TypeEdgeG***) malloc(g->sizeTable*sizeof(TypeEdgeG**));
	g->present = (int**) malloc(g->sizeTable*sizeof(int*));
	for(t=0; t<g->sizeTable; t++) {
		g->present[t] = (int*) malloc(g->sizeGraph*sizeof(int));
		g->edge[t] = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
		for(n=0; n<g->sizeGraph; n++) {
			g->present[t][n] = 1;
			g->edge[t][n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
			g->edge[t][n][n] = 0;
			for(m=0; m<n; m++) {
				g->edge[t][n][m] = getEdge(model->pi[class[n]][class[m]]);
				g->edge[t][m][n] = g->edge[t][n][m];
			}
		}
	}
	return g;
}

/*simulate a random multi graph following a unique model*/
TypeMultiGraph *getRandomMultiGraphUni(TypeERMG *model, int sizeTable, TypeSizeG sizeGraph, TypePartition *part) {
	TypeMultiGraph *g;
	TypeSizeG n, m;
	char tmp[100];
	int *class, t, bad = 0, *test, c;
	class = (int*) malloc(sizeGraph*sizeof(int));
	test = (int*) malloc(model->sizeAtom*sizeof(int));
	do {
		bad = 0;
		for(c=0; c<model->sizeAtom; c++)
			test[c] = 0;
		for(n=0; n<sizeGraph; n++) {
			class[n] = getClass(model->alpha, model->sizeAtom);
			test[class[n]]++;
		}
		for(c=0; c<model->sizeAtom && !bad; c++)
			bad = (test[c] < 2);	
	} while(bad);
	free((void*) test);
	g = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
	g->sizeGraph = sizeGraph;
	g->sizeTable = sizeTable;
	g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
	for(n=0; n<g->sizeGraph; n++) {
		char tmp[100];
		sprintf(tmp, "%d", n+1);
		g->name[n] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
		strcpy(g->name[n], tmp);
	}
	g->edge = (TypeEdgeG***) malloc(g->sizeTable*sizeof(TypeEdgeG**));
	g->present = (int**) malloc(g->sizeTable*sizeof(int*));
	for(t=0; t<g->sizeTable; t++) {
		g->present[t] = (int*) malloc(g->sizeGraph*sizeof(int));
		g->edge[t] = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
		for(n=0; n<g->sizeGraph; n++) {
			g->present[t][n] = 1;
			g->edge[t][n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
			g->edge[t][n][n] = 0;
			for(m=0; m<n; m++) {
				g->edge[t][n][m] = getEdge(model->pi[class[n]][class[m]]);
				g->edge[t][m][n] = g->edge[t][n][m];
			}
		}
	}
//	monfree((void*) class);
	part->sizeAtom = model->sizeAtom;
	part->sizeItem = g->sizeGraph;
	part->atom = class;
	return g;
}

/*simulate a random multi graph following model*/
TypeMultiGraph *getRandomMultiGraph(TypeMultiERMG *model, TypeSizeG sizeGraph, TypePartition *part) {
	TypeMultiGraph *g;
	TypeSizeG n, m;
	char tmp[100];
	int t;

	g = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
	g->sizeGraph = sizeGraph;
	g->sizeTable = model->sizeTable;
	g->name = NULL;
	g->edge = (TypeEdgeG***) malloc(g->sizeTable*sizeof(TypeEdgeG**));
	g->present = (int**) malloc(g->sizeTable*sizeof(int*));
	for(t=0; t<g->sizeTable; t++) {
		g->present[t] = (int*) malloc(g->sizeGraph*sizeof(int));
		g->edge[t] = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
		for(n=0; n<g->sizeGraph; n++) {
			g->present[t][n] = 1;
			g->edge[t][n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
			g->edge[t][n][n] = 0;
			for(m=0; m<n; m++) {
				g->edge[t][n][m] = getEdge(model->pi[t][part->atom[n]][part->atom[m]]);
				g->edge[t][m][n] = g->edge[t][n][m];
			}
		}
	}
	return g;
}


/*simulate a random multi graph following model*/
TypeMultiGraph *getRandomMultiGraphMissing(TypeMultiERMG *model, TypeSizeG sizeGraph, TypePartition *part, double ppres) {
	TypeMultiGraph *g;
	TypeSizeG n, m;
	char tmp[100];
	int t;

	g = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
	g->sizeGraph = sizeGraph;
	g->sizeTable = model->sizeTable;
	g->name = NULL;
	g->edge = (TypeEdgeG***) malloc(g->sizeTable*sizeof(TypeEdgeG**));
	g->present = (int**) malloc(g->sizeTable*sizeof(int*));
	for(t=0; t<g->sizeTable; t++) {
		g->present[t] = getPresent(g->sizeGraph, ppres);
		g->edge[t] = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
		for(n=0; n<g->sizeGraph; n++) {
			g->edge[t][n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
			g->edge[t][n][n] = 0;
			if(g->present[t][n]) {
				for(m=0; m<n; m++) {
					if(g->present[t][m]) {
						g->edge[t][n][m] = getEdge(model->pi[t][part->atom[n]][part->atom[m]]);
						g->edge[t][m][n] = g->edge[t][n][m];
					} else {
						g->edge[t][n][m] = 0;
						g->edge[t][m][n] = 0;
					}
				}
			} else {
				for(m=0; m<n; m++) {
					g->edge[t][n][m] = 0;
					g->edge[t][m][n] = 0;
				}
			}
		}
	}
	return g;
}

int *getRandomClass(TypeMultiERMG *model, int sizeGraph) {
	TypeSizeG n, m;
	int *class, t, bad = 0, *test, c;
	class = (int*) malloc(sizeGraph*sizeof(int));
	test = (int*) malloc(model->sizeAtom*sizeof(int));
	do {
		bad = 0;
		for(c=0; c<model->sizeAtom; c++)
			test[c] = 0;
		for(n=0; n<sizeGraph; n++) {
			class[n] = getClass(model->alpha, model->sizeAtom);
			test[class[n]]++;
		}
		for(c=0; c<model->sizeAtom && !bad; c++)
			bad = (test[c] < 2);	
	} while(bad);
	free((void*) test);
	return class;
}

int *getBalancedClass(int sizeAtom, int sizeGraph) {
	int *class, t, bad = 0, *test, c;
	TypeSizeG n, m;
	class = (int*) malloc(sizeGraph*sizeof(int));
	for(n=0; n<sizeGraph; n++) {
		class[n] = (n*sizeAtom)/sizeGraph;
	}
	return class;
}

void setRandomPartition(TypeMultiERMG *model, int sizeGraph, TypePartition *part) {
	part->sizeAtom = model->sizeAtom;
	part->sizeItem = sizeGraph;
	part->atom = getRandomClass(model, sizeGraph);
}

void setBalancedPartition(int sizeAtom, int sizeGraph, TypePartition *part) {
	part->sizeAtom = sizeAtom;
	part->sizeItem = sizeGraph;
	part->atom = getBalancedClass(sizeAtom, sizeGraph);
}

/*simulate a random multi graph following model*/
TypeMultiGraph *getMultiGraphWeighted(TypeMultiERMG *model, TypeSizeG sizeGraph, TypePartition *part) {
	TypeMultiGraph *g;
	TypeSizeG n, m;
	char tmp[100];
	int *class, t, bad = 0, *test, c;
	class = part->atom;
	g = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
	g->sizeGraph = sizeGraph;
	g->sizeTable = model->sizeTable;
	g->name = NULL;
	g->edge = (TypeEdgeG***) malloc(g->sizeTable*sizeof(TypeEdgeG**));
	g->present = (int**) malloc(g->sizeTable*sizeof(int*));
	for(t=0; t<g->sizeTable; t++) {
		g->present[t] = (int*) malloc(g->sizeGraph*sizeof(int));
		g->edge[t] = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
		for(n=0; n<g->sizeGraph; n++) {
			g->present[t][n] = 1;
			g->edge[t][n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
			g->edge[t][n][n] = 0;
			for(m=0; m<n; m++) {
				g->edge[t][n][m] = model->pi[t][class[n]][class[m]];
				g->edge[t][m][n] = g->edge[t][n][m];
			}
		}
	}
	return g;
}

/*simulate a random graph following model*/
TypeGraph *getRandomGraph(TypeERMG *model, TypeSizeG sizeGraph, TypePartition *part) {
	TypeGraph *g;
	TypeSizeG n, m;
	char tmp[100];

	if(part->atom == NULL) {
		int bad = 0, *test, c;
		part->atom = (int*) malloc(sizeGraph*sizeof(int));
		test = (int*) malloc(model->sizeAtom*sizeof(int));
		for(c=0; c<model->sizeAtom; c++)
			test[c] = 0;
		do {
			for(n=0; n<sizeGraph; n++) {
				part->atom[n] = getClass(model->alpha, model->sizeAtom);
				test[part->atom[n]]++;
			}
			for(c=0; c<model->sizeAtom && !bad; c++)
				bad = (test[c] == 0);	
		} while(bad);
		free((void*) test);
		part->sizeAtom = model->sizeAtom;
		part->sizeItem = sizeGraph;
	} else {
		if(part->sizeAtom != model->sizeAtom || part->sizeItem != sizeGraph) {
			printf("Bad partition in random graph generation\n");
			exit(1);
		}
	}
	g = (TypeGraph*) malloc(sizeof(TypeGraph));
	g->sizeGraph = sizeGraph;
	g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
	g->edge = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
	for(n=0; n<g->sizeGraph; n++) {
		char tmp[100];
		sprintf(tmp, "%d", n+1);
		g->name[n] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
		strcpy(g->name[n], tmp);
		g->edge[n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
		g->edge[n][n] = 0;
		for(m=0; m<n; m++) {
			g->edge[n][m] = getEdge(model->pi[part->atom[n]][part->atom[m]]);
			g->edge[m][n] = g->edge[n][m];
		}
	}
//	monfree((void*) class);
	return g;
}
