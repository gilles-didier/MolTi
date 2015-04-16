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




#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Graph.h"


#define INFTY 9.9E99;
#define BASIC_TMP_SIZE 1000
#define BASIC_INC_BUFFER 20
#define INC_SIZE_BUF 1000
#define MAX_NAME_SIZE 3000

double getModularity(TypeGraph *g, TypePartition *part) {
	double res = 0., *k, m;
	TypePartitionList *pl;
	int c, i, j;
	
	m = 0.;
	k = (double*) malloc(g->sizeGraph*sizeof(double));
	for(i=0; i<g->sizeGraph; i++) {
		k[i] = 0.;
		for(j=0; j<g->sizeGraph; j++)
			k[i] += g->edge[i][j];
		m += k[i];
	}
	m /= 2.;
	pl = getPartitionList(part);
	for(c=0; c<part->sizeAtom; c++) {
		int e, f;
		for(e=pl->start[c]; e>=0; e=pl->next[e])
			for(f=pl->start[c]; f>=0; f=pl->next[f])
				res += ((double) g->edge[e][f]-k[e]*k[f]/(2.*m))/(2.*m);
	}
	free((void*)k);
	freePartitionList(pl);
	return res;
}


double getGraphDensity(TypeGraph *g) {
	double res = 0.,m = 0;
	int i, j;

	for(i=0; i<g->sizeGraph; i++)
		for(j=0; j<i; j++) {
			res += g->edge[i][j];
			m++;
		}
	res /= m;
	return res;
}

double getModularityMulti(TypeMultiGraph *g, TypePartition *part) {
	double res = 0., *k, m;
	TypePartitionList *pl;
	int t, c, i, j;
	
	k = (double*) malloc(g->sizeGraph*sizeof(double));
	for(t=0; t<g->sizeTable; t++) {
		m = 0.;
		for(i=0; i<g->sizeGraph; i++) {
			k[i] = 0.;
			for(j=0; j<g->sizeGraph; j++)
				k[i] += g->edge[t][i][j];
			m += k[i];
		}
		m /= 2.;
		pl = getPartitionList(part);
		for(c=0; c<part->sizeAtom; c++) {
			int e, f;
			for(e=pl->start[c]; e>=0; e=pl->next[e])
				for(f=pl->next[e]; f>=0; f=pl->next[f])
					res += ((double) g->edge[t][e][f]-k[e]*k[f]/(2.*m))/(2.*m);
		}
	}
	free((void*)k);
	freePartitionList(pl);
	return res;
}

double getFullModularityMulti(TypeMultiGraph *g) {
	double res = 0., *k, m, pi;
	int t, c, i, j;
	
	k = (double*) malloc(g->sizeGraph*sizeof(double));
	for(t=0; t<g->sizeTable; t++) {
		double pi;
		m = 0.;
		for(i=0; i<g->sizeGraph; i++) {
			k[i] = 0.;
			for(j=0; j<g->sizeGraph; j++)
				k[i] += g->edge[t][i][j];
			m += k[i];
		}
		m /= 2.;
		pi = ((double)m)/(((double)g->sizeGraph)*((double)g->sizeGraph-1));
		for(i=0; i<g->sizeGraph; i++)
			for(j=0; j<i; j++)
				res += ((double) g->edge[t][i][j]-k[i]*k[j]/(2.*m))/(2.*m);
	}
	free((void*)k);
	return res;
}

double getTestMulti(TypeMultiGraph *g) {
	double res = 0., *k, m, pi, N;
	int t, c, i, j;
	
	k = (double*) malloc(g->sizeGraph*sizeof(double));
	for(t=0; t<g->sizeTable; t++) {
		double pi;
		m = 0.;
		for(i=0; i<g->sizeGraph; i++) {
			k[i] = 0.;
			for(j=0; j<g->sizeGraph; j++)
				k[i] += g->edge[t][i][j];
			m += k[i];
		}
		m /= 2.;
		pi = ((double)m)/(((double)g->sizeGraph)*((double)g->sizeGraph-1));
		N = (double)g->sizeGraph;
		for(i=0; i<g->sizeGraph; i++)
			for(j=0; j<i; j++)
				res += (((N-2.)*(N-3.)/2.)*((double)k[i])*((double)k[i])+(N-1.)*((double)m-k[i]-k[j]+1)*(N-1.-k[i]-k[j]))/(((N)*(N-1.)/2.)*((double)k[i])*((double)k[i])+(N-1.)*((double)m)*(N-1.-k[i]-k[j]));
	}
	free((void*)k);
	return res;
}

double expERExp(double pi, int ki, int kj, int N) {
	return ((1.-pi)*((double)ki)*((double)kj))/(((double)ki)*((double)kj)+pi*(((double)N)-1.)*(((double)N)-1.-((double)ki)-((double)kj)));
}

double getERModularity(TypeGraph *g, TypePartition *part) {
	double res = 0., *k, m, pi;
	TypePartitionList *pl;
	int c, i, j;
	
	m = 0.;
	k = (double*) malloc(g->sizeGraph*sizeof(double));
	for(i=0; i<g->sizeGraph; i++) {
		k[i] = 0.;
		for(j=0; j<g->sizeGraph; j++)
			k[i] += g->edge[i][j];
		m += k[i];
	}
	m /= 2.;
	pi = ((double)m)/(((double)g->sizeGraph)*((double)g->sizeGraph-1));
	pl = getPartitionList(part);
	for(c=0; c<part->sizeAtom; c++) {
		int e, f;
		for(e=pl->start[c]; e>=0; e=pl->next[e])
			for(f=pl->next[e]; f>=0; f=pl->next[f])
				res += ((double) g->edge[e][f])-expERExp(pi,k[e],k[f],g->sizeGraph);
	}
	free((void*)k);
	freePartitionList(pl);
	return res;
}

void orderGraphName(TypeGraph *g) {
	size_t *index;
	int i, j;
	TypeEdgeG **new;
	if(g->name == NULL)
		return;
	new = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
	for(i=0; i<g->sizeGraph; i++)
		new[i] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
	index = qsortindex(g->name, g->sizeGraph, sizeof(char*), compareString);
	for(i=0; i<g->sizeGraph; i++) {
		for(j=0; j<g->sizeGraph; j++) {
			new[index[i]][index[j]] = g->edge[i][j];
		}
	}
	for(i=0; i<g->sizeGraph; i++)
		free((void*)g->edge[i]);
	free((void*)g->edge);
	free((void*)index);
	g->edge = new;
}
		


double getERModularityMulti(TypeMultiGraph *g, TypePartition *part) {
	double res = 0., *k, m, pi;
	TypePartitionList *pl;
	int t, c, i, j;
	
	k = (double*) malloc(g->sizeGraph*sizeof(double));
	for(t=0; t<g->sizeTable; t++) {
		double pi;
		m = 0.;
		for(i=0; i<g->sizeGraph; i++) {
			k[i] = 0.;
			for(j=0; j<g->sizeGraph; j++)
				k[i] += g->edge[t][i][j];
			m += k[i];
		}
		m /= 2.;
		pi = ((double)m)/(((double)g->sizeGraph)*((double)g->sizeGraph-1));
		pl = getPartitionList(part);
//	printf("X pi %.2lE, ke %.0lE, kf %.0lE, %.2lE\n", pi, k[0], k[2], expERExp(pi,k[0],k[2],g->sizeGraph));
		for(c=0; c<part->sizeAtom; c++) {
			int e, f;
			for(e=pl->start[c]; e>=0; e=pl->next[e]) {
				for(f=pl->next[e]; f>=0; f=pl->next[f]) {
//printf("X e %d f %d %d %.2lE\n", e, f, g->edge[t][e][f], ((double) g->edge[t][e][f])-expERExp(pi,k[e],k[f],g->sizeGraph));
						res += ((double) g->edge[t][e][f])-expERExp(pi,k[e],k[f],g->sizeGraph);
				}
			}
		}
	}
	free((void*)k);
	freePartitionList(pl);
	return res;
}

double getERFullModularityMulti(TypeMultiGraph *g) {
	double res = 0., *k, m, pi;
	int t, c, i, j;
	
	k = (double*) malloc(g->sizeGraph*sizeof(double));
	for(t=0; t<g->sizeTable; t++) {
		double pi;
		m = 0.;
		for(i=0; i<g->sizeGraph; i++) {
			k[i] = 0.;
			for(j=0; j<g->sizeGraph; j++)
				k[i] += g->edge[t][i][j];
			m += k[i];
		}
		m /= 2.;
		pi = ((double)m)/(((double)g->sizeGraph)*((double)g->sizeGraph-1));
		for(i=0; i<g->sizeGraph; i++)
			for(j=0; j<i; j++)
				res += ((double) g->edge[t][i][j])-expERExp(pi,k[i],k[j],g->sizeGraph);
	}
	free((void*)k);
	return res;
}

double getERFullModularityMultiOne(int t, TypeMultiGraph *g) {
	double res = 0., *k, m, pi;
	int c, i, j;
	
	k = (double*) malloc(g->sizeGraph*sizeof(double));
	m = 0.;
	for(i=0; i<g->sizeGraph; i++) {
		k[i] = 0.;
		for(j=0; j<g->sizeGraph; j++)
			k[i] += g->edge[t][i][j];
		m += k[i];
	}
	m /= 2.;
	pi = ((double)m)/(((double)g->sizeGraph)*((double)g->sizeGraph-1));
	for(i=0; i<g->sizeGraph; i++)
		for(j=0; j<i; j++)
			res += ((double) g->edge[t][i][j])-expERExp(pi,k[i],k[j],g->sizeGraph);
	free((void*)k);
	return res;
}

TypeGraph *readMatrixGraph(FILE *f) {
	char c, tmpS[MAX_NAME_SIZE+1], tmpE[MAX_NAME_SIZE+1];
	int sizeTmp = 0, sizeBuffer, n, m;
	TypeEdgeG *buffer;
	TypeGraph *g;

	g = (TypeGraph*) malloc(sizeof(TypeGraph));
	sizeBuffer = BASIC_INC_BUFFER;
	buffer = (TypeEdgeG *) monmalloc(sizeBuffer*sizeof(TypeEdgeG));
	for(c=getc(f); c!=EOF && issepline(c); c = fgetc(f));
	if(c != EOF) { //read first line
		int i, indexS, indexE;
		if(c == '\'' || c == '"') {
			c = fgetc(f);
			for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
				tmpS[i] = c;
				c = fgetc(f);
			}
			if(c == '\'' || c == '"')
				c = fgetc(f);
			else
				exitProg(ErrorReading, "Missing closing \" or '...");			
		} else {
			for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issepline(c); i++) {
				tmpS[i] = c;
				c = fgetc(f);
			}
		}
		if(i == 0)
			exitProg(ErrorExec, "Empty name ...");
		if(i == MAX_NAME_SIZE)
			exitProg(ErrorExec, "Name too much long...");
		tmpS[i++] = '\0';
//printf("[%s]\n", tmpS);
		for(; c != EOF && issep(c); c = fgetc(f));
		while(c != EOF && !isline(c)) {
			for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issepline(c); i++) {
				tmpE[i] = c;
				c = fgetc(f);
			}
			if(i == 0)
				exitProg(ErrorExec, "Empty value ...");
			if(i == MAX_NAME_SIZE)
				exitProg(ErrorExec, "Value too much long...");
			tmpE[i++] = '\0';
//printf("[%s]\n", tmpE);
			if(sizeTmp >= sizeBuffer) {
				sizeBuffer += BASIC_INC_BUFFER;
				buffer = (TypeEdgeG *) monrealloc((void *) buffer, sizeBuffer*sizeof(TypeEdgeG));
			}
			buffer[sizeTmp++] = (TypeEdgeG) atof(tmpE);
			for(; c != EOF && issep(c); c = fgetc(f));
		}
	}
	g->sizeGraph = sizeTmp;
	g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
	g->edge = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
	g->name[0] = (char*) malloc((strlen(tmpS)+1)*sizeof(char));
	strcpy(g->name[0], tmpS);
	g->edge[0] = (TypeEdgeG*) realloc(buffer, sizeTmp*sizeof(TypeEdgeG));
	for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
	for(n=1; c != EOF && n<g->sizeGraph; n++) {
		int i;
		if(c == '\'' || c == '"') {
			c = fgetc(f);
			for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
				tmpS[i] = c;
				c = fgetc(f);
			}
			if(c == '\'' || c == '"')
				c = fgetc(f);
			else
				exitProg(ErrorReading, "Missing closing \" or '...");			
		} else {
			for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issepline(c); i++) {
				tmpS[i] = c;
				c = fgetc(f);
			}
		}
		if(i == 0)
			exitProg(ErrorExec, "Empty name ...");
		if(i == MAX_NAME_SIZE)
			exitProg(ErrorExec, "Name too much long...");
		tmpS[i++] = '\0';
//printf("[%s]\n", tmpS);
		g->name[n] = (char*) malloc((strlen(tmpS)+1)*sizeof(char));
		strcpy(g->name[n], tmpS);
		g->edge[n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
		for(; c != EOF && issep(c); c = fgetc(f));
		for(m=0; c != EOF && !isline(c) && m<g->sizeGraph; m++) {
			double valTmp;
			for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issepline(c); i++) {
				tmpE[i] = c;
				c = fgetc(f);
			}
			if(i == 0)
				exitProg(ErrorExec, "Empty value ...");
			if(i == MAX_NAME_SIZE)
				exitProg(ErrorExec, "Value too much long...");
			tmpE[i++] = '\0';
//printf("[%s]\n", tmpE);
			g->edge[n][m] = (TypeEdgeG) atof(tmpE);
			for(; c != EOF && issep(c); c = fgetc(f));
		}
		if(m<g->sizeGraph)
			exitProg(ErrorExec, "Missing value...");
		for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
	}
	if(n<g->sizeGraph) {
		printf("%d\t%d\n", m, g->sizeGraph);
		exitProg(ErrorExec, "Missing line...");
	}
	return g;
}


/*read graph in standard format*/
TypeGraph *readGraph(FILE *f) {
	char c;
	TypeGraph *g;
	TypeSizeG sizebuf = INC_SIZE_BUF, sizeGraph, *start, *end, n;
	TypeLexiTree *dict;
	dict = newLexiTree();
	g = (TypeGraph*) malloc(sizeof(TypeGraph));
	start = (TypeSizeG*) malloc(sizebuf*sizeof(TypeSizeG));
	end = (TypeSizeG*) malloc(sizebuf*sizeof(TypeSizeG));
	sizeGraph = 0;
	for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
	while(c != EOF) {
		char tmpS[MAX_NAME_SIZE+1], tmpE[MAX_NAME_SIZE+1];
		int i, indexS, indexE;
		if(c == '\'' || c == '"') {
			c = fgetc(f);
			for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
				tmpS[i] = c;
				c = fgetc(f);
			}
			if(c == '\'' || c == '"')
				c = fgetc(f);
			else
				exitProg(ErrorReading, "Missing closing \" or '...");			
		} else {
			for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issepline(c); i++) {
				tmpS[i] = c;
				c = fgetc(f);
			}
		}
/*		if(i == 0)
			exitProg(ErrorExec, "Empty name...");
*/		if(i == MAX_NAME_SIZE)
			exitProg(ErrorExec, "Name too much long...");
		if(i>0) {
			tmpS[i++] = '\0';
//printf("%s-%d\t->\t", tmp, indexS);
			for(; c != EOF && issep(c); c = fgetc(f));
			if(c == '\'' || c == '"') {
				c = fgetc(f);
				for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
					tmpE[i] = c;
					c = fgetc(f);
				}
				if(c == '\'' || c == '"')
					c = fgetc(f);
				else
					exitProg(ErrorReading, "Missing closing \" or '...");			
			} else {
				for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issepline(c); i++) {
					tmpE[i] = c;
					c = fgetc(f);
				}
			}
/*			if(i == 0)
				exitProg(ErrorExec, "Empty name...");
*/			if(i == MAX_NAME_SIZE)
				exitProg(ErrorExec, "Name too much long...");
			if(i>0) {
				tmpE[i++] = '\0';
				indexS = addWordLexiTree(tmpS, dict);
				indexE = addWordLexiTree(tmpE, dict);
//printf("%s-%d\n", tmp, indexE);
				if(sizeGraph >= sizebuf) {
					sizebuf += INC_SIZE_BUF;
					start = (TypeSizeG*) realloc((void*)start, sizebuf*sizeof(TypeSizeG));
					end = (TypeSizeG*) realloc((void*)end, sizebuf*sizeof(TypeSizeG));
				}
				start[sizeGraph] = indexS;
				end[sizeGraph] = indexE;
				sizeGraph++;
			}
		}
		for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
	}
	g->sizeGraph = dict->sizeWord;
	g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
	fillLexiTree(g->name, dict);
	g->edge = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
	for(n=0; n<g->sizeGraph; n++) {
		TypeSizeG m;
		g->edge[n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
		for(m=0; m<g->sizeGraph; m++)
			g->edge[n][m] = 0;
	}
	for(n=0; n<sizeGraph; n++) {
		g->edge[start[n]][end[n]] = 1;
		g->edge[end[n]][start[n]] = 1;
	}
	free((void*)start);
	free((void*)end);
	freeLexiTree(dict);
	return g;
}

/*cpy multi graph in standard format*/	
TypeMultiGraph *cpyMultiGraph(TypeMultiGraph *h) {
	TypeMultiGraph *g;
	TypeSizeG i, n, m;;
	int t;
	g = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
	g->sizeGraph = h->sizeGraph;
	g->sizeTable = h->sizeTable;
	g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
	for(n=0; n<g->sizeGraph; n++) {
		g->name[n] = (char*) malloc((strlen(h->name[n])+1)*sizeof(char));
		strcpy(g->name[n], h->name[n]);
	}
	g->edge = (TypeEdgeG***) malloc(g->sizeTable*sizeof(TypeEdgeG**));
	g->present = (int**) malloc(g->sizeTable*sizeof(int*));
	for(t=0; t<g->sizeTable; t++) {
		g->edge[t] = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
		g->present[t] = (int*) malloc(g->sizeGraph*sizeof(int));
		for(n=0; n<g->sizeGraph; n++) {
			g->present[t][n] = h->present[t][n];
			g->edge[t][n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
			for(m=0; m<g->sizeGraph; m++)
				g->edge[t][n][m] = h->edge[t][n][m];
		}
	}
	return g;
}

/*free multi graph in standard format*/	
void freeMultiGraph(TypeMultiGraph *g) {
	TypeSizeG i;
	int t;
	for(i=0; i<g->sizeGraph; i++)
		free((void*)g->name[i]);
	free((void*)g->name);
	for(t=0; t<g->sizeTable; t++) {
		free((void*)g->present[t]);
		for(i=0; i<g->sizeGraph; i++)
			free((void*)g->edge[t][i]);
		free((void*)g->edge[t]);
	}
	free((void*)g->present);
	free((void*)g->edge);
	free((void*) g);
}

/*free  graph in standard format*/	
void freeGraph(TypeGraph *g) {
	TypeSizeG i;
	for(i=0; i<g->sizeGraph; i++)
		free((void*)g->name[i]);
	free((void*)g->name);
	for(i=0; i<g->sizeGraph; i++)
		free((void*)g->edge[i]);
	free((void*)g->edge);
}
/*print graph in standard format*/	
void fprintGraph(FILE *f, TypeGraph *g) {
	TypeSizeG n, m;
	for(n=0; n<g->sizeGraph; n++)
		for(m=0; m<n; m++) {
			if(g->edge[n][m])
				fprintf(f, "%s\t%s\n", g->name[n], g->name[m]);
		}
}

/*print graph in standard format*/	
void fprintGraphTable(FILE *f, TypeGraph *g) {
	TypeSizeG n, m;
	for(n=0; n<g->sizeGraph; n++) {
		fprintf(f, "%d", (int) g->edge[n][0]);
		for(m=1; m<g->sizeGraph; m++)
		 fprintf(f, " %d", (int) g->edge[n][m]);
		fprintf(f, "\n");
	}
}

/*print graph in standard format*/	
void fprintGraphInci(FILE *f, TypeGraph *g) {
	TypeSizeG n, m;
	for(n=0; n<g->sizeGraph; n++) {
		fprintf(f, "%d", (int) g->edge[n][0]);
		for(m=1; m<g->sizeGraph; m++)
			fprintf(f, "\t%d", (int) g->edge[n][m]);
		fprintf(f, "\n");
	}
}

/*print ith graph in standard format*/	
void fprintMultiGraph(FILE *f, int i, TypeMultiGraph *g) {
	TypeSizeG n, m;
	for(n=0; n<g->sizeGraph; n++) {
		for(m=0; m<g->sizeGraph; m++) {
			if(g->edge[i][n][m])
				fprintf(f, "%s\t%s\n", g->name[n], g->name[m]);
		}
	}
}

/*print ith graph in standard format*/	
void fprintMultiGraphTable(FILE *f, int t, TypeMultiGraph *g) {
	TypeSizeG n, m;
	for(n=0; n<g->sizeGraph; n++) {
			fprintf(f, "%d", (int) g->edge[t][n][0]);
			for(m=1; m<g->sizeGraph; m++)
		 	fprintf(f, " %d", (int) g->edge[t][n][m]);
			fprintf(f, "\n");
		}
}
/*print ith graph in standard format*/	
void fprintMultiGraphDebug(FILE *f, TypeMultiGraph *g) {
	TypeSizeG n, m;
	int t;
	for(t=0; t<g->sizeTable; t++) {
		for(n=0; n<g->sizeGraph; n++) {
			fprintf(f, "%d", (int) g->edge[t][n][0]);
			for(m=1; m<g->sizeGraph; m++)
		 	fprintf(f, " %d", (int) g->edge[t][n][m]);
			fprintf(f, "\n");
		}
		fprintf(f, "\n");
	}
}

int countEdge(TypeMultiGraph *g) {
		TypeSizeG n, m;
		int t, tot = 0;
		for(n=0; n<g->sizeGraph; n++)
			for(m=0; m<n; m++) {
				for(t=0; t<g->sizeTable && g->edge[t][n][m] == 0; t++);
				if(t < g->sizeTable)
					tot++;
			}
	return tot;
}


TypeMultiGraph *toMultiGraph(TypeGraph **graph, int sizeTable) {
	int i, t;
	TypeMultiGraph *multi;
	TypeLexiTree *dict;
	
	multi = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
	dict = newLexiTree();
	for(t=0; t<sizeTable; t++) {
		TypeSizeG j;
		for(j=0; j<graph[t]->sizeGraph; j++)
			addWordLexiTree(graph[t]->name[j], dict);
	}
	multi->sizeTable = sizeTable;
	multi->sizeGraph = dict->sizeWord;
	multi->name = (char**) malloc(multi->sizeGraph*sizeof(char*));
	fillLexiTree(multi->name, dict);
	multi->edge = (TypeEdgeG***) malloc(multi->sizeTable*sizeof(TypeEdgeG**));
	multi->present = (int**) malloc(multi->sizeTable*sizeof(int*));
	for(t=0; t<multi->sizeTable; t++) {
		TypeSizeG j;
		multi->edge[t] = (TypeEdgeG**) malloc(multi->sizeGraph*sizeof(TypeEdgeG*));
		multi->present[t] = (int*) malloc(multi->sizeGraph*sizeof(int));
		for(i=0; i<multi->sizeGraph; i++) {
			TypeSizeG k;
			multi->present[t][i] = 0;
			multi->edge[t][i] = (TypeEdgeG*) malloc(multi->sizeGraph*sizeof(TypeEdgeG));
			for(j=0; j<multi->sizeGraph; j++)
				multi->edge[t][i][j] = 0;
		}
		for(j=0; j<graph[t]->sizeGraph; j++) {
			TypeSizeG k, indexJ = findWordLexi(graph[t]->name[j], dict);
			if(indexJ>=0) {
				multi->present[t][indexJ] = 1;
				for(k=0; k<j; k++) {
					TypeSizeG indexK = findWordLexi(graph[t]->name[k], dict);
					if(indexK>=0) {
						multi->edge[t][indexJ][indexK] = graph[t]->edge[j][k];
						multi->edge[t][indexK][indexJ] = multi->edge[t][indexJ][indexK];
					}
				}
			}
		}
	}
	freeLexiTree(dict);
	return multi;
}

TypeGraph *sumMultiGraph(TypeMultiGraph *multi) {
	TypeGraph *g;
	TypeSizeG n, m;
	char tmp[100];
	
	g = (TypeGraph*) malloc(sizeof(TypeGraph));
	g->sizeGraph = multi->sizeGraph;
	g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
	g->edge = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
	for(n=0; n<g->sizeGraph; n++) {
		int t;
		g->name[n] = (char*) malloc((strlen(multi->name[n])+1)*sizeof(char));
		strcpy(g->name[n], multi->name[n]);
		g->edge[n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
		for(m=0; m<g->sizeGraph; m++) {
			g->edge[n][m] = 0;
			for(t=0; t<multi->sizeTable; t++)
				g->edge[n][m] += multi->edge[t][n][m];
		}
	}
	return g;
}


TypeGraph *interMultiGraph(TypeMultiGraph *multi) {
	TypeGraph *g;
	TypeSizeG n, m;
	char tmp[100];
	
	g = (TypeGraph*) malloc(sizeof(TypeGraph));
	g->sizeGraph = multi->sizeGraph;
	g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
	g->edge = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
	for(n=0; n<g->sizeGraph; n++) {
		int t;
		g->name[n] = (char*) malloc((strlen(multi->name[n])+1)*sizeof(char));
		strcpy(g->name[n], multi->name[n]);
		g->edge[n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
		for(m=0; m<g->sizeGraph; m++) {
			for(t=0; t<multi->sizeTable && multi->edge[t][n][m]>0; t++);
			if(t == multi->sizeTable)
				g->edge[n][m] = 1;
			else
				g->edge[n][m] = 0;
		}
	}
	return g;
}

TypeGraph *unionMultiGraph(TypeMultiGraph *multi) {
	TypeGraph *g;
	TypeSizeG n, m;
	char tmp[100];
	
	g = (TypeGraph*) malloc(sizeof(TypeGraph));
	g->sizeGraph = multi->sizeGraph;
	g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
	g->edge = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
	for(n=0; n<g->sizeGraph; n++) {
		int t;
		g->name[n] = (char*) malloc((strlen(multi->name[n])+1)*sizeof(char));
		strcpy(g->name[n], multi->name[n]);
		g->edge[n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
		for(m=0; m<g->sizeGraph; m++) {
			for(t=0; t<multi->sizeTable && multi->edge[t][n][m]==0; t++);
			if(t == multi->sizeTable)
				g->edge[n][m] = 0;
			else
				g->edge[n][m] = 1;
		}
	}
	return g;
}

/*print ith graph in standard format*/	
void fixEdgeMultiGraph(TypeMultiGraph *g) {
	TypeSizeG n, m;
	int t;
	for(t=0; t<g->sizeTable; t++)
		for(n=0; n<g->sizeGraph; n++) {
			for(m=0; m<g->sizeGraph; m++)
				g->edge[t][n][m] = (g->edge[t][n][m]!=0)?1:0;
			g->edge[t][n][n] = 0;
		}
}

/*print ith graph in standard format*/	
void fixEdgeGraph(TypeGraph *g) {
	TypeSizeG n, m;
		for(n=0; n<g->sizeGraph; n++) {
			for(m=0; m<g->sizeGraph; m++)
				g->edge[n][m] = (g->edge[n][m]!=0)?1:0;
			g->edge[n][n] = 0;
		}
}



TypeMultiGraph *sumMultiGraphMulti(TypeMultiGraph *multi) {
	TypeMultiGraph *g;
	TypeSizeG n, m;
	
	g = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
	g->sizeTable = 1;
	g->sizeGraph = multi->sizeGraph;
	if(multi->name != NULL) {
		g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
		for(n=0; n<multi->sizeGraph; n++)
			if(multi->name[n] != NULL) {
				g->name[n] = (char*) malloc((strlen(multi->name[n])+1)*sizeof(char));
				strcpy(g->name[n], multi->name[n]);
			} else
				g->name[n] = NULL;
	} else
		g->name = NULL;
	g->edge = (TypeEdgeG***) malloc(sizeof(TypeEdgeG**));
	g->present = (int**) malloc(sizeof(int*));
	g->edge[0] = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
	g->present[0] = (int*) malloc(g->sizeGraph*sizeof(int));
	for(n=0; n<g->sizeGraph; n++) {
		int t;
		for(t=0; t<multi->sizeTable && multi->present[t][n] == 0; t++);
		g->present[0][n] = t<multi->sizeTable;
		g->edge[0][n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
		for(m=0; m<g->sizeGraph; m++) {
			g->edge[0][n][m] = 0;
			for(t=0; t<multi->sizeTable; t++)
				g->edge[0][n][m] += multi->edge[t][n][m];
		}
	}
	return g;
}

TypeMultiGraph *interMultiGraphMulti(TypeMultiGraph *multi) {
	TypeMultiGraph *g;
	TypeSizeG n, m;
	
	g = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
	g->sizeTable = 1;
	g->sizeGraph = multi->sizeGraph;
	if(multi->name != NULL) {
		g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
		for(n=0; n<multi->sizeGraph; n++)
			if(multi->name[n] != NULL) {
				g->name[n] = (char*) malloc((strlen(multi->name[n])+1)*sizeof(char));
				strcpy(g->name[n], multi->name[n]);
			} else
				g->name[n] = NULL;
	} else
		g->name = NULL;
	g->edge = (TypeEdgeG***) malloc(sizeof(TypeEdgeG**));
	g->present = (int**) malloc(sizeof(int*));
	g->edge[0] = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
	g->present[0] = (int*) malloc(g->sizeGraph*sizeof(int));
	for(n=0; n<g->sizeGraph; n++) {
		int t;
		for(t=0; t<multi->sizeTable && multi->present[t][n] == 0; t++);
		g->present[0][n] = t<multi->sizeTable;
		g->edge[0][n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
		g->edge[0][n][n] = 0;
		for(m=0; m<n; m++) {
			if(g->present[0][m]) {
				for(t=0; t<multi->sizeTable && multi->edge[t][n][m]>0; t++);
				if(t == multi->sizeTable)
					g->edge[0][n][m] = 1;
				else
					g->edge[0][n][m] = 0;
				g->edge[0][m][n] = g->edge[0][n][m];
			} else {
				g->edge[0][n][m] = 0;
				g->edge[0][m][n] = 0;
			}
		}
	}
	return g;
}

TypeMultiGraph *unionMultiGraphMulti(TypeMultiGraph *multi) {
	TypeMultiGraph *g;
	TypeSizeG n, m;
	
	g = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
	g->sizeTable = 1;
	g->sizeGraph = multi->sizeGraph;
	if(multi->name != NULL) {
		g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
		for(n=0; n<multi->sizeGraph; n++)
			if(multi->name[n] != NULL) {
				g->name[n] = (char*) malloc((strlen(multi->name[n])+1)*sizeof(char));
				strcpy(g->name[n], multi->name[n]);
			} else
				g->name[n] = NULL;
	} else
		g->name = NULL;
	g->edge = (TypeEdgeG***) malloc(sizeof(TypeEdgeG**));
	g->present = (int**) malloc(sizeof(int*));
	g->edge[0] = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
	g->present[0] = (int*) malloc(g->sizeGraph*sizeof(int));
	for(n=0; n<g->sizeGraph; n++) {
		int t;
		for(t=0; t<multi->sizeTable && multi->present[t][n] == 0; t++);
		g->present[0][n] = t<multi->sizeTable;
		g->edge[0][n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
		g->edge[0][n][n] = 0;
		for(m=0; m<n; m++) {
			if(g->present[0][m]) {
				for(t=0; t<multi->sizeTable && multi->edge[t][n][m]==0; t++);
				if(t == multi->sizeTable)
					g->edge[0][n][m] = 0;
				else
					g->edge[0][n][m] = 1;
				g->edge[0][m][n] = g->edge[0][n][m];
			} else {
				g->edge[0][n][m] = 0;
				g->edge[0][m][n] = 0;
			}
		}
	}
	return g;
}
