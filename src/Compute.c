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
#include <stdio.h>
#include <math.h>

#include "Utils.h"
#include "Graph.h"
#include "Partition.h"
#include "Louvain.h"
#include "EdgesComposition.h"

#ifdef DO_PS
#endif

#define STRING_SIZE 300
#define INC_CUT 5
#define SEQ_SIZE 30
#define NAME_OUTPUT "outFile"
#define EXT_OUTPUT "cls"
#define EXT_OUTPUT_TREE "_tree.newick"
#define EXT_OUTPUT_FOSSIL "_fossil.csv"
#define INTMIN 0.6
#define INTMAX 0.7
#define EXTMIN 0.4
#define EXTMAX 0.5
#define NTABLE 3
#define NCLASS 2
#define NITER 4
#define NTRIALS 2

#define HELPMESSAGE "\nusage: multinet [options] <input file1> <input file2> ...\n\noptions are\n\t-o <file name>\t\tset the output file name prefix\n\t-p <real number>\tset the Newman modularity parameter\n\t-s\t\t\tcompute partition on each graph individually and on the sum graph\n\t-h\t\t\tdisplay help\n"

void fillMultiOne(TypeGraph *g, TypeMultiGraph *m) {
	m->name = g->name;
	m->edge[0] = g->edge;
	m->sizeGraph = g->sizeGraph;
}


int main(int argc, char **argv) {		
	char **inputNameTable, **name, inputFileName[STRING_SIZE], outputFileName[STRING_SIZE], effectifFileName[STRING_SIZE], outputPrefix[STRING_SIZE], outputFileNameRand[STRING_SIZE], buffer[STRING_SIZE], inputFileNameModel[STRING_SIZE], *tmp, option[256];
	FILE *fi, *fo, *fs;
	int i, j, t, sizeTable, singleF = 0;
	TypeMultiGraph *graph, *gtmp;
	TypeGraph **tableGraph, *g;
	TypePartition  part, *tablePart;
	double alpha = 1.;

	for(i=0; i<256; i++)
		option[i] = 0;
	   
	sprintf(outputFileName, "%s.%s", NAME_OUTPUT, EXT_OUTPUT);
	tableGraph = (TypeGraph**) malloc((argc+3)*sizeof(TypeGraph*));
	inputNameTable = (char**) malloc((argc+3)*sizeof(char*));
	sizeTable = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['o']) {
			option['o'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%s", outputFileName) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -f");
		}
		if(option['s']) {
			option['s'] = 0;
			singleF = 1;
		}
		if(option['p']) {
			option['p'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &alpha) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a real number is required after option -p");
		}
		if(option['m']) {
			option['m'] = 0;
			if(!(sscanf(argv[i+1], "%s", inputFileName) == 1))
				exitProg(ErrorArgument, "wrong file name");
			i++;
			inputNameTable[sizeTable] = (char*) malloc((strlen(inputFileName)+1)*sizeof(char));
			strcpy(inputNameTable[sizeTable], inputFileName);
			printf("Reading file %s\n", inputNameTable[sizeTable]);
			if(fi = fopen(inputNameTable[sizeTable], "r")) {
				tableGraph[sizeTable++] = readMatrixGraph(fi);
				fclose(fi);
			} else
				exitProg(ErrorReading, inputNameTable[sizeTable]);
		}
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	for(j=i; j<argc; j++) {
		if(!(sscanf(argv[j], "%s", inputFileName) == 1))
			exitProg(ErrorArgument, "wrong file name");
		inputNameTable[sizeTable] = (char*) malloc((strlen(inputFileName)+1)*sizeof(char));
		strcpy(inputNameTable[sizeTable], inputFileName);
printf("Reading file %s\n", inputNameTable[sizeTable]);
		if(fi = fopen(inputNameTable[sizeTable], "r")) {
			tableGraph[sizeTable] = readGraph(fi);
			fclose(fi);
		} else
			exitProg(ErrorReading, inputNameTable[sizeTable]);
		fixEdgeGraph(tableGraph[sizeTable]);
printf("%d nodes\n", tableGraph[sizeTable]->sizeGraph);
		sizeTable++;
	}
	if(sizeTable <= 0)
		exitProg(ErrorArgument, "at least one graph is required.");
	strcpy(outputPrefix, outputFileName);
	if((tmp = strrchr(outputPrefix, '.')) != NULL)
		tmp[0] = '\0';
	graph = toMultiGraph(tableGraph, sizeTable);
	part = getPartition(graph, LouvainType, &alpha);
	tableGraph[sizeTable] = sumMultiGraph(graph);
	inputNameTable[sizeTable] = (char*) malloc((strlen("Sum")+1)*sizeof(char));
	strcpy(inputNameTable[sizeTable], "Sum");
	sizeTable++;
	tableGraph[sizeTable] = unionMultiGraph(graph);
	inputNameTable[sizeTable] = (char*) malloc((strlen("Union")+1)*sizeof(char));
	strcpy(inputNameTable[sizeTable], "Union");
	sizeTable++;
	tableGraph[sizeTable] = interMultiGraph(graph);
	inputNameTable[sizeTable] = (char*) malloc((strlen("Intersection")+1)*sizeof(char));
	strcpy(inputNameTable[sizeTable], "Intersection");
	sizeTable++;
	name = (char**) malloc(sizeTable*sizeof(char*));
	for(t=0; t<sizeTable; t++) {
		char *tmp;
		if((tmp = strrchr(inputNameTable[t], '.')) != NULL)
			tmp[0] = '\0';
		if((tmp=strrchr(inputNameTable[t], '/')) == NULL)
			tmp = inputNameTable[t];
		else
			tmp++;
		name[t] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
		strcpy(name[t], tmp);
//		printf("name[%d]\t%s\n", t, name[t]);
	}
	tablePart = (TypePartition*) malloc(sizeTable*sizeof(TypePartition));
	if(singleF) {
		TypeMultiGraph *gtmp;
		gtmp = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
		gtmp->sizeTable = 1;
		gtmp->edge = (TypeEdgeG***) malloc(sizeof(TypeEdgeG**));
		gtmp->present = (int**) malloc(sizeof(int*));
		gtmp->present[0] = (int*) malloc(graph->sizeGraph*sizeof(int));
		for(i=0; i<graph->sizeGraph; i++)
			gtmp->present[0][i] = 1;
		for(t=0; t<sizeTable; t++) {
			char output[STRING_SIZE], *tmp;
			fillMultiOne(tableGraph[t], gtmp);
			sprintf(output, "%s_%s.%s", outputPrefix, name[t], EXT_OUTPUT);
printf("Computing %s\n", output);
			tablePart[t] = getPartition(gtmp, LouvainType, &alpha);
			if(fo = fopen(output, "w")) {
				fprintPartitionClustNSee(fo, &(tablePart[t]), gtmp->name);
				fclose(fo);
			} else
				exitProg(ErrorWriting, output);
		}
		free((void*)gtmp->present[0]);
		free((void*)gtmp->present);
		free((void*)gtmp->edge);
		free((void*)gtmp);
	}
	if(fo = fopen(outputFileName, "w")) {
		fprintPartitionClustNSee(fo, &part, graph->name);
		fclose(fo);
	} else
		exitProg(ErrorWriting, outputFileName);
	sprintf(effectifFileName, "%s_effectif.csv", outputPrefix);
	if(fo = fopen(effectifFileName, "w")) {
		int **eff, c, t, *cs;
		eff = getEdgesNumbers(&part, graph);
		cs = getClassSize(&part);
		fprintf(fo, "\tSize");
		for(t=0; t<graph->sizeTable; t++)
			fprintf(fo, "\t%s", name[t]);
		fprintf(fo, "\n");
		for(c=0; c<part.sizeAtom; c++) {
			fprintf(fo, "ClusterID:%d", c+1);
			fprintf(fo, "\t%d", cs[c]);
			for(t=0; t<graph->sizeTable; t++) {
				fprintf(fo, "\t%d", eff[c][t]);
			}
			fprintf(fo, "\n");
		}
		fclose(fo);
		for(c=0; c<part.sizeAtom; c++)
			free((void*)eff[c]);
		free((void*)eff);
		free((void*)cs);
	} else
		exitProg(ErrorWriting, effectifFileName);
		
	if(singleF) {
		for(t=sizeTable-3; t<sizeTable; t++) {
			sprintf(effectifFileName, "%s_%s_effectif.csv", outputPrefix, name[t]);
			if((fo = fopen(effectifFileName, "w"))) {
				int **eff, c, t, *cs;
				eff = getEdgesNumbers(&(tablePart[sizeTable-1]), graph);
				cs = getClassSize(&(tablePart[sizeTable-1]));
				fprintf(fo, "\tSize");
				for(t=0; t<graph->sizeTable; t++)
					fprintf(fo, "\t%s", name[t]);
				fprintf(fo, "\n");
				for(c=0; c<tablePart[sizeTable-1].sizeAtom; c++) {
					fprintf(fo, "ClusterID:%d", c+1);
					fprintf(fo, "\t%d", cs[c]);
					for(t=0; t<graph->sizeTable; t++) {
						fprintf(fo, "\t%d", eff[c][t]);
					}
					fprintf(fo, "\n");
				}
				fclose(fo);
				for(c=0; c<tablePart[sizeTable-1].sizeAtom; c++)
					free((void*)eff[c]);
				free((void*)eff);
				free((void*)cs);
			} else
				exitProg(ErrorWriting, effectifFileName);
			sprintf(effectifFileName, "%s_%s_graph.csv", outputPrefix, name[t]);
			if((fo = fopen(effectifFileName, "w"))) {
				fprintGraph(fo, tableGraph[t]);
			} else
				exitProg(ErrorWriting, effectifFileName);
		}
	}
	if(singleF) {
		int tl, tc;
		char *tmp;
		sprintf(outputFileNameRand, "%s_Rand.csv", outputPrefix);
		if(fo = fopen(outputFileNameRand, "w")) {
			char *tmp;
			fprintf(fo,"All\t%lf\n", comparePartDiff(correctedRandIndex, &(part), graph->name, &(part), graph->name));
			for(tl=0; tl<sizeTable; tl++) {
				fprintf(fo, "%s\t%lf", name[tl], comparePartDiff(correctedRandIndex, &(part), graph->name, &(tablePart[tl]), tableGraph[tl]->name));
				for(tc=0; tc<=tl; tc++)
					fprintf(fo, "\t%lf", comparePartDiff(correctedRandIndex, &(tablePart[tc]),  tableGraph[tc]->name, &(tablePart[tl]), tableGraph[tl]->name));
				fprintf(fo, "\n");
			}
			fprintf(fo, "\tAll");
			for(tl=0; tl<sizeTable; tl++) {
				fprintf(fo, "\t%s", name[tl]);
			}
			fprintf(fo, "\n");
			fclose(fo);
		} else
			exitProg(ErrorWriting, outputFileNameRand);
	}
	for(t=0; t<sizeTable; t++) {
		freeGraph(tableGraph[t]);
		free((void*) tableGraph[t]);
		free((void*) name[t]);
		free((void*) inputNameTable[t]);
	}
	free((void*) tableGraph);
	free((void*) name);
	free((void*) inputNameTable);
	exit(0);
	return 0;
}
