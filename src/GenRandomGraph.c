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
#include <time.h>

#include "Utils.h"
#include "Graph.h"
#include "ERMG.h"
#include "RandomGraph.h"
#include "Partition.h"

#ifdef DO_PS
#endif

#define STRING_SIZE 300
#define INC_CUT 5
#define SEQ_SIZE 30
#define NAME_OUTPUT "out"
#define EXT_OUTPUT "csv"
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
#define NP_MAX 20


#define INFTY 1E99
#define TRIAL 10
#define MAX_ITER 1000



//./test -i 1000 -t 10 -c 20 -p 0.2 0.1 -p 0.9 0.8 -p 0.6 0.4 -g 1000
//./tx -i 100 -t 10 -c 20 -p 0.2 0.1 -p 0.9 0.8 -p 0.5 0.4 -g 400
//./tx -i 100 -t 10 -c 2 -p 0.2 0.1 -p 0.9 0.8 -p 0.5 0.4 -g 10
//./tx -i 1 -t 10 -c 2 -p 0.2 0.1 -p 0.9 0.8 -p 0.5 0.4 -g 8

static int trial=TRIAL, maxIter=MAX_ITER;

#define HELPMESSAGE "\ntest first simulates random multiplex networks with g vertices and from 1 to t layers with a balanced community structure of c communities (i simulations for each number of layers). It next detects communities by using aggregation and multiplex-modularity approaches and computes the adjusted Rand index (and the normalized mutual information) between the reference community structure and the detected one. It writes a '.csv' file containing the means and the standard deviations of the adjusted Rand indexes for each method, with the following format:\n Column 1: number of layers\n Column 2 and 3: mean and standard deviation obtained with the multiplex-modularity approach\n Column 4 and 5: mean and standard deviation obtained with the sum-aggregation approach\n Column 6 and 7: mean and standard deviation obtained with the intersection\n Column 6 and 7: mean and standard deviation obtained with the union aggregation approach. \n\nusage: test <options> <output file name>\n\nOptions:\n	-g <number>	set the number of vertices of the random graphs (50 as default)\n	-t <number>	set the max number of random graphs (10 as default)\n	-c <number>	set the number of classes (3 as default)\n	-i <number>	set the number of iterations (4 as default)\n	-p <prob intra> <prob inter>	add a new pair of probas\n	-a  <number>	set the modularity parameter (1 as default)\n	-h	display help message\n"


static TypePartition getPartitionGenLouvain(TypeMultiGraph *graph, void *info, char *genlouvainpath, char *id);


int main(int argc, char **argv) {		
	char inputFileName[STRING_SIZE], inputFileNameModel[STRING_SIZE], option[256], kind = 'S', *name;
	FILE *fi, *fo, *fs;
	int i, j, t, nTable = 10, maxTable = 10, size = 50, niter = NITER, sizeAtom = 3, nsample = 50, oneModel = 0, reference = 0, np = 0;
	TypeMultiERMG *model, *model2, **modelA;
	TypeMultiGraph *graph;
	double *alpha, scale = 1., pin[NP_MAX], pex[NP_MAX], alp = 1., *pi, *pe, p, thre = 0., pres = 1.;
	TypePartition part1, part2;
	
//	srand(time(NULL));
	for(i=0; i<256; i++)
		option[i] = 0;
	name = "A";
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['g']) {
			option['g'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &size) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -g");
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &maxTable) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -f");
		}
		if(option['c']) {
			option['c'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &sizeAtom) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -f");
		}
		if(option['n']) {
			option['n'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%s", name) == 1)
				i++;
			else
				exitProg(ErrorArgument, "an ident is required after option -f");
		}

		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &niter) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -f");
		}
		if(option['p']) {
			option['p'] = 0;
			if(np >= NP_MAX) {
				fprintf(stderr, "No more than %d probas\n", NP_MAX-1);
				exit(1);
			}
			if((i+1)<argc && sscanf(argv[i+1], "%lf", pin+np) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -f");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", pex+np) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -f");
			np++;
		}
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	if(np == 0) {
		pin[0] = 0.2;
		pex[0] = 0.1;
		np = 1;
	}
	alpha = (double*) malloc(sizeAtom*sizeof(double));
	for(i=0;i<sizeAtom; i++)
		alpha[i] = 1./((double)sizeAtom);
	model = getBasicMultiERMGX(sizeAtom, nTable, alpha, pin, pex, np);
	setBalancedPartition(sizeAtom, size, &part1);
	graph = getRandomMultiGraphMissing(model, size, &part1, 1.);
	freeMultiERMG(model);
	part2 = getPartitionGenLouvain(graph, NULL, "../", "coucou");
	for(i=0;i<graph->sizeGraph; i++)
		printf("%d\t%d\n", i, part2.atom[i]);
	return 0;
	for(t=0; t<nTable; t++) {
		char outputFileName[STRING_SIZE], ident[STRING_SIZE];
		sprintf(ident, "%s_%d", name, t+1);
		sprintf(outputFileName, "%s.mat", ident);
		if(fo = fopen(outputFileName, "w")) {
			fprintMultiGraphOctaveMatrix(fo, t, graph, ident);
			fclose(fo);
		} else
			fprintf(stderr, "Can't write the file %s", outputFileName);
	}
	return 0;
}


TypePartition getPartitionGenLouvain(TypeMultiGraph *graph, void *info, char *genlouvainpath, char *id) {
	char fileName[STRING_SIZE], ident[STRING_SIZE], command[STRING_SIZE];
	int t;
	FILE *f;
	TypePartition part;
	
	for(t=0; t<graph->sizeTable; t++) {
		sprintf(ident, "%s_%d", id, t);
		sprintf(fileName, "%s_%d.mat", id, t);
		if(f = fopen(fileName, "w")) {
			fprintMultiGraphOctaveMatrix(f, t, graph, ident);
			fclose(f);
		} else {
			fprintf(stderr, "Can't write the file %s", fileName);
			exit(1);
		}
	}
	sprintf(fileName, "%s_octave.m", id);
	if(f = fopen(fileName, "w")) {
		fprintf(f, "addpath(\"%s\");load genlouvain.m;\n", genlouvainpath);
		for(t=0; t<graph->sizeTable; t++)
			fprintf(f, "load %s_%d.mat;\n", id, t);
		fprintf(f, "A = cell(%d);\nA={%s_%d", graph->sizeTable, id, 0);
		for(t=1; t<graph->sizeTable; t++)
			fprintf(f, ",%s_%d", id, t);
		fprintf(f, "};\nN=length(A{1});\nT=length(A);\ngamma = 1;\nomega = T+1\nB=spalloc(N*T,N*T,N*N*T+2*N*T);\ntwomu=0;\nfor s=1:T\nk=sum(A{s});\ntwom=sum(k);\ntwomu=twomu+twom;\nindx=[1:N]+(s-1)*N;\nB(indx,indx)=A{s}-gamma*k'*k/twom;\nend\ntwomu=twomu+2*omega*N*(T-1);\nB = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);\n[S,Q] = genlouvain(B);\nQ = Q/twomu\nS = reshape(S,N,T);\n");
		fprintf(f, "filename = \"%s_result.parti\";\nfid = fopen (filename, \"w\");\nfor s=1:N\nfprintf(fid,\"%%d\\t%%d\\n\", s-1, S(s,1)-1);\nend\nfclose(fid)", id);
		fclose(f);
	} else {
		fprintf(stderr, "Can't write the file %s", fileName);
		exit(1);
	}
	sprintf(command, "octave %s_octave.m", id);
	if(system(command)<0) {
		fprintf(stderr, "Execution error on command '%s'\n", command);
		exit(1);
	}
	sprintf(fileName, "%s_result.parti", id);
	if(f = fopen(fileName, "r")) {
		part = readPartitionNumber(f);
	} else {
		fprintf(stderr, "Can't read the file %s", fileName);
		exit(1);
	}
/*clean up stuff*/
	for(t=0; t<graph->sizeTable; t++) {
		sprintf(fileName, "%s_%d.mat", id, t);
		remove(fileName);
	}
	sprintf(fileName, "%s_octave.m", id);
	remove(fileName);
	sprintf(fileName, "%s_result.parti", id);
	remove(fileName);
	return part;
}
