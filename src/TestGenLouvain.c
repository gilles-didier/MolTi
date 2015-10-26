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
//#include "PartitionConsensus.h"
#include "Louvain.h"

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

//./testgl -i 100 -t 10 -c 5 -p 0.001 0.01 -p 0.5 0.2 -g 100

//./testgl -i 100 -t 5 -c 5  -p 0.4 0.1 -p 0.01 0.001 -g 200 -l "../" CC
//./test -i 1000 -t 10 -c 20 -p 0.2 0.1 -p 0.9 0.8 -p 0.6 0.4 -g 1000
//./tx -i 100 -t 10 -c 20 -p 0.2 0.1 -p 0.9 0.8 -p 0.5 0.4 -g 400
//./tx -i 100 -t 10 -c 2 -p 0.2 0.1 -p 0.9 0.8 -p 0.5 0.4 -g 10
//./tx -i 1 -t 10 -c 2 -p 0.2 0.1 -p 0.9 0.8 -p 0.5 0.4 -g 8

static int trial=TRIAL, maxIter=MAX_ITER;
static char *genLouvainPath;

#define HELPMESSAGE "\ntestgl first simulates random multiplex networks with g vertices and from 1 to t layers with a balanced community structure of c communities (i simulations for each number of layers). It next detects communities by using aggregation and multiplex-modularity approaches and computes the adjusted Rand index (and the normalized mutual information) between the reference community structure and the detected one. It writes a '.csv' file containing the means and the standard deviations of the adjusted Rand indexes for each method, with the following format:\n Column 1: number of layers\n Column 2 and 3: mean and standard deviation obtained with the multiplex-modularity approach\n Column 4 and 5: mean and standard deviation obtained with the sum-aggregation approach\n Column 6 and 7: mean and standard deviation obtained with the intersection\n Column 6 and 7: mean and standard deviation obtained with the union aggregation approach. \n\nusage: test <options> <output file name>\n\nOptions:\n	-g <number>	set the number of vertices of the random graphs (50 as default)\n	-t <number>	set the max number of random graphs (10 as default)\n	-c <number>	set the number of classes (3 as default)\n	-i <number>	set the number of iterations (4 as default)\n	-p <prob intra> <prob inter>	add a new pair of probas\n	-a  <number>	set the modularity parameter (1 as default)\n	-h	display help message\n"


static TypePartition getPartitionGenLouvain(TypeMultiGraph *graph, void *info, char *genlouvainpath, char *id);


int main(int argc, char **argv) {		
	char inputFileName[STRING_SIZE], outputFileName[STRING_SIZE], outputFileNameTree[STRING_SIZE], outputFileNameFossil[STRING_SIZE], inputFileNameModel[STRING_SIZE], option[256], kind = 'S';
	FILE *fi, *fo, *fs;
	int i, j, width=1000, height=500, Y_MAX=height-50, absolute = 1, modelF=0, outF=0, nTable, maxTable = 10, size = 50, niter = NITER, sizeAtom = 3, nsample = 50, oneModel = 0, reference = 0, np = 0;
	TypeMultiERMG *model, *model2, **modelA;
	TypeMultiGraph *graph;
	TypePartition part1, part2, part;
	double **tau, ***tauA, like, *alpha, scale = 1., pin[NP_MAX], pex[NP_MAX], alp = 1., *pi, *pe, p, thre = 0., pres = 1.;
	TypePartitionMethod type = LouvainType;
	
	genLouvainPath = "./";
	
//	srand(time(NULL));
	for(i=0; i<256; i++)
		option[i] = 0;
	   
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
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &thre) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a real number is required after option -f");
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
		if(option['a']) {
			option['a'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &alp) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a real number is required after option -a");
		}
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &pres) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a real number is required after option -x");
		}
		if(option['l']) {
			option['l'] = 0;
			genLouvainPath = argv[++i];
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

	if(!(i<argc && sscanf(argv[i++], "%s", outputFileName) == 1)) {
		sprintf(outputFileName, "M_%c_%s_%d_%d_%d_%d.%s", kind, NAME_OUTPUT, sizeAtom, size, maxTable, niter, EXT_OUTPUT);
	}
	alpha = (double*) malloc(sizeAtom*sizeof(double));
	for(i=0;i<sizeAtom; i++)
		alpha[i] = 1./((double)sizeAtom);


	if(fo = fopen(outputFileName, "w")) {
		double *randInd, mean, std;
		double *randIndS, meanS, stdS;
		double *randIndI, meanI, stdI;
		double *randIndU, meanU, stdU;
		double *randIndG, meanG, stdG;
		randInd = (double*) malloc(niter*sizeof(double));
		randIndS = (double*) malloc(niter*sizeof(double));
		randIndI = (double*) malloc(niter*sizeof(double));
		randIndU = (double*) malloc(niter*sizeof(double));
		randIndG = (double*) malloc(niter*sizeof(double));
		pe = (double*) malloc((maxTable+1)*sizeof(double));
		pi = (double*) malloc((maxTable+1)*sizeof(double));
		setBalancedPartition(sizeAtom, size, &part1);
		p = getPartitionDensity(&part1);
		for(nTable=1; nTable<maxTable; nTable++) {
			mean = 0.;
			meanS = 0.;
			meanI = 0.;
			meanU = 0.;
			meanG = 0.;
			printf("simulating %d graph\n", nTable);
			i = 0;
			while(i<niter) {
				double dist, x;
				TypeMultiGraph *gtmp;
				int t;
				char ident[STRING_SIZE];
				sprintf(ident, "A_%d_%d", nTable, i);

				model = getBasicMultiERMGX(sizeAtom, nTable, alpha, pin, pex, np);
				for(t=0; t<nTable; t++) {
					pe[t] = model->pi[t][0][1];
					pi[t] = model->pi[t][0][0];
				}
//				graph = getRandomMultiGraph(model, size, &part1);
				graph = getRandomMultiGraphMissing(model, size, &part1, pres);
		 		freeMultiERMG(model);
		 		
				part2 = getPartitionGenLouvain(graph, NULL, genLouvainPath, ident);
				if(part2.sizeItem>0) {
					randIndG[i] = comparePart(correctedRandIndex, &part1, &part2);
					free((void*) part2.atom);
					meanG += randIndG[i];
					
					part2 = getPartition(graph, LouvainType, &alp);
					randInd[i] = comparePart(correctedRandIndex, &part1, &part2);
					free((void*) part2.atom);
					mean += randInd[i];
					
					gtmp = sumMultiGraphMulti(graph);
					part2 = getPartition(gtmp, LouvainType, &alp);
					freeMultiGraph(gtmp);
					randIndS[i] = comparePart(correctedRandIndex, &part1, &part2);
					free((void*) part2.atom);
					meanS += randIndS[i];

					gtmp = interMultiGraphMulti(graph);
					part2 = getPartition(gtmp, LouvainType, &alp);
					freeMultiGraph(gtmp);
					randIndI[i] = comparePart(correctedRandIndex, &part1, &part2);;
					free((void*) part2.atom);
					meanI += randIndI[i];
					
					gtmp = unionMultiGraphMulti(graph);
					part2 = getPartition(gtmp, LouvainType, &alp);
					freeMultiGraph(gtmp);
					randIndU[i] = comparePart(correctedRandIndex, &part1, &part2);
					free((void*) part2.atom);
					meanU += randIndU[i];
					i++;
				}
				freeMultiGraph(graph);
			}
			mean /= (double) niter;
			std = 0.;
			for(i=0;i<niter; i++)
				std += pow(randInd[i]-mean, 2.);
			std /= (double) niter -1.;
			std = sqrt(std);
			
			meanS /= (double) niter;
			stdS = 0.;
			for(i=0;i<niter; i++)
				stdS += pow(randIndS[i]-mean, 2.);
			stdS /= (double) niter -1.;
			stdS = sqrt(stdS);
			
			meanI /= (double) niter;
			stdI = 0.;
			for(i=0;i<niter; i++)
				stdI += pow(randIndI[i]-mean, 2.);
			stdI /= (double) niter -1.;
			stdI = sqrt(stdI);
			
			meanU /= (double) niter;
			stdU = 0.;
			for(i=0;i<niter; i++)
				stdU += pow(randIndU[i]-mean, 2.);
			stdU /= (double) niter -1.;
			stdU = sqrt(stdU);
			
			meanG /= (double) niter;
			stdG = 0.;
			for(i=0;i<niter; i++)
				stdG += pow(randIndU[i]-mean, 2.);
			stdG /= (double) niter -1.;
			stdG = sqrt(stdG);

			fprintf(fo, "%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE", nTable, mean, std, meanS, stdS, meanI, stdI, meanU, stdU, meanG, stdG);
			fprintf(fo, "\n");
			fprintf(stdout, "%d\t%lE\t%lE\t%lE\t%lE\t%lE\\n", nTable, mean, meanS, meanI, meanU, meanG);

		}
		free((void*) part1.atom);
		free((void*)randInd);
		free((void*)randIndS);
		free((void*)randIndI);
		free((void*)randIndU);
		fclose(fo);
	} else
		exitProg(ErrorWriting, "Can't read the file");
	exit(0);
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
			part.sizeItem = 0;
			part.atom = NULL;
			return part;
		}
	}
	sprintf(fileName, "%s_octave.m", id);
	if(f = fopen(fileName, "w")) {
		fprintf(f, "addpath(\"%s\");\nload genlouvain.m;\n", genlouvainpath);
		if(graph->sizeTable <= 1) {
			fprintf(f, "load %s_0.mat;\n", id);
			fprintf(f, "N=length(%s_0);gamma = 1;\nk = full(sum(%s_0));\ntwom = sum(k);\nB = full(%s_0 - gamma*k'*k/twom);\n[S,Q] = genlouvain(B);\nQ = Q/twom\n", id, id, id);
		} else {
			for(t=0; t<graph->sizeTable; t++)
				fprintf(f, "load %s_%d.mat;\n", id, t);
			fprintf(f, "A = cell(%d);\nA={%s_%d", graph->sizeTable, id, 0);
			for(t=1; t<graph->sizeTable; t++)
				fprintf(f, ",%s_%d", id, t);
//			fprintf(f, "};\nN=length(A{1});\nT=length(A);\ngamma = 1;\nomega = T+1\nB=spalloc(N*T,N*T,N*N*T+2*N*T);\ntwomu=0;\nfor s=1:T\nk=sum(A{s});\ntwom=sum(k);\ntwomu=twomu+twom;\nindx=[1:N]+(s-1)*N;\nB(indx,indx)=A{s}-gamma*k'*k/twom;\nend\ntwomu=twomu+2*omega*N*(T-1);\nB = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);\n[S,Q] = genlouvain(B, 10000, 0, 0);\nQ = Q/twomu\nS = reshape(S,N,T);\n");
			fprintf(f, "};gamma = 1;\nomega = 1;\nN=length(A{1});\nT=length(A);\nB=spalloc(N*T,N*T,(N+T)*N*T);\ntwomu=0;\nfor s=1:T\n k=sum(A{s});\n twom=sum(k);\n twomu=twomu+twom;\n    indx=[1:N]+(s-1)*N;\n B(indx,indx)=A{s}-gamma*k'*k/twom;\nend\ntwomu=twomu+T*omega*N*(T-1);\nall2all = N*[(-T+1):-1,1:(T-1)];\nB = B + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);\n[S,Q] = genlouvain(B, 10000, 0, 0);\nQ = Q/twomu\nS = reshape(S,N,T);\n");
		}
		fprintf(f, "filename = \"%s_result.parti\";\nfid = fopen (filename, \"w\");\nfor s=1:N\nfprintf(fid,\"%%d\\t%%d\\n\", s-1, S(s,1)-1);\nend\nfclose(fid)", id);
		fclose(f);
	} else {
		fprintf(stderr, "Can't write the file %s", fileName);
		exit(1);
		part.sizeItem = 0;
		part.atom = NULL;
		return part;
	}
fprintf(stderr, "Computing %s\n", id);
//	sprintf(command, "octave %s_octave.m > /dev/null", id);
	sprintf(command, "octave %s_octave.m", id);
	if(system(command)<0) {
		fprintf(stderr, "Execution error on command '%s'\n", command);
		part.sizeItem = 0;
		part.atom = NULL;
		return part;
	}
	sprintf(fileName, "%s_result.parti", id);
	if(f = fopen(fileName, "r")) {
		part = readPartitionNumber(f);
	} else {
		fprintf(stderr, "Can't read the file %s\n", fileName);
		part.sizeItem = 0;
		part.atom = NULL;
		return part;
	}
/*clean up stuff*/
fprintf(stderr, "Clean files\n");
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
