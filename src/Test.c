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
#include <nlopt.h>
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
#define TOLERANCE_CONSTRAINT 0.001
#define TOLERANCE_OPTIM 0.001
#define NLOPT_ALGO NLOPT_LN_BOBYQA
//#define NLOPT_ALGO NLOPT_LN_COBYLA
//#define NLOPT_ALGO NLOPT_GN_ISRES

//./test -i 1000 -t 10 -c 20 -p 0.2 0.1 -p 0.9 0.8 -p 0.6 0.4 -g 1000
//./tx -i 100 -t 10 -c 20 -p 0.2 0.1 -p 0.9 0.8 -p 0.5 0.4 -g 400
//./tx -i 100 -t 10 -c 2 -p 0.2 0.1 -p 0.9 0.8 -p 0.5 0.4 -g 10
//./tx -i 1 -t 10 -c 2 -p 0.2 0.1 -p 0.9 0.8 -p 0.5 0.4 -g 8

static int trial=TRIAL, maxIter=MAX_ITER;
static double tolConstraint = TOLERANCE_CONSTRAINT, tolOptim = TOLERANCE_OPTIM;
static nlopt_algorithm algorithm = NLOPT_ALGO;

#define HELPMESSAGE "\nusage: draw <input file> [<output file>]\n\nThe input file has to be in Diagonal format, \nreturns a figure of the N-map in <input file> in eps format.\n"




int main(int argc, char **argv) {		
	char inputFileName[STRING_SIZE], outputFileName[STRING_SIZE], outputFileNameTree[STRING_SIZE], outputFileNameFossil[STRING_SIZE], inputFileNameModel[STRING_SIZE], option[256], kind = 'S';
	FILE *fi, *fo, *fs;
	int i, j, width=1000, height=500, Y_MAX=height-50, absolute = 1, modelF=0, outF=0, nTable, maxTable = 10, size = 50, niter = NITER, sizeAtom = 3, nsample = 50, oneModel = 0, reference = 0, np = 0;
	TypeMultiERMG *model, *model2, **modelA;
	TypeMultiGraph *graph;
	TypePartition part1, part2, part;
	double **tau, ***tauA, like, *alpha, scale = 1., pin[NP_MAX], pex[NP_MAX], alp = 1., *pi, *pe, p, thre = 0.;
	TypePartitionMethod type = LouvainType;
	
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
		randInd = (double*) malloc(niter*sizeof(double));
		randIndS = (double*) malloc(niter*sizeof(double));
		randIndI = (double*) malloc(niter*sizeof(double));
		randIndU = (double*) malloc(niter*sizeof(double));
		pe = (double*) malloc((maxTable+1)*sizeof(double));
		pi = (double*) malloc((maxTable+1)*sizeof(double));
		setBalancedPartition(sizeAtom, size, &part1);
		p = getPartitionDensity(&part1);
		for(nTable=1; nTable<maxTable; nTable++) {
			mean = 0.;
			meanS = 0.;
			meanI = 0.;
			meanU = 0.;
			printf("simulating %d graph\n", nTable);
			for(i=0;i<niter; i++) {
				double dist, x;
				TypeMultiGraph *gtmp;
				int t;
				model = getBasicMultiERMGX(sizeAtom, nTable, alpha, pin, pex, np);
				for(t=0; t<nTable; t++) {
					pe[t] = model->pi[t][0][1];
					pi[t] = model->pi[t][0][0];
				}
				graph = getRandomMultiGraph(model, size, &part1);
		 		freeMultiERMG(model);
		 		
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

			fprintf(fo, "%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE", nTable, mean, std, meanS, stdS, meanI, stdI, meanU, stdU);
			fprintf(fo, "\n");
			fprintf(stdout, "%d\t%lE\t%lE\t%lE\t%lE\n", nTable, mean, meanS, meanI, meanU);

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
