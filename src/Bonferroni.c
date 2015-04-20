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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "Utils.h"
#include "Annotation.h"
#include "Partition.h"
#include "StatAnnotationPartition.h"

//./annot -i 100000 -t 0.0004 Net_BP_HGNC.Annot p10.cls 


#define STRING_SIZE 300

#define HELPMESSAGE "\nusage: bonf [options] <annotation file> <partition file> [<output file>]\n\noptions are\n\t-o <ontology file>\t\tload ontology descriptors\n\t-t <real number>\tset the threshold\n\t-f c or f\t\t\tindicate the partition format\n\t-h\t\t\tdisplay help\n"


int main(int argc, char **argv) {		
	char inputFileName[STRING_SIZE], ontoFileName[STRING_SIZE], outputFileName[STRING_SIZE], option[256], *out, format='s', **name, type='b';
	FILE *fi, *fo;
	int i, j, t, s, niter = 100;
	TypeAnnotation *an;
	TypePartition part, perm;
	double *max, threshold = 0.001;
	TypeStatAnnotationPartition *sap;
	TypeStatFisher sf;
	TypeSignificantTable sig;
	TypeOntologyInfo *info;
	
	info = NULL;
	for(i=0; i<256; i++)
		option[i] = 0;
	   
//	sprintf(outputFileName, "%s.%s", NAME_OUTPUT, EXT_OUTPUT);
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &niter) == 1)
				i++;
			else
				exitProg(ErrorArgument, "an integer number is required after option -i");
		}
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &format) == 1)
				i++;
			else
				exitProg(ErrorArgument, "'c' or 'f' are expected after option -f");
		}
		if(option['o']) {
			option['o'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%s", ontoFileName) == 1) {
				i++;
				if((fi = fopen(ontoFileName, "r"))) {
					info = readOntologyInfo(fi);
					fclose(fi);
				} else
					exitProg(ErrorReading, ontoFileName);
			} else
				exitProg(ErrorArgument, "'c' or 'f' are expected after option -f");
		}
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &type) == 1)
				i++;
			else
				exitProg(ErrorArgument, "'c' or 'f' are expected after option -f");
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &threshold) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a float number is required after option -t");
		}
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}

	
	if(i>=argc)
		exitProg(ErrorArgument, "wrong or absent file name");
printf("\n\nReading Annotation file\n");
	if((fi = fopen(argv[i], "r"))) {
		an = readAnnotation(fi);
		fclose(fi);
	} else
		exitProg(ErrorReading, argv[i]);
	i++;
printf("done\n");
printf("\n\nReading Partition file\n");
	if((fi = fopen(argv[i], "r"))) {
		switch(format) {
			case 'c':
				part = readPartition(fi, &name);
				break;
			case 's':
			default:
				part = readPartitionLine(fi, &name);
		}
		fclose(fi);
	} else
		exitProg(ErrorReading, argv[i]);
	i++;
printf("done\n");
printf("\n\nIntegrating data\n");
	sap = newStatAnnotationPartition(an, &part, name);
printf("done\n");
printf("\n\nComputing Correction\n");
	switch(type) {
		case 'h':
			sig = getHochberg(sap);
			break;
		case 'b':
		default:
			sig = getSignificantTableBonferroni(sap, threshold);
			//getBonferroni(sap);
	}
printf("done\n");
	if(i<argc) {
		out = argv[i];
	} else {
		out = "out.csv";
	}
	if((fo = fopen(out, "w"))) {
		fprintSignificantTableEmpiricalAnais(fo, sig, sap->sizeG, sap->nameA, info);
//		fprintSignificantTableEmpirical(fo, sig, sap->nameA);
//		fprintAnnotation(fo, an);
		fclose(fo);
	} else
		exitProg(ErrorWriting, out);
	freeAnnotation(an);
	freeOntologyInfo(info);
	return 0;
}

