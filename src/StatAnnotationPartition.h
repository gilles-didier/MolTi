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




#ifndef StatAnnotationPartitionF
#define StatAnnotationPartitionF

#include <stdlib.h>
#include <stdio.h>

#include "Annotation.h"
#include "Partition.h"

#define INC_SIZE 100
#define MAX_CURRENT 1000000

#ifdef __cplusplus
extern "C" {
#endif

typedef struct TYPE_STAT_ANNOTATION_PARTITION {
	char **nameG, **nameA;
	int sizeG, sizeA, sizeC, **mat, *atom, *startC, *itemC, *startG, *itemG, *startA, *itemA;
} TypeStatAnnotationPartition;

typedef struct TYPE_SIGNIFICANT {
	int atom, ontology, effC, effO, effCO;
	double fisher, corrected;
} TypeSignificant;

typedef struct TYPE_SIGNIFICANT_TABLE {
	int size, total, *start;
	TypeSignificant *table;
} TypeSignificantTable;

typedef struct TYPE_STAT_FISHER {
	int size;
	double *fisher, *corrected;
} TypeStatFisher;

TypeStatAnnotationPartition *newStatAnnotationPartition(TypeAnnotation *annot, TypePartition *part, char **name);
int *maxAnnotation(TypeStatAnnotationPartition *sap);
double *maxPartition(TypeStatAnnotationPartition *sap);
void freeStatAnnotationPartition(TypeStatAnnotationPartition *sap);
double getFisher(int n1, int n, int k, int t);
double *maxFisher(TypeStatAnnotationPartition *sap);
TypeStatFisher statFisher(TypeStatAnnotationPartition *sap, int niter);
double minFisher(TypeStatAnnotationPartition *sap);
TypeSignificantTable getSignificantTable(TypeStatAnnotationPartition *sap, double threshold);
TypeSignificantTable getSignificantTableEmpirical(TypeStatAnnotationPartition *sap, double threshold, int niter);
void fprintSignificantTableEmpirical(FILE *f, TypeSignificantTable sig, char **name);
void testPartition(int *part, int size);
void removeClass(int cl, TypeStatAnnotationPartition *sap);
TypeSignificantTable getSignificantTableEmpiricalBis(TypeStatAnnotationPartition *sap, double threshold, int niter);
TypeSignificantTable getSignificantTableEmpiricalTer(TypeStatAnnotationPartition *sap, double threshold, int niter);
double minFisherBis(TypeStatAnnotationPartition *sap, int *which);
TypeSignificantTable getSignificantTableEmpiricalX(TypeStatAnnotationPartition *sap, double threshold, int niter);
double getFisherClass(TypeStatAnnotationPartition *sap, int i);
TypeSignificantTable getBonferroni(TypeStatAnnotationPartition *sap);
double *allFisher(TypeStatAnnotationPartition *sap);
TypeSignificantTable getHochberg(TypeStatAnnotationPartition *sap);
TypeSignificantTable getSignificantTableBonferroni(TypeStatAnnotationPartition *sap, double threshold);
void fprintSignificantTableEmpiricalAnais(FILE *f, TypeSignificantTable sig, int totalGene, char **name, TypeOntologyInfo *info);
void freeStatAnnotationPartition(TypeStatAnnotationPartition *sap);
void fprintSignificantTableEmpiricalAnaisBis(FILE *f, TypeSignificantTable sig, int totalGene, char **name, TypeOntologyInfo *info, double threshold);

#ifdef __cplusplus
}
#endif

#endif
