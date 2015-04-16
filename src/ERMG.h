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




#ifndef ERMGF
#define ERMGF

#include <stdlib.h>
#include <stdio.h>
#include "Graph.h"
#include "Partition.h"

typedef struct MODEL_ERMG {
	int sizeAtom;
	double *alpha, **pi;
} TypeERMG;

typedef struct MULTI_MODEL_ERMG {
	int sizeAtom, sizeTable;
	double *alpha, ***pi;
} TypeMultiERMG;

/*return a random ermg model*/
TypeERMG *getEmptyERMG(int nclass);
/*return a random ermg model*/
TypeMultiERMG *getEmptyMultiERMG(int sizeTable, int sizeAtom);
/*return a random ermg model*/
TypeMultiERMG *getBasicMultiERMGX(int nclass, int nTable, double *alpha, double *intra, double *inter, int p);
/*return a random ermg model*/
void freeMultiERMG(TypeMultiERMG *model);
/*free a random ermg model*/
void freeERMG(TypeERMG *model);
/*return a random ermg model*/
TypeMultiERMG *copyMultiERMG(TypeMultiERMG *src);
TypeMultiERMG *copySpecialMultiERMG(TypeMultiERMG *src, int sizeT);
/*return a random ermg model*/
//TypeMultiERMG *getRandomMultiERMG(int nTable, int nClass, double intmin, double intmax, double extmin, double extmax);
TypeMultiERMG *getRandomMultiERMG(int nTable, int nClass, double *alpha, double intmin, double intmax, double extmin, double extmax);
/*return a random ermg model*/
TypeERMG *copyERMG(TypeERMG *src);
/*return a random ermg model*/
TypeERMG *getRandomERMG(int nclass, double *alpha, double intmin, double intmax, double extmin, double extmax);
/*print a multiERMG model*/
void fprintMultiERMG(FILE *f, TypeMultiERMG *model);
/*read a multi ERMG model*/
TypeMultiERMG *freadMultiERMG(FILE *f);
/*read an ERMG model*/
TypeERMG *freadERMG(FILE *f);
/*print an ERMG model in a more compact way*/
void fprintERMGCompact(FILE *f, TypeERMG *model);
/*return a basic ermg model*/
TypeERMG *getBasicERMG(int nclass, double *alpha, double intra, double inter);
void scaleERGM(TypeERMG *model, double scale);
TypeERMG *getRandomSpecialERMG(int nclass, double *alpha, double intmin, double intmax);
double distERMG(TypeERMG *model0, TypeERMG *model1);
/*return a basic ermg model*/
TypeMultiERMG *getBasicMultiERMG(int nclass, int nTable, double *alpha, double intra, double inter);
TypeMultiERMG *getCopyMultiERMG(TypeERMG *muni, int nTable);

#define PROBIN(x,p) ((x)>0?(p):(1-p))
#define LOG_PROBIN(x,p) ((x)>0?(log(p)):(log(1-p)))

#endif
