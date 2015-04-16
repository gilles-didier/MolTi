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




#include <string.h>
#include <math.h>
#include "ERMG.h"
#include "Utils.h"

#define MAX_NAME_SIZE 3000




double distERMG(TypeERMG *model0, TypeERMG *model1) {
	int q, l;
	double dist = 0.;
	for(q=0; q<model0->sizeAtom; q++)
		for(l=0; l<model0->sizeAtom; l++)
			dist += pow(model0->pi[q][l]-model1->pi[q][l], 2.);
	return sqrt(dist);
}

void scaleERGM(TypeERMG *model, double scale) {
	int i, j;
	for(i=0; i<model->sizeAtom; i++)
		for(j=0; j<=i; j++) {
			model->pi[i][j] *= scale;
			model->pi[j][i] *= scale;
		}
}

/*return a random ermg model*/
TypeERMG *getEmptyERMG(int nclass) {
	TypeERMG *model;
	int i,j, k;
	model = (TypeERMG*) malloc(sizeof(TypeERMG));
	model->sizeAtom =nclass;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	model->pi = (double**) malloc(model->sizeAtom*sizeof(double*));
	for(i=0; i<model->sizeAtom; i++) {
		model->pi[i] = (double*) malloc(model->sizeAtom*sizeof(double));
		for(j=0; j<model->sizeAtom; j++)
			model->pi[i][j] = 0.;
	}
	return model;
}

/*return a random ermg model*/
TypeMultiERMG *getEmptyMultiERMG(int sizeTable, int sizeAtom) {
	TypeMultiERMG *model;
	int i,j, t;
	model = (TypeMultiERMG*) malloc(sizeof(TypeMultiERMG));
	model->sizeAtom =sizeAtom;
	model->sizeTable =sizeTable;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	model->pi = (double***) malloc(model->sizeTable*sizeof(double**));
	for(t=0; t<model->sizeTable; t++) {
		model->pi[t] = (double**) malloc(model->sizeAtom*sizeof(double*));
		for(i=0; i<model->sizeAtom; i++) {
			model->pi[t][i] = (double*) malloc(model->sizeAtom*sizeof(double));
			for(j=0; j<model->sizeAtom; j++)
				model->pi[t][i][j] = 0.;
		}
	}
	return model;
}

/*free a random ermg model*/
void freeERMG(TypeERMG *model) {
	int i;
	for(i=0; i<model->sizeAtom; i++)
		free((void*)model->pi[i]);
	free((void*)model->pi);
	free((void*)model->alpha);
	free((void*)model);
}

/*free a random ermg model*/
void freeMultiERMG(TypeMultiERMG *model) {
	int t;
	for(t=0; t<model->sizeTable; t++) {
		int i;
		for(i=0; i<model->sizeAtom; i++)
			free((void*)model->pi[t][i]);
		free((void*)model->pi[t]);
	}
	free((void*)model->pi);
	free((void*)model->alpha);
	free((void*)model);
}

/*return a random ermg model*/
TypeMultiERMG *copySpecialMultiERMG(TypeMultiERMG *src, int sizeT) {
	TypeMultiERMG *model;
	int i,j, t;
	model = (TypeMultiERMG*) malloc(sizeof(TypeMultiERMG));
	model->sizeAtom = src->sizeAtom;
	model->sizeTable = sizeT;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	for(i=0; i<model->sizeAtom; i++)
		model->alpha[i] = src->alpha[i];
	model->pi = (double***) malloc(model->sizeTable*sizeof(double**));
	for(t=0; t<model->sizeTable; t++) {
		model->pi[t] = (double**) malloc(model->sizeAtom*sizeof(double*));
		for(i=0; i<model->sizeAtom; i++) {
			model->pi[t][i] = (double*) malloc(model->sizeAtom*sizeof(double));
			for(j=0; j<model->sizeAtom; j++) 
				model->pi[t][i][j] = src->pi[t][i][j];
		}
	}
	return model;
}


/*return a random ermg model*/
TypeMultiERMG *copyMultiERMG(TypeMultiERMG *src) {
	TypeMultiERMG *model;
	int i,j, t;
	model = (TypeMultiERMG*) malloc(sizeof(TypeMultiERMG));
	model->sizeAtom =src->sizeAtom;
	model->sizeTable =src->sizeTable;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	for(i=0; i<model->sizeAtom; i++)
		model->alpha[i] = src->alpha[i];
	model->pi = (double***) malloc(model->sizeTable*sizeof(double**));
	for(t=0; t<model->sizeTable; t++) {
		model->pi[t] = (double**) malloc(model->sizeAtom*sizeof(double*));
		for(i=0; i<model->sizeAtom; i++) {
			model->pi[t][i] = (double*) malloc(model->sizeAtom*sizeof(double));
			for(j=0; j<model->sizeAtom; j++) 
				model->pi[t][i][j] = src->pi[t][i][j];
		}
	}
	return model;
}

/*return a random ermg model*/
TypeMultiERMG *getRandomMultiERMG(int nTable, int nClass, double *alpha, double intmin, double intmax, double extmin, double extmax) {
	TypeMultiERMG *model;
	double sumAlpha = 0.;
	int i,j, k, t;
	model = (TypeMultiERMG*) malloc(sizeof(TypeMultiERMG));
	model->sizeAtom =nClass;
	model->sizeTable =nTable;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	if(alpha == NULL) {
		for(i=0; i<model->sizeAtom; i++) {
			model->alpha[i] = (UNIF_RAND+1.)/2.;
			sumAlpha += model->alpha[i];
		}
		for(i=0; i<model->sizeAtom; i++)
			model->alpha[i] /= sumAlpha;
	} else
		for(i=0; i<model->sizeAtom; i++)
			model->alpha[i] = alpha[i];
	model->pi = (double***) malloc(model->sizeTable*sizeof(double**));
	for(t=0; t<model->sizeTable; t++) {
		model->pi[t] = (double**) malloc(model->sizeAtom*sizeof(double*));
		for(i=0; i<model->sizeAtom; i++) {
			model->pi[t][i] = (double*) malloc(model->sizeAtom*sizeof(double));
			model->pi[t][i][i] = intmin+UNIF_RAND*(intmax-intmin);
			for(j=0; j<i; j++) {
				model->pi[t][i][j] = extmin+UNIF_RAND*(extmax-extmin);
				model->pi[t][j][i] = model->pi[t][i][j];
			}
		}
	}
	return model;
}
/*return a basic ermg model*/
TypeMultiERMG *getBasicMultiERMG(int nclass, int nTable, double *alpha, double intra, double inter) {
	TypeMultiERMG *model;
	double sumAlpha = 0.;
	int i,j,k,t;
	model = (TypeMultiERMG*) malloc(sizeof(TypeMultiERMG));
	model->sizeAtom =nclass;
	model->sizeTable =nTable;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	if(alpha == NULL)
		for(i=0; i<model->sizeAtom; i++)
			model->alpha[i] = 1./((double)nclass);
	else
		for(i=0; i<model->sizeAtom; i++)
			model->alpha[i] = alpha[i];
	model->pi = (double***) malloc(model->sizeTable*sizeof(double**));
	for(t=0; t<model->sizeTable; t++) {
		model->pi[t] = (double**) malloc(model->sizeAtom*sizeof(double*));
		for(i=0; i<model->sizeAtom; i++) {
			model->pi[t][i] = (double*) malloc(model->sizeAtom*sizeof(double));
			model->pi[t][i][i] = intra;
			for(j=0; j<i; j++) {
				model->pi[t][i][j] = inter;
				model->pi[t][j][i] = model->pi[t][i][j];
			}
		}
	}
	return model;
}

/*return a basic ermg model*/
TypeMultiERMG *getBasicMultiERMGX(int nclass, int nTable, double *alpha, double *intra, double *inter, int p) {
	TypeMultiERMG *model;
	double sumAlpha = 0.;
	int i,j,k,t;
	model = (TypeMultiERMG*) malloc(sizeof(TypeMultiERMG));
	model->sizeAtom =nclass;
	model->sizeTable =nTable;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	if(alpha == NULL)
		for(i=0; i<model->sizeAtom; i++)
			model->alpha[i] = 1./((double)nclass);
	else
		for(i=0; i<model->sizeAtom; i++)
			model->alpha[i] = alpha[i];
	model->pi = (double***) malloc(model->sizeTable*sizeof(double**));
	for(t=0; t<model->sizeTable; t++) {
		int l = RANGE_RAND(p);
		model->pi[t] = (double**) malloc(model->sizeAtom*sizeof(double*));
		for(i=0; i<model->sizeAtom; i++) {
			model->pi[t][i] = (double*) malloc(model->sizeAtom*sizeof(double));
			model->pi[t][i][i] = intra[l];
			for(j=0; j<i; j++) {
				model->pi[t][i][j] = inter[l];
				model->pi[t][j][i] = model->pi[t][i][j];
			}
		}
	}
	return model;
}

/*return a copy ermg model*/
TypeMultiERMG *getCopyMultiERMG(TypeERMG *muni, int nTable) {
	TypeMultiERMG *model;
	double sumAlpha = 0.;
	int i,j,k,t;
	model = (TypeMultiERMG*) malloc(sizeof(TypeMultiERMG));
	model->sizeAtom =muni->sizeAtom;
	model->sizeTable =nTable;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	for(i=0; i<model->sizeAtom; i++)
		model->alpha[i] = muni->alpha[i];
	model->pi = (double***) malloc(model->sizeTable*sizeof(double**));
	for(t=0; t<model->sizeTable; t++) {
		model->pi[t] = (double**) malloc(model->sizeAtom*sizeof(double*));
		for(i=0; i<model->sizeAtom; i++) {
			model->pi[t][i] = (double*) malloc(model->sizeAtom*sizeof(double));
			for(j=0; j<model->sizeAtom; j++) {
				model->pi[t][i][j] = muni->pi[i][j];
			}
		}
	}
	return model;
}

/*return a random ermg model*/
/*
TypeMultiERMG *getRandomMultiERMG(int nTable, int nClass, double intmin, double intmax, double extmin, double extmax) {
	TypeMultiERMG *model;
	double sumAlpha = 0.;
	int i,j, k, t;
	model = (TypeMultiERMG*) malloc(sizeof(TypeMultiERMG));
	model->sizeAtom =nClass;
	model->sizeTable =nTable;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	for(i=0; i<model->sizeAtom; i++) {
		model->alpha[i] = UNIF_RAND;
		sumAlpha += model->alpha[i];
	}
	for(i=0; i<model->sizeAtom; i++)
		model->alpha[i] /= sumAlpha;
	model->pi = (double***) malloc(model->sizeTable*sizeof(double**));
	for(t=0; t<model->sizeTable; t++) {
		model->pi[t] = (double**) malloc(model->sizeAtom*sizeof(double*));
		for(i=0; i<model->sizeAtom; i++) {
			model->pi[t][i] = (double*) malloc(model->sizeAtom*sizeof(double));
			model->pi[t][i][i] = intmin+UNIF_RAND*(intmax-intmin);
			for(j=0; j<i; j++) {
				model->pi[t][i][j] = extmin+UNIF_RAND*(extmax-extmin);
				model->pi[t][j][i] = model->pi[t][i][j];
			}
		}
	}
	return model;
}
*/

/*return a random ermg model*/
TypeERMG *copyERMG(TypeERMG *src) {
	TypeERMG *model;
	int i,j;
	model = (TypeERMG*) malloc(sizeof(TypeERMG));
	model->sizeAtom =src->sizeAtom;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	model->pi = (double**) malloc(model->sizeAtom*sizeof(double*));
	for(i=0; i<model->sizeAtom; i++) {
		model->alpha[i] = src->alpha[i];
		model->pi[i] = (double*) malloc(model->sizeAtom*sizeof(double));
		for(j=0; j<model->sizeAtom; j++) 
			model->pi[i][j] = src->pi[i][j];
	}
	return model;
}

/*return a random ermg model*/
TypeERMG *getRandomERMG(int nclass, double *alpha, double intmin, double intmax, double extmin, double extmax) {
	TypeERMG *model;
	double sumAlpha = 0.;
	int i,j, k;
	model = (TypeERMG*) malloc(sizeof(TypeERMG));
	model->sizeAtom =nclass;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	if(alpha == NULL) {
		for(i=0; i<model->sizeAtom; i++) {
			model->alpha[i] = (UNIF_RAND+1.)/2.;
			sumAlpha += model->alpha[i];
		}
		for(i=0; i<model->sizeAtom; i++)
			model->alpha[i] /= sumAlpha;
	} else
		for(i=0; i<model->sizeAtom; i++)
			model->alpha[i] = alpha[i];
	model->pi = (double**) malloc(model->sizeAtom*sizeof(double*));
	for(i=0; i<model->sizeAtom; i++) {
		model->pi[i] = (double*) malloc(model->sizeAtom*sizeof(double));
		model->pi[i][i] = intmin+UNIF_RAND*(intmax-intmin);
		for(j=0; j<i; j++) {
			model->pi[i][j] = extmin+UNIF_RAND*(extmax-extmin);
			model->pi[j][i] = model->pi[i][j];
		}
	}
	return model;
}

/*return a random ermg model*/
TypeERMG *getRandomSpecialERMG(int nclass, double *alpha, double intmin, double intmax) {
	TypeERMG *model;
	double sumAlpha = 0.;
	int i,j, k;
	model = (TypeERMG*) malloc(sizeof(TypeERMG));
	model->sizeAtom =nclass;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	if(alpha == NULL) {
		for(i=0; i<model->sizeAtom; i++) {
			model->alpha[i] = (UNIF_RAND+1.)/2.;
			sumAlpha += model->alpha[i];
		}
		for(i=0; i<model->sizeAtom; i++)
			model->alpha[i] /= sumAlpha;
	} else
		for(i=0; i<model->sizeAtom; i++)
			model->alpha[i] = alpha[i];
	model->pi = (double**) malloc(model->sizeAtom*sizeof(double*));
	for(i=0; i<model->sizeAtom; i++) {
		model->pi[i] = (double*) malloc(model->sizeAtom*sizeof(double));
		model->pi[i][i] = intmin+UNIF_RAND*(intmax-intmin);
		for(j=0; j<i; j++) {
			model->pi[i][j] = 0.;
			model->pi[j][i] = 0.;
		}
	}
	return model;
}

/*return a basic ermg model*/
TypeERMG *getBasicERMG(int nclass, double *alpha, double intra, double inter) {
	TypeERMG *model;
	double sumAlpha = 0.;
	int i,j, k;
	model = (TypeERMG*) malloc(sizeof(TypeERMG));
	model->sizeAtom =nclass;
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	if(alpha == NULL)
		for(i=0; i<model->sizeAtom; i++)
			model->alpha[i] = 1./((double)nclass);
	else
		for(i=0; i<model->sizeAtom; i++)
			model->alpha[i] = alpha[i];
	model->pi = (double**) malloc(model->sizeAtom*sizeof(double*));
	for(i=0; i<model->sizeAtom; i++) {
		model->pi[i] = (double*) malloc(model->sizeAtom*sizeof(double));
		model->pi[i][i] = intra;
		for(j=0; j<i; j++) {
			model->pi[i][j] = inter;
			model->pi[j][i] = model->pi[i][j];
		}
	}
	return model;
}

/*print an ERMG model in a more compact way*/
void fprintERMGCompact(FILE *f, TypeERMG *model) {
	int i,j;
	if(model->sizeAtom == 0)
		return;
	fprintf(f, "%lf", model->alpha[0]);
	for(i=1; i<model->sizeAtom; i++)
		fprintf(f, "\t%lf", model->alpha[i]);
	fprintf(f, "\n");
	for(i=0; i<model->sizeAtom; i++) {
		fprintf(f, "%lf", model->pi[i][0]);
		for(j=1; j<model->sizeAtom; j++)
			fprintf(f, "\t%lf", model->pi[i][j]);
		fprintf(f, "\n");
	}
}
/*print an ERMG model*/
void fprintERMG(FILE *f, TypeERMG *model) {
	int i,j;
	fprintf(f, "%d\n\n", model->sizeAtom);
	if(model->sizeAtom == 0)
		return;
	fprintf(f, "%lf", model->alpha[0]);
	for(i=1; i<model->sizeAtom; i++)
		fprintf(f, "\t%lf", model->alpha[i]);
	fprintf(f, "\n\n");
	for(i=0; i<model->sizeAtom; i++) {
		fprintf(f, "%lf", model->pi[i][0]);
		for(j=1; j<model->sizeAtom; j++)
			fprintf(f, "\t%lf", model->pi[i][j]);
		fprintf(f, "\n");
	}
}

/*print a multiERMG model*/
void fprintMultiERMG(FILE *f, TypeMultiERMG *model) {
	int i,j, t;
	fprintf(f, "%d\t%d\n\n", model->sizeTable, model->sizeAtom);
	if(model->sizeAtom == 0)
		return;
	fprintf(f, "%lf", model->alpha[0]);
	for(i=1; i<model->sizeAtom; i++)
		fprintf(f, "\t%lf", model->alpha[i]);
	for(t=0; t<model->sizeTable; t++) {
		fprintf(f, "\n\n");
		for(i=0; i<model->sizeAtom; i++) {
			fprintf(f, "%lf", model->pi[t][i][0]);
			for(j=1; j<model->sizeAtom; j++)
				fprintf(f, "\t%lf", model->pi[t][i][j]);
			fprintf(f, "\n");
		}
	}
}

/*read a multi ERMG model*/
TypeMultiERMG *freadMultiERMG(FILE *f) {
	TypeMultiERMG *model;
	char c, tmp[MAX_NAME_SIZE+1];
	int i, j, k, t;
	model = (TypeMultiERMG*) malloc(sizeof(TypeMultiERMG));
	for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
	for(k=0; k<MAX_NAME_SIZE && c != EOF  && !issepline(c); k++) {
		tmp[k] = c;
		c = fgetc(f);
	}
	tmp[k++] = '\0';
	model->sizeTable = atoi(tmp);
	for(c = fgetc(f); c != EOF && issep(c); c = fgetc(f));
	for(k=0; k<MAX_NAME_SIZE && c != EOF  && !issepline(c); k++) {
		tmp[k] = c;
		c = fgetc(f);
	}
	tmp[k++] = '\0';
	model->sizeAtom = atoi(tmp);
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	model->pi = (double***) malloc(model->sizeTable*sizeof(double*));
	for(t=0; t<model->sizeTable; t++) {
		model->pi[t] = (double**) malloc(model->sizeAtom*sizeof(double*));
		for(i=0; i<model->sizeAtom; i++)
			model->pi[t][i] = (double*) malloc(model->sizeAtom*sizeof(double));
	}
	for(; c != EOF && issepline(c); c = fgetc(f));
	i = 0;
	while( c != EOF && !isline(c)) {
		for(k=0; k<MAX_NAME_SIZE && c != EOF && !issepline(c); k++) {
			tmp[k] = c;
			c = fgetc(f);
		}
		tmp[k++] = '\0';
		model->alpha[i++] = atof(tmp);
		for(; c != EOF && issep(c); c = fgetc(f));
	}
	if(i!= model->sizeAtom)
		exitProg(ErrorReading, "Not enough alphas");
	for(; c != EOF && issepline(c); c = fgetc(f));
	t = 0;
	while(c != EOF && t<model->sizeTable) {
		i = 0;
		while(c != EOF && !isline(c) && i<model->sizeAtom) {
			j = 0;
			while(!isline(c)) {
				for(k=0; k<MAX_NAME_SIZE && c != EOF &&  !issepline(c); k++) {
					tmp[k] = c;
					c = fgetc(f);
				}
				tmp[k++] = '\0';
				model->pi[t][i][j++] = atof(tmp);
				for(; c != EOF && issep(c); c = fgetc(f));
			}
			if(j!= model->sizeAtom)
				exitProg(ErrorReading, "Missing pi");
			i++;
			for(; c != EOF && issepline(c); c = fgetc(f));
		}
		if(i!= model->sizeAtom)
			exitProg(ErrorReading, "Not enough pi");
	}
	if(t!= model->sizeTable)
		exitProg(ErrorReading, "Not enough pi");
	return model;
}

/*read an ERMG model*/
TypeERMG *freadERMG(FILE *f) {
	TypeERMG *model;
	char c, tmp[MAX_NAME_SIZE+1];
	int i,j, k;
	model = (TypeERMG*) malloc(sizeof(TypeERMG));
	for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
	for(k=0; k<MAX_NAME_SIZE && c != EOF  && !issepline(c); k++) {
		tmp[k] = c;
		c = fgetc(f);
	}
	tmp[k++] = '\0';
	model->sizeAtom = atoi(tmp);
	model->alpha = (double*) malloc(model->sizeAtom*sizeof(double));
	model->pi = (double**) malloc(model->sizeAtom*sizeof(double*));
	for(i=0; i<model->sizeAtom; i++)
		model->pi[i] = (double*) malloc(model->sizeAtom*sizeof(double));
	i = 0;
	for(; c != EOF && issepline(c); c = fgetc(f));
	while(!isline(c) && i<model->sizeAtom) {
		for(k=0; k<MAX_NAME_SIZE && c != EOF && !issepline(c); k++) {
			tmp[k] = c;
			c = fgetc(f);
		}
		tmp[k++] = '\0';
		model->alpha[i++] = atof(tmp);
		for(; c != EOF && issep(c); c = fgetc(f));
	}
	if(i != model->sizeAtom)
		exitProg(ErrorReading, "Not enough alphas");
	for(; c != EOF && issepline(c); c = fgetc(f));
	i = 0;
	while(c != EOF && i<model->sizeAtom) {
		j = 0;
		while(!isline(c)) {
			for(k=0; k<MAX_NAME_SIZE && c != EOF &&  !issepline(c); k++) {
				tmp[k] = c;
				c = fgetc(f);
			}
			tmp[k++] = '\0';
			model->pi[i][j++] = atof(tmp);
			for(; c != EOF && issep(c); c = fgetc(f));
		}
		if(j!= model->sizeAtom)
			exitProg(ErrorReading, "Missing pi");
		i++;
		for(; c != EOF && issepline(c); c = fgetc(f));
	}
	if(i!= model->sizeAtom)
		exitProg(ErrorReading, "Not enough pi");
	return model;
}
