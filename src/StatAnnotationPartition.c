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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>

#include "Utils.h"
#include "StatAnnotationPartition.h"

static double getMin(double *tab, int size);
static void shuffle(const gsl_rng * r, int *elt, int size);
static int compareSignificant(const void* a, const void* b);


TypeStatAnnotationPartition *newStatAnnotationPartition(TypeAnnotation *annot, TypePartition *part, char **name) {
    int *indexA, *indexG, *index, i, j, size, n, c, *start, *next, ind, sum;
    TypeLexiTree *dict;
    TypeStatAnnotationPartition *sap;

    dict = getDictFromTable(annot->nameG, annot->sizeG);

    indexA = (int*) malloc(annot->sizeG*sizeof(int));
    indexG = (int*) malloc(part->sizeItem*sizeof(int));
    index = (int*) malloc(UTILS_MIN(annot->sizeG, part->sizeItem)*sizeof(int));
    size = 0;
    for(i=0; i<part->sizeItem; i++)
        indexG[i] = -1;
    for(i=0; i<annot->sizeG; i++)
        indexA[i] = -1;
    for(i=0; i<part->sizeItem; i++)
        if((n = findWordLexi(name[i], dict)) >= 0) {
            index[size] = n;
            indexA[n] = size;
            indexG[i] = size;
            size++;
        }
    sap = (TypeStatAnnotationPartition*) malloc(sizeof(TypeStatAnnotationPartition));
    sap->sizeG = size;
    sap->sizeA = annot->sizeA;
    sap->sizeC = part->sizeAtom;
    sap->atom = (int*) malloc(sap->sizeG*sizeof(int));
    next = (int*) malloc(sap->sizeG*sizeof(int));
    start = (int*) malloc(sap->sizeC*sizeof(int));
    sap->startC = (int*) malloc((sap->sizeC+1)*sizeof(int*));
    sap->itemC = (int*) malloc((sap->sizeG)*sizeof(int*));
    sap->mat = (int**) malloc(sap->sizeG*sizeof(int*));
    sap->nameG = (char**) malloc(sap->sizeG*sizeof(char*));
    sap->nameA = (char**) malloc(sap->sizeA*sizeof(char*));
    for(i=0; i<sap->sizeG; i++)
        sap->mat[i] = (int*) malloc(sap->sizeA*sizeof(int));
    for(i=0; i<annot->sizeG; i++)
        if(indexA[i] != -1) {
            sap->nameG[indexA[i]] = strdup(annot->nameG[index[indexA[i]]]);
            sap->mat[indexA[i]] = (int*) malloc(sap->sizeA*sizeof(int));
            for(j=0; j<sap->sizeA; j++)
                sap->mat[indexA[i]][j] = annot->mat[i][j];
        }
    for(i=0; i<part->sizeItem; i++)
        if(indexG[i] != -1)
            sap->atom[indexG[i]] = part->atom[i];

    for(c=0; c<sap->sizeC; c++)
        start[c] = -1;
    for(i=0; i<sap->sizeG; i++) {
        next[i] = start[sap->atom[i]];
        start[sap->atom[i]] = i;
    }
    ind = 0;
    for(c=0; c<sap->sizeC; c++) {
        sap->startC[c] = ind;
        for(i=start[c]; i>=0; i=next[i])
            sap->itemC[ind++] = i;
    }
    sap->startC[sap->sizeC] = ind;
    for(i=0; i<annot->sizeA; i++) {
        sap->nameA[i] = (char*) malloc((strlen(annot->nameA[i])+1)*sizeof(char));
        strcpy(sap->nameA[i], annot->nameA[i]);
    }
    sum = 0;
    for(i=0; i<sap->sizeG; i++)
        for(j=0; j<sap->sizeA; j++)
            if(sap->mat[i][j])
                sum++;

    sap->startG = (int*) malloc((sap->sizeG+1)*sizeof(int*));
    sap->itemG = (int*) malloc(sum*sizeof(int*));
    ind = 0;
    for(i=0; i<sap->sizeG; i++) {
        sap->startG[i] = ind;
        for(j=0; j<sap->sizeA; j++)
            if(sap->mat[i][j])
                sap->itemG[ind++] = j;
    }
    sap->startG[sap->sizeG] = ind;

    sap->startA = (int*) malloc((sap->sizeA+1)*sizeof(int*));
    sap->itemA = (int*) malloc(sum*sizeof(int*));
    ind = 0;
    for(j=0; j<sap->sizeA; j++) {
        sap->startA[j] = ind;
        for(i=0; i<sap->sizeG; i++)
            if(sap->mat[i][j])
                sap->itemA[ind++] = i;
    }
    sap->startA[sap->sizeA] = ind;
    freeLexiTree(dict);
    free((void*)start);
    free((void*)index);
    free((void*)indexA);
    free((void*)indexG);
    return sap;
}

void testPartition(int *part, int size) {
    int *present, i;

    present = (int*) malloc(size*sizeof(int));
    for(i=0; i<size; i++)
        present[i] = -1;
    for(i=0; i<size; i++)
        if(part[i]>=0 && part[i]<size) {
            if(present[part[i]] >= 0) {
                printf("Duplicate part[%d] = %d (part[%d] = %d) (present = %d)\n", i, part[i], present[part[i]], part[present[part[i]]], present[part[i]]);
                exit(0);
            } else
                present[part[i]] = i;
        } else {
            printf("Overflow part[%d] = %d\n", i, part[i]);
            exit(0);
        }
    free((void*)present);
}

void fprintStatAnnotationPartitionInfo(FILE *f, TypeStatAnnotationPartition *sap) {
    fprintf(f, "%d genes %d annotations %d atomes\n", sap->sizeG, sap->sizeA, sap->sizeC);
}

int *maxAnnotation(TypeStatAnnotationPartition *sap) {
    int *eff, *max, i, j, c;

    max = (int*) malloc(sap->sizeA*sizeof(int));
    eff = (int*) malloc(sap->sizeC*sizeof(int*));
    for(i=0; i<sap->sizeA; i++) {
        for(c=0; c<sap->sizeC; c++)
            eff[c] = 0;
        for(j=0; j<sap->sizeG; j++)
            if(sap->mat[j][i])
                eff[sap->atom[j]]++;
        max[i] = 0;
        for(c=0; c<sap->sizeC; c++)
            if(eff[c]>max[i])
                max[i] = eff[c];
    }
    free((void*)eff);
    return max;
}

double *maxPartition(TypeStatAnnotationPartition *sap) {
    int *eff, i, j, c;
    double *max;

    max = (double*) malloc(sap->sizeC*sizeof(double));
    eff = (int*) malloc(sap->sizeA*sizeof(int*));
    for(i=0; i<sap->sizeC; i++) {
        for(c=0; c<sap->sizeA; c++)
            eff[c] = 0;
        for(j=sap->startC[i]; j<sap->startC[i+1]; j++)
            for(c=0; c<sap->sizeA; c++)
                if(sap->mat[sap->itemC[j]][c])
                    eff[c]++;
        max[i] = 0;
        for(c=0; c<sap->sizeA; c++)
            if(eff[c]>max[i])
                max[i] = (double) eff[c];
        max[i] /= (double) sap->startC[i+1]-sap->startC[i];
    }
    free((void*)eff);
    return max;
}

/*Fisher C(n_1, k) C(n-n_1, t - k) / C(n_1 + n_2, t) n_1 number of elt with ontologie o in the atom, n card of the atom, k card of o, t total*/

double getFisher(int k, int n1, int n2, int t) {
    int i = n1;
    double res = 0.;

    res = gsl_cdf_hypergeometric_Q(k-1, n1, n2, t);
    if(res == 0.)
        printf("Probable overflow : k %d n1 %d n2 %d t %d %lE\n", k, n1, n2, t,res);
    return res;
}

double *maxFisher(TypeStatAnnotationPartition *sap) {
    int *eff, i, j, k, c;
    double *fisher;

    fisher = (double*) malloc(sap->sizeC*sizeof(double));
    eff = (int*) malloc(sap->sizeA*sizeof(int*));
    for(i=0; i<sap->sizeC; i++) {
        for(c=0; c<sap->sizeA; c++)
            eff[c] = 0;
        for(j=sap->startC[i]; j<sap->startC[i+1]; j++)
            for(k=sap->startG[sap->itemC[j]]; k<sap->startG[sap->itemC[j]+1]; k++)
                eff[sap->itemG[k]]++;
        fisher[i] = 1.;
        for(c=0; c<sap->sizeA; c++) {
            if(eff[c]) {
                double tmp = getFisher(eff[c], sap->startA[c+1]-sap->startA[c], sap->sizeG-(sap->startA[c+1]-sap->startA[c]), sap->startC[i+1]-sap->startC[i]);
                if(tmp<fisher[i]) {
                    fisher[i] = tmp;
                }
                if(tmp == 0. || tmp == 1.) {
                    int t0 = 0, t1 = 0;
                    printf("Possible problem on atom %d of size %d (%d to %d) with annotation %d of size %d -> effectif %d fisher %.2lE\n", i, sap->startC[i+1]-sap->startC[i], sap->startC[i], sap->startC[i+1], c, sap->startA[c+1]-sap->startA[c], eff[c], tmp);
                }
            }
        }
    }
    free((void*)eff);
    return fisher;
}


double *allFisher(TypeStatAnnotationPartition *sap) {
    int *eff, i, j, k, c;
    double *fisher;

    fisher = (double*) malloc(sap->sizeC*sap->sizeA*sizeof(double));
    eff = (int*) malloc(sap->sizeA*sizeof(int*));
    for(i=0; i<sap->sizeC; i++) {
        for(c=0; c<sap->sizeA; c++)
            eff[c] = 0;
        for(j=sap->startC[i]; j<sap->startC[i+1]; j++)
            for(k=sap->startG[sap->itemC[j]]; k<sap->startG[sap->itemC[j]+1]; k++)
                eff[sap->itemG[k]]++;
        for(c=0; c<sap->sizeA; c++) {
            fisher[i*sap->sizeA+c] = getFisher(eff[c], sap->startA[c+1]-sap->startA[c], sap->sizeG-(sap->startA[c+1]-sap->startA[c]), sap->startC[i+1]-sap->startC[i]);
        }
    }
    free((void*)eff);
    return fisher;
}

double minFisher(TypeStatAnnotationPartition *sap) {
    int *eff, i, j, k, c;
    double min = 1.;

    eff = (int*) malloc(sap->sizeA*sizeof(int*));
    for(i=0; i<sap->sizeC; i++) {
        for(c=0; c<sap->sizeA; c++)
            eff[c] = 0;
        for(j=sap->startC[i]; j<sap->startC[i+1]; j++)
            for(k=sap->startG[sap->itemC[j]]; k<sap->startG[sap->itemC[j]+1]; k++)
                eff[sap->itemG[k]]++;
        for(c=0; c<sap->sizeA; c++) {
            if(eff[c]) {
//				double tmp = ((double)sap->sizeA)*getFisher(eff[c], sap->startA[c+1]-sap->startA[c], sap->sizeG-(sap->startA[c+1]-sap->startA[c]), sap->startC[i+1]-sap->startC[i]);
                double tmp = getFisher(eff[c], sap->startA[c+1]-sap->startA[c], sap->sizeG-(sap->startA[c+1]-sap->startA[c]), sap->startC[i+1]-sap->startC[i]);
                if(tmp<min) {
                    min = tmp;
                }
                if(tmp == 0.) {
                    int t0 = 0, t1 = 0;
                    printf("Problem: atom %d of size %d (%d to %d) with annotation %d of size %d -> effectif %d\n", i, sap->startC[i+1]-sap->startC[i], sap->startC[i], sap->startC[i+1], c, sap->startA[c+1]-sap->startA[c], eff[c]);
                    for(j=sap->startC[i]; j<sap->startC[i+1]; j++)
                        for(k=sap->startG[sap->itemC[j]]; k<sap->startG[sap->itemC[j]+1]; k++) {
                            if(sap->itemG[k] == c) {
                                printf("gene%d %d annot%d %d\n", j, sap->itemC[j], k, sap->itemG[k]);
                                t0++;
                            }
                        }
                }
            }
        }
    }
    free((void*)eff);
    return min;
}

double getFisherClass(TypeStatAnnotationPartition *sap, int i) {
    int *eff, j, k, c;
    double min = 1.;

    eff = (int*) malloc(sap->sizeA*sizeof(int*));
    for(c=0; c<sap->sizeA; c++)
        eff[c] = 0;
    for(j=sap->startC[i]; j<sap->startC[i+1]; j++)
        for(k=sap->startG[sap->itemC[j]]; k<sap->startG[sap->itemC[j]+1]; k++)
            eff[sap->itemG[k]]++;
    for(c=0; c<sap->sizeA; c++) {
        if(eff[c]) {
            double tmp = getFisher(eff[c], sap->startA[c+1]-sap->startA[c], sap->sizeG-(sap->startA[c+1]-sap->startA[c]), sap->startC[i+1]-sap->startC[i]);
            if(tmp<min) {
                min = tmp;
            }
            if(tmp == 0.) {
                int t0 = 0, t1 = 0;
                printf("Problem: atom %d of size %d (%d to %d) with annotation %d of size %d -> effectif %d\n", i, sap->startC[i+1]-sap->startC[i], sap->startC[i], sap->startC[i+1], c, sap->startA[c+1]-sap->startA[c], eff[c]);
                for(j=sap->startC[i]; j<sap->startC[i+1]; j++)
                    for(k=sap->startG[sap->itemC[j]]; k<sap->startG[sap->itemC[j]+1]; k++) {
                        if(sap->itemG[k] == c) {
                            printf("gene%d %d annot%d %d\n", j, sap->itemC[j], k, sap->itemG[k]);
                            t0++;
                        }
                    }
            }
        }
    }
    free((void*)eff);
    return min;
}

double minFisherBis(TypeStatAnnotationPartition *sap, int *which) {
    int *eff, i, j, k, c;
    double min = 1.;
    *which = -1;
    eff = (int*) malloc(sap->sizeA*sizeof(int*));
    for(i=0; i<sap->sizeC; i++) {
        for(c=0; c<sap->sizeA; c++)
            eff[c] = 0;
        for(j=sap->startC[i]; j<sap->startC[i+1]; j++)
            for(k=sap->startG[sap->itemC[j]]; k<sap->startG[sap->itemC[j]+1]; k++)
                eff[sap->itemG[k]]++;
        for(c=0; c<sap->sizeA; c++) {
            if(eff[c]) {
                double tmp = getFisher(eff[c], sap->startA[c+1]-sap->startA[c], sap->sizeG-(sap->startA[c+1]-sap->startA[c]), sap->startC[i+1]-sap->startC[i]);
                if(tmp<min) {
                    min = tmp;
                    *which = i;
                }
                if(tmp == 0.) {
                    int t0 = 0, t1 = 0;
                    printf("Possible problem: atom %d of size %d (%d to %d) with annotation %d of size %d -> effectif %d\n", i, sap->startC[i+1]-sap->startC[i], sap->startC[i], sap->startC[i+1], c, sap->startA[c+1]-sap->startA[c], eff[c]);
                }
            }
        }
    }
    free((void*)eff);
    return min;
}

int compareSignificant(const void* a, const void* b) {
    return compareDouble(&(((TypeSignificant*)a)->corrected), &(((TypeSignificant*)b)->corrected));
}

TypeSignificantTable getSignificantTableBonferroni(TypeStatAnnotationPartition *sap, double threshold) {
    int *eff, i, j, k, c, inc_buffer = sap->startG[sap->sizeG]/10, size_buffer = inc_buffer, ind;
    TypeSignificantTable res;
    double bonf, *min;
    TypeSignificant *tmp;
    size_t *index;

    res.table = (TypeSignificant*) malloc(size_buffer*sizeof(TypeSignificant));
    res.size = 0;
    res.total = sap->sizeG;
    res.start = (int*) malloc((sap->sizeC+1)*sizeof(int*));
    bonf = ((double) sap->sizeC)*((double) sap->sizeA);
    min = (double*) malloc(sap->sizeC*sizeof(double*));
    eff = (int*) malloc(sap->sizeA*sizeof(int*));
    for(i=0; i<sap->sizeC; i++) {
        for(c=0; c<sap->sizeA; c++)
            eff[c] = 0;
        for(j=sap->startC[i]; j<sap->startC[i+1]; j++)
            for(k=sap->startG[sap->itemC[j]]; k<sap->startG[sap->itemC[j]+1]; k++)
                eff[sap->itemG[k]]++;
        res.start[i] = res.size;
        for(c=0; c<sap->sizeA; c++) {
            if(eff[c]) {
                double tmp = getFisher(eff[c], sap->startA[c+1]-sap->startA[c], sap->sizeG-(sap->startA[c+1]-sap->startA[c]), sap->startC[i+1]-sap->startC[i]), corrected = tmp*bonf;
                if(corrected<=threshold) {
                    if(res.size>=size_buffer) {
                        size_buffer += inc_buffer;
                        res.table = (TypeSignificant*) realloc(res.table, size_buffer*sizeof(TypeSignificant));
                    }
                    res.table[res.size].effC = sap->startC[i+1]-sap->startC[i];
                    res.table[res.size].effO = sap->startA[c+1]-sap->startA[c];
                    res.table[res.size].effCO = eff[c];
                    res.table[res.size].atom = i;
                    res.table[res.size].ontology = c;
                    res.table[res.size].fisher = tmp;
                    res.table[res.size].corrected = corrected;
                    res.size++;
                }
            }
        }
        if(res.size>res.start[i]) {
            qsort(&(res.table[res.start[i]]), res.size-res.start[i], sizeof(TypeSignificant), compareSignificant);
            min[i] = res.table[res.start[i]].corrected;
        } else
            min[i] = 1.;
    }
    res.start[sap->sizeC] = res.size;
    free((void*)min);
    free((void*)eff);
    return res;
}



TypeSignificantTable getSignificantTable(TypeStatAnnotationPartition *sap, double threshold) {
    int *eff, i, j, k, c, inc_buffer = sap->startG[sap->sizeG]/10, size_buffer = inc_buffer;
    TypeSignificantTable res;
    res.table = (TypeSignificant*) malloc(size_buffer*sizeof(TypeSignificant));
    res.size = 0;
    eff = (int*) malloc(sap->sizeA*sizeof(int*));
    for(i=0; i<sap->sizeC; i++) {
        for(c=0; c<sap->sizeA; c++)
            eff[c] = 0;
        for(j=sap->startC[i]; j<sap->startC[i+1]; j++)
            for(k=sap->startG[sap->itemC[j]]; k<sap->startG[sap->itemC[j]+1]; k++)
                eff[sap->itemG[k]]++;
        for(c=0; c<sap->sizeA; c++) {
            if(eff[c]) {
//				double tmp = ((double)sap->sizeA)*getFisher(eff[c], sap->startA[c+1]-sap->startA[c], sap->sizeG-(sap->startA[c+1]-sap->startA[c]), sap->startC[i+1]-sap->startC[i]);
                double tmp = getFisher(eff[c], sap->startA[c+1]-sap->startA[c], sap->sizeG-(sap->startA[c+1]-sap->startA[c]), sap->startC[i+1]-sap->startC[i]);
                if(tmp<threshold) {
                    if(res.size>=size_buffer) {
                        size_buffer += inc_buffer;
                        res.table = (TypeSignificant*) realloc(res.table, size_buffer*sizeof(TypeSignificant));
                    }
                    res.table[res.size].effC = sap->startC[i+1]-sap->startC[i];
                    res.table[res.size].effO = sap->startA[c+1]-sap->startA[c];
                    res.table[res.size].effCO = eff[c];
                    res.table[res.size].atom = i;
                    res.table[res.size].ontology = c;
                    res.table[res.size].fisher = tmp;
                    res.table[res.size].corrected = 0.;
                    res.size++;
                }
            }
        }
    }
    free((void*)eff);
printf("size %d\n", res.size);
    return res;
}

TypeSignificantTable getSignificantTableEmpirical(TypeStatAnnotationPartition *sap, double threshold, int niter) {
    int *eff, i, j, k, c, inc_buffer = sap->startG[sap->sizeG]/10, size_buffer = inc_buffer, ind;
    TypeSignificantTable res = getSignificantTable(sap, 0.0001);
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    for(i=0; i<niter; i++) {
        double min;
        shuffle(r, sap->itemC, sap->sizeG);
        min = minFisher(sap);
        for(j=0; j<res.size; j++)
            if(min<=res.table[j].fisher)
                res.table[j].corrected++;
    }
    gsl_rng_free (r);
    for(j=0; j<res.size; j++)
        res.table[j].corrected /= (double) niter;

    for(j=0; j<res.size; j++)
        res.table[j].corrected /= (double) niter;
    ind = 0;
    for(j=0; j<res.size; j++)
        if(res.table[j].corrected <= threshold)
            res.table[ind++] = res.table[j];
    res.size = ind;
    res.table = (TypeSignificant*) realloc(res.table, res.size*sizeof(TypeSignificant));
    return res;
}

void fprintSignificantTableEmpirical(FILE *f, TypeSignificantTable sig, char **name) {
    int i, j, last, g;

    if(sig.size == 0)
        return;
    last = sig.table[0].atom;
    g = 1;
    fprintf(f, "1st atom %d (size %d)\n", sig.table[0].atom+1, sig.table[0].effC);
    for(i=0; i<sig.size; i++) {
        if(sig.table[i].atom != last) {
            fprintf(f, "\n%dth atom %d (size %d)\n", ++g, sig.table[i].atom+1, sig.table[i].effC);
            last = sig.table[i].atom;
        }
        fprintf(f, "%s\t%.3lE\t%.3lE (%d elts from ont. atom of size %d)\n", name[sig.table[i].ontology], sig.table[i].fisher, sig.table[i].corrected, sig.table[i].effCO, sig.table[i].effO);
    }
}

void fprintSignificantTableEmpiricalAnais(FILE *f, TypeSignificantTable sig, int totalGene, char **name, TypeOntologyInfo *info) {
    int i, j, last, g;
    TypeLexiTree *dict;

    if(sig.size == 0)
        return;
    if(info != NULL)
        dict = 	getDictFromTable(info->id, info->size);
    else
        dict = NULL;
    last = sig.table[0].atom;
    g = 1;
    fprintf(f, "Community %d (annotated elements %d)\n", sig.table[0].atom+1, sig.table[0].effC);
    for(i=0; i<sig.size; i++) {
        int n;
        if(sig.table[i].atom != last) {
            fprintf(f, "\natom %d (size %d)\n", sig.table[i].atom+1, sig.table[i].effC);
            last = sig.table[i].atom;
        }
//		fprintf(f, "%s\t%.3lE\t%.3lE\t(%d %d %d %d)\n", name[sig.table[i].ontology], sig.table[i].fisher, sig.table[i].corrected, sig.table[i].effCO, sig.table[i].effC-sig.table[i].effCO, sig.table[i].effO, totalGene-sig.table[i].effO);
        fprintf(f, "%s\t%.3lE\t%.3lE\t(%d %d %d %d)", name[sig.table[i].ontology], sig.table[i].fisher, sig.table[i].corrected, sig.table[i].effCO, sig.table[i].effC-sig.table[i].effCO, sig.table[i].effO-sig.table[i].effCO, totalGene-sig.table[i].effO-sig.table[i].effC+sig.table[i].effCO);
        if(info != NULL && (n = findWordLexi(name[sig.table[i].ontology], dict)) >= 0)
            fprintf(f, "\t%s", info->name[n]);
        fprintf(f, "\n");
    }
    if(dict != NULL)
        freeLexiTree(dict);
}


void fprintSignificantTableEmpiricalAnaisBis(FILE *f, TypeSignificantTable sig, int totalGene, char **name, TypeOntologyInfo *info, double threshold) {
    int i, j, last, g;
    TypeLexiTree *dict;

    if(sig.size == 0)
        return;
    if(info != NULL)
        dict = 	getDictFromTable(info->id, info->size);
    else
        dict = NULL;
    last = -1;
    for(i=0; i<sig.size; ) {
        int n;
        for(; i<sig.size && sig.table[i].corrected >= threshold; i++);
        if(i<sig.size) {
            fprintf(f, "\natom %d (size %d)\n", sig.table[i].atom, sig.table[i].effC);
            last = sig.table[i].atom;
            for(; i<sig.size && sig.table[i].corrected <= threshold && sig.table[i].atom == last; i++) {
                fprintf(f, "%s\t%.3lE\t%.3lE\t(%d %d %d %d)", name[sig.table[i].ontology], sig.table[i].fisher, sig.table[i].corrected, sig.table[i].effCO, sig.table[i].effC-sig.table[i].effCO, sig.table[i].effO-sig.table[i].effCO, totalGene-sig.table[i].effO-sig.table[i].effC+sig.table[i].effCO);
                if(info != NULL && (n = findWordLexi(name[sig.table[i].ontology], dict)) >= 0)
                    fprintf(f, "\t%s", info->name[n]);
                fprintf(f, "\n");
            }
        }
    }
    if(dict != NULL)
        freeLexiTree(dict);
}

void shuffle(const gsl_rng *r, int *elt, int size) {
    int i, tmp;
    for(i=0; i<size; i++) {
        int tmp;
//		int j = (rand() % (size-i))+i;
        int j = gsl_rng_uniform_int(r,  (size-i))+i;
        tmp = elt[i];
        elt[i] = elt[j];
        elt[j] = tmp;
    }
}

double getMin(double *tab, int size) {
    int i;
    double min = tab[0];
    for(i=1; i<size; i++)
        if(tab[i]<min)
            min = tab[i];
    return min;
}

TypeStatFisher statFisher(TypeStatAnnotationPartition *sap, int niter) {
    double *tmp;
    int i, j, c;
    TypeStatFisher sf;

    sf.corrected = (double*) malloc(sap->sizeC*sizeof(double));
    for(c=0; c<sap->sizeC; c++)
        sf.corrected[c] = 0.;
    sf.fisher = maxFisher(sap);
//	qsort(sf.fisher, sap->sizeC, sizeof(double), compareDouble);
    for(i=0; i<niter; i++) {
        double min;
//		shuffle(sap->itemC, sap->sizeG);
        tmp = maxFisher(sap);
//		qsort(tmp, sap->sizeC, sizeof(double), compareDouble);
        min = getMin(tmp, sap->sizeC);
//printf("min %lE\n", min);
/*		for(c=0; c<sap->sizeC; c++)
            if(tmp[c]<=sf.corrected[c])
                sf.corrected[c] = tmp[c];
*/
        for(c=0; c<sap->sizeC; c++)
            if(min<=sf.fisher[c])
                sf.corrected[c]++;
        free((void*)tmp);
    }
    for(c=0; c<sap->sizeC; c++)
        sf.corrected[c] /= (double) niter;
    return sf;
}

TypeSignificantTable getSignificantTableEmpiricalX(TypeStatAnnotationPartition *sap, double threshold, int niter) {
    int *eff, i, j, k, c, inc_buffer = sap->startG[sap->sizeG]/10, size_buffer = inc_buffer, ind, prec, *itemCtmp, *save;
    double *fisher, *fisherBis, effTmp, fishTmp;
    const gsl_rng_type *T;
    gsl_rng *r;
    TypeSignificantTable res;
    size_t *index;
    res.table = (TypeSignificant*) malloc(sap->sizeC*sizeof(TypeSignificant));
    itemCtmp = (int*) malloc(sap->sizeG*sizeof(int));
    res.size = 0;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);


/*	fisher = maxFisher(sap);
    prec = 0;
    for(j=1; j<sap->sizeC; j++)
        if(fisher[prec]>fisher[j])
            prec = j;
printf("atom min %d\t%.2lE\n", prec, fisher[prec]);
    removeClass(prec, sap);
    fisherBis = maxFisher(sap);
    for(j=0; j<sap->sizeC; j++)
        printf("%d\t%.2lE\t%.2lE\n", j, fisher[j], fisherBis[j]);

    fishTmp = minFisherBis(sap, &j);
printf("%dth atom\t%.2lE\t%.2lE\n", j, fishTmp, effTmp);
exit(0);
*/
    fishTmp = minFisherBis(sap, &j);
printf("%d\t%lE\n", j, fishTmp);
    effTmp = 0.;
    for(k=0; k<sap->sizeG; k++)
        itemCtmp[k] = sap->itemC[k];
    if(fishTmp>0.000000001) {
        save = sap->itemC;
        sap->itemC = itemCtmp;
        for(i=0; i<niter; i++) {
            double min;
            shuffle(r, sap->itemC, sap->sizeG);
            min = minFisher(sap);
            if(min<=fishTmp)
                effTmp++;
        }
        sap->itemC = save;
    }
    effTmp /= (double) niter;
printf("%dth atom %d \t%.2lE\t%.2lE\n", res.size+1, j, fishTmp, effTmp);
    if(effTmp <= threshold) {
        res.table[res.size].corrected = effTmp;
        res.table[res.size].fisher = fishTmp;
        res.table[res.size].atom = j;
        res.table[res.size].ontology = j;
        res.table[res.size].effC = j;
        res.table[res.size].effO = j;
        res.table[res.size].effCO = j;
        res.size++;
    }
    prec = j;
    while(res.size<sap->sizeC && effTmp<=threshold) {
        removeClass(prec, sap);
        fishTmp = minFisherBis(sap, &j);
        effTmp = 0.;
        for(k=0; k<sap->sizeG; k++)
            itemCtmp[k] = sap->itemC[k];
        if(fishTmp>0.000000001) {
            save = sap->itemC;
            sap->itemC = itemCtmp;
            for(i=0; i<niter; i++) {
                double min;
                shuffle(r, sap->itemC, sap->sizeG);
                min = minFisher(sap);
                if(min<=fishTmp)
                    effTmp++;
            }
            sap->itemC = save;
        }
        effTmp /= (double) niter;
printf("%dth atom %d \t%.2lE\t%.2lE\n", res.size+1, j, fishTmp, effTmp);
        if(effTmp <= threshold) {
            res.table[res.size].corrected = effTmp;
            res.table[res.size].fisher = fishTmp;
            res.table[res.size].atom = j;
            res.table[res.size].ontology = j;
            res.table[res.size].effC = j;
            res.table[res.size].effO = j;
            res.table[res.size].effCO = j;
            res.size++;
        }
        prec = j;
    }
    gsl_rng_free (r);
    res.table = (TypeSignificant*) realloc(res.table, res.size*sizeof(TypeSignificant));
    return res;
}

TypeSignificantTable getSignificantTableEmpiricalTer(TypeStatAnnotationPartition *sap, double threshold, int niter) {
    int *eff, i, j, k, c, inc_buffer = sap->startG[sap->sizeG]/10, size_buffer = inc_buffer, ind;
    double *fisher, effTmp, fishTmp;
    const gsl_rng_type *T;
    gsl_rng *r;
    TypeSignificantTable res;
    size_t *index;
    res.table = (TypeSignificant*) malloc(sap->sizeC*sizeof(TypeSignificant));
    res.size = 0;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    fisher = maxFisher(sap);
    for(j=0; j<sap->sizeC; j++) {
        res.table[j].fisher = fisher[j];
        res.table[j].corrected = 0.;
    }
    index = qsortIndirect(fisher, sap->sizeC, sizeof(double), compareDouble);
    for(j=0; j<sap->sizeC && fisher[index[j]]<threshold*0.000001; j++)
        removeClass(index[j], sap);
    fishTmp = getFisherClass(sap, j);
    effTmp = 0.;
    for(i=0; i<niter; i++) {
        double min;
        shuffle(r, sap->itemC, sap->sizeG);
        min = minFisher(sap);
        if(min<=fishTmp)
            effTmp++;
    }
    effTmp /= (double) niter;
printf("%dth atom %d\t%.2lE\t%.2lE\t%.2lE\n", j+1, (int) index[j], fisher[index[j]], fishTmp, effTmp);
    if(effTmp <= threshold) {
        res.table[res.size].corrected = effTmp;
        res.table[res.size].fisher = fisher[index[j]];
        res.table[res.size].atom = index[j];
        res.table[res.size].ontology = j;
        res.table[res.size].effC = j;
        res.table[res.size].effO = j;
        res.table[res.size].effCO = j;
        res.size++;
    }
    j++;
    for(; j<sap->sizeC && effTmp<threshold; j++) {
        removeClass(index[j-1], sap);
        fishTmp = getFisherClass(sap, j);
        effTmp = 0.;
        for(i=0; i<niter; i++) {
            double min;
            shuffle(r, sap->itemC, sap->sizeG);
            min = minFisher(sap);
            if(min<=fishTmp)
                effTmp++;
        }
        effTmp /= (double) niter;
printf("%dth atom %d\t%.2lE\t%.2lE\t%.2lE\n", j+1, (int) index[j], fisher[index[j]], fishTmp, effTmp);
        if(effTmp <= threshold) {
            res.table[res.size].corrected = effTmp;
            res.table[res.size].fisher = fisher[index[0]];
            res.table[res.size].atom = index[0];
            res.table[res.size].ontology = 0;
            res.table[res.size].effC = 0;
            res.table[res.size].effO = 0;
            res.table[res.size].effCO = 0;
            res.size++;
        }
    }
    free((void*)index);
    free((void*)fisher);
    gsl_rng_free (r);
    res.table = (TypeSignificant*) realloc(res.table, res.size*sizeof(TypeSignificant));
    return res;
}



/*TypeSignificantTable getBenjaminiHochberg(TypeStatAnnotationPartition *sap, double threshold, int niter) {
    int *eff, i, j, k, c, inc_buffer = sap->startG[sap->sizeG]/10, size_buffer = inc_buffer, ind;
    double *fisher, effTmp;
    const gsl_rng_type *T;
    TypeSignificantTable res;
    size_t *index;
    res.table = (TypeSignificant*) malloc(sap->sizeC*sizeof(TypeSignificant));
    res.size = 0;
    fisher = maxFisher(sap);
    for(j=0; j<sap->sizeC; j++) {
        res.table[j].fisher = fisher[j];
        res.table[j].corrected = 0.;
    }
    index = qsortIndirect(fisher, sap->sizeC, sizeof(double), compareDouble);
    effTmp = 0.;
    for(j=0; j<sap->sizeC; j++) {
        res.table[index[j]].corrected = ((double) j)*;
    }
    free((void*)index);
    free((void*)fisher);
    res.table = (TypeSignificant*) realloc(res.table, res.size*sizeof(TypeSignificant));
    return res;
}
*/

TypeSignificantTable getHochberg(TypeStatAnnotationPartition *sap) {
    int *yet, i, j, k, c, inc_buffer = sap->startG[sap->sizeG]/10, size_buffer = inc_buffer, ind;
    double *fisher, *corrected, effTmp, bonf;
    TypeSignificantTable res;
    size_t *index;
    bonf = (double) sap->sizeC*sap->sizeA+1.;
    res.table = (TypeSignificant*) malloc(sap->sizeC*sizeof(TypeSignificant));
    res.size = 0;
    bonf = ((double) sap->sizeC)*((double) sap->sizeA);
    fisher = allFisher(sap);
    index = qsortIndirect(fisher, sap->sizeC*sap->sizeA, sizeof(double), compareDouble);
    corrected = (double*) malloc(sap->sizeC*sap->sizeA*sizeof(double));
    for(j=0; j<sap->sizeC*sap->sizeA; j++)
        corrected[index[j]] = fisher[index[j]]*(bonf-((double)j));
    yet = (int*) malloc(sap->sizeC*sizeof(int));
    for(j=0; j<sap->sizeC; j++)
        yet[j] = 0;
    ind = 0;
    for(j=0; j<sap->sizeC*sap->sizeA; j++) {
        int atom = index[j]/sap->sizeA;
        if(yet[atom] == 0) {
            yet[atom] = 1;
            res.table[ind].atom = atom;
            res.table[ind].ontology = index[j]%sap->sizeA;
            res.table[ind].fisher = fisher[index[j]];
            res.table[ind].corrected = corrected[index[j]];
            res.table[ind].effC = 0;
            res.table[ind].effO = 0;
            res.table[ind].effCO = 0;
            ind++;
        }
    }
    free((void*)yet);
    free((void*)index);
    free((void*)fisher);
    free((void*)corrected);
    return res;
}

TypeSignificantTable getBonferroni(TypeStatAnnotationPartition *sap) {
    int *eff, i, j, k, c, inc_buffer = sap->startG[sap->sizeG]/10, size_buffer = inc_buffer, ind;
    double *fisher, effTmp, bonf;
    TypeSignificantTable res;
    size_t *index;
    res.table = (TypeSignificant*) malloc(sap->sizeC*sizeof(TypeSignificant));
    res.size = sap->sizeC;
    bonf = ((double) sap->sizeC)*((double) sap->sizeA);
    fisher = maxFisher(sap);
    index = qsortIndirect(fisher, sap->sizeC, sizeof(double), compareDouble);
    for(j=0; j<sap->sizeC; j++) {
        res.table[j].atom = index[j];
        res.table[j].ontology = 0;
        res.table[j].fisher = fisher[index[j]];
        res.table[j].corrected = fisher[index[j]]*bonf;
        res.table[j].effC = 0;
        res.table[j].effO = 0;
        res.table[j].effCO = 0;
    }
    free((void*)index);
    free((void*)fisher);
    return res;
}

TypeSignificantTable getSignificantTableEmpiricalBis(TypeStatAnnotationPartition *sap, double threshold, int niter) {
    int *eff, i, j, k, c, inc_buffer = sap->startG[sap->sizeG]/10, size_buffer = inc_buffer, ind;
    double *fisher, effTmp;
    const gsl_rng_type *T;
    gsl_rng *r;
    TypeSignificantTable res;
    size_t *index;
    res.table = (TypeSignificant*) malloc(sap->sizeC*sizeof(TypeSignificant));
    res.size = 0;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    fisher = maxFisher(sap);
    for(j=0; j<sap->sizeC; j++) {
        res.table[j].fisher = fisher[j];
        res.table[j].corrected = 0.;
    }
    index = qsortIndirect(fisher, sap->sizeC, sizeof(double), compareDouble);
    effTmp = 0.;
    for(i=0; i<niter; i++) {
        double min;
        shuffle(r, sap->itemC, sap->sizeG);
        min = minFisher(sap);
        if(min<=res.table[index[0]].fisher)
            effTmp++;
    }
    effTmp /= (double) niter;
printf("atom %d\t%.2lE\t%.2lE\n", (int) index[0], fisher[index[0]], effTmp);
    if(effTmp <= threshold) {
printf("atom %d\t%.2lE\t%.2lE\n", (int) index[0], fisher[index[0]], effTmp);
        res.table[res.size].corrected = effTmp;
        res.table[res.size].fisher = fisher[index[0]];
        res.table[res.size].atom = index[0];
        res.table[res.size].ontology = 0;
        res.table[res.size].effC = 0;
        res.table[res.size].effO = 0;
        res.table[res.size].effCO = 0;
        res.size++;
    }
    for(j=1; j<sap->sizeC && effTmp<threshold; j++) {
        removeClass(index[j-1], sap);
        for(i=0; i<niter; i++) {
            double min;
            shuffle(r, sap->itemC, sap->sizeG);
            min = minFisher(sap);
            if(min<=res.table[index[j]].fisher)
                effTmp++;
        }
        effTmp /= (double) niter;
printf("atom %d\t%.2lE\t%.2lE\n", (int) index[j], fisher[index[j]], effTmp);
        if(effTmp <= threshold) {
            res.table[res.size].corrected = effTmp;
            res.table[res.size].fisher = fisher[index[0]];
            res.table[res.size].atom = index[0];
            res.table[res.size].ontology = 0;
            res.table[res.size].effC = 0;
            res.table[res.size].effO = 0;
            res.table[res.size].effCO = 0;
            res.size++;
        }
    }
    free((void*)index);
    free((void*)fisher);
    gsl_rng_free (r);
    res.table = (TypeSignificant*) realloc(res.table, res.size*sizeof(TypeSignificant));
    return res;
}




void freeStatAnnotationPartition(TypeStatAnnotationPartition *sap) {
    int n;
    if(sap->nameG != NULL) {
        for(n=0; n<sap->sizeG; n++)
            if(sap->nameG[n] != NULL)
                free((void*)sap->nameG[n]);
        free((void*)sap->nameG);
    }
    if(sap->nameA != NULL) {
        for(n=0; n<sap->sizeA; n++)
            if(sap->nameA[n] != NULL)
                free((void*)sap->nameA[n]);
        free((void*)sap->nameA);
    }
    if(sap->mat != NULL) {
        for(n=0; n<sap->sizeG; n++)
            if(sap->mat[n] != NULL)
                free((void*)sap->mat[n]);
        free((void*)sap->mat);
    }
    if(sap->atom != NULL) {
        free((void*)sap->atom);
    }
    if(sap->startC != NULL)
        free((void*) sap->startC);
    if(sap->itemC != NULL)
        free((void*) sap->itemC);
    if(sap->startG != NULL)
        free((void*) sap->startG);
    if(sap->itemG != NULL)
        free((void*) sap->itemG);
    if(sap->startA != NULL)
        free((void*) sap->startA);
    if(sap->itemA != NULL)
        free((void*) sap->itemA);
    free((void*)sap);
}


void removeClass(int cl, TypeStatAnnotationPartition *sap) {
    int *index, *start, i, j, newSize, ind, c;

    index = (int*) malloc(sap->sizeG*sizeof(int));
    for(i=0; i<sap->sizeG; i++)
        index[i] = 0;
    for(i=sap->startC[cl]; i<sap->startC[cl+1]; i++)
        index[i] = -1;
    ind = 0;
    for(i=0; i<sap->sizeG; i++)
        if(index[i]>=0)
            index[i] = ind++;
    newSize = ind;
    start = (int*) malloc((sap->sizeC+1)*sizeof(int));
    ind = 0;
    for(c=0; c<sap->sizeC; c++) {
        start[c] = ind;
        for(i=sap->startC[c]; i<sap->startC[c+1]; i++)
            if(index[sap->itemC[i]]>=0)
                sap->itemC[ind++] = index[sap->itemC[i]];
    }
    start[sap->sizeC] = ind;
    for(c=0; c<=sap->sizeC; c++)
        sap->startC[c] = start[c];
    free((void*) start);
    start = (int*) malloc((newSize+1)*sizeof(int));
    ind = 0;
    for(i=0; i<sap->sizeG; i++)
        if(index[i]>=0) {
            start[index[i]] = ind;
            for(j=sap->startG[i]; j<sap->startG[i+1]; j++)
                sap->itemG[ind++] = sap->itemG[j];
        }
    start[newSize] = ind;
    for(i=0; i<=newSize; i++)
        sap->startG[i] = start[i];
    free((void*) start);
    sap->sizeG = newSize;
    start = (int*) malloc((sap->sizeA+1)*sizeof(int));
    ind = 0;
    for(i=0; i<sap->sizeA; i++) {
        start[i] = ind;
        for(j=sap->startA[i]; j<sap->startA[i+1]; j++)
            if(index[sap->itemA[j]]>=0)
                sap->itemA[ind++] = index[sap->itemA[j]];
    }
    start[sap->sizeA] = ind;
    for(i=0; i<=sap->sizeA; i++)
        sap->startA[i] = start[i];
    free((void*) start);
}


