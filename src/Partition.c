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




#include <assert.h>
#include <string.h>
#include <math.h>
#include "Partition.h"
#include "Utils.h"
#define CLUSTERID "ClusterID"

void freePartition(TypePartition *part) {
    if(part != NULL) {
       if(part->atom != NULL)
           free((void*) part->atom);
       free((void*) part);
    }
}

double logProbPart(TypePartition *part, double **distri) {
    int i;
    double sum = 0.;

    for(i=0; i<part->sizeItem; i++)
        if(distri[i][part->atom[i]] > 0.)
            sum += log(distri[i][part->atom[i]]);
        else
            return 0.;
    return sum;
}

/*print partition*/
void fprintPartition(FILE *f, TypePartition *part, char **name) {
    if(name) {
        int i;
        for(i=0; i<part->sizeItem; i++)
            fprintf(f, "%s\t%d\n", name[i], part->atom[i]);
    } else {
        int i;
        for(i=0; i<part->sizeItem; i++)
            fprintf(f, "%d\t%d\n", i, part->atom[i]);
    }
}

/*print partition*/
void fprintPartitionList(FILE *f, TypePartition *part, char **name) {
    TypePartitionList *pl = getPartitionList(part);
    if(name) {
        int c,i;
        for(c=0; c<pl->sizeAtom; c++) {
            fprintf(f, ">>Class %d\n", c+1);
            for(i=pl->start[c]; i>=0; i=pl->next[i])
                fprintf(f, "%s ", name[i]);
            fprintf(f, "\n");
        }
    } else {
        int c,i;
        for(c=0; c<pl->sizeAtom; c++) {
            fprintf(f, ">>Class %d\n", c+1);
            for(i=pl->start[c]; i>=0; i=pl->next[i])
                fprintf(f, "%d ", i+1);
            fprintf(f, "\n");
        }
    }
    freePartitionList(pl);
}

int *getClassSize(TypePartition *part) {
    int *eff, i, c;

    eff = (int*) malloc(part->sizeAtom*sizeof(int));
    for(c=0; c<part->sizeAtom; c++)
        eff[c] = 0;
    for(i=0; i<part->sizeItem; i++)
        eff[part->atom[i]]++;
    return eff;
}

/*permute partition*/
TypePartition getPermutedPartition(TypePartition *part) {
    TypePartitionList *pl = getPartitionList(part);
    TypePartition res;
    int i, j, c, *elt, ind, *start;

    start = (int*) malloc((part->sizeAtom+1)*sizeof(int));
    elt = (int*) malloc(part->sizeItem*sizeof(int));
    ind = 0;
    for(c=0; c<pl->sizeAtom; c++) {
        start[c] = ind;
        for(i=pl->start[c]; i>=0; i=pl->next[i])
            elt[ind++] = i;
    }
    start[pl->sizeAtom] = part->sizeItem+1;
/*shuffle*/
    for(i=0; i<pl->sizeItem; i++) {
        int j = (rand() % (pl->sizeItem-i))+i, tmp;
        tmp = elt[i];
        elt[i] = elt[j];
        elt[j] = tmp;
    }
    res.sizeItem = part->sizeItem;
    res.sizeAtom = part->sizeAtom;
    res.atom = (int*) malloc(res.sizeItem*sizeof(int));
    for(c=0; c<pl->sizeAtom; c++)
        for(i=start[c]; i<start[c+1]; i++)
            res.atom[elt[i]] = c;
    return res;
}

/*print partition*/
void fprintPartitionClustNSee(FILE *f, TypePartition *part, char **name) {
    TypePartitionList *pl = getPartitionList(part);
    fprintf(f, "#ClustnSee analysis export\n");
    if(name) {
        int c,i;
        for(c=0; c<pl->sizeAtom; c++) {
            fprintf(f, "ClusterID:%d||\n", c+1);
            for(i=pl->start[c]; i>=0; i=pl->next[i])
                fprintf(f, "%s\n", name[i]);
            fprintf(f, "\n");
        }
    } else {
        int c,i;
        for(c=0; c<pl->sizeAtom; c++) {
            fprintf(f, "ClusterID:%d||\n", c+1);
            for(i=pl->start[c]; i>=0; i=pl->next[i])
                fprintf(f, "%d\n", i+1);
            fprintf(f, "\n");
        }
    }
    freePartitionList(pl);
}

/*print partition*/
void fprintPartitionStandard(FILE *f, TypePartition *part, char **name) {
    TypePartitionList *pl = getPartitionList(part);
    if(name) {
        int c,i;
        for(c=0; c<pl->sizeAtom; c++) {
            for(i=pl->start[c]; i>=0; i=pl->next[i])
                fprintf(f, "%s\t", name[i]);
            fprintf(f, "\n");
        }
    } else {
        int c,i;
        for(c=0; c<pl->sizeAtom; c++) {
            for(i=pl->start[c]; i>=0; i=pl->next[i])
                fprintf(f, "%d\t", i+1);
            fprintf(f, "\n");
        }
    }
    freePartitionList(pl);
}

double entropyPartition(TypePartition *part) {
    double *p, ent;
    int c, i;
    p = (double*) malloc(part->sizeAtom*sizeof(double));
    for(c=0; c<part->sizeAtom; c++)
        p[c] = 0.;
    for(i=0; i<part->sizeItem; i++)
        p[part->atom[i]]++;
    for(c=0; c<part->sizeAtom; c++)
        p[c] /= (double) part->sizeItem;
    ent = 0.;
    for(c=0; c<part->sizeAtom; c++)
        ent -= p[c]*log(p[c]);
    free((void*)p);
    return ent;
}
double mutualInformationPartition(TypePartition *part1, TypePartition *part2) {
    double **p12, *p1, *p2, mi;
    int c, d, i;

    if(part1->sizeItem != part2->sizeItem)
        return -1;
    p1 = (double*) malloc(part1->sizeAtom*sizeof(double));
    p2 = (double*) malloc(part2->sizeAtom*sizeof(double));
    p12 = (double**) malloc(part1->sizeAtom*sizeof(double*));
    for(c=0; c<part1->sizeAtom; c++) {
        p12[c] = (double*) malloc(part2->sizeAtom*sizeof(double));
        for(d=0; d<part2->sizeAtom; d++)
            p12[c][d] = 0.;
    }
    for(c=0; c<part1->sizeAtom; c++)
        p1[c] = 0.;
    for(d=0; d<part2->sizeAtom; d++)
        p2[d] = 0.;
    for(i=0; i<part1->sizeItem; i++) {
        p1[part1->atom[i]]++;
        p2[part2->atom[i]]++;
        p12[part1->atom[i]][part2->atom[i]]++;
    }
    for(c=0; c<part1->sizeAtom; c++)
        p1[c] /= (double) part1->sizeItem;
    for(d=0; d<part2->sizeAtom; d++)
        p2[d] /= (double) part2->sizeItem;
    for(c=0; c<part1->sizeAtom; c++)
        for(d=0; d<part2->sizeAtom; d++)
            p12[c][d] /= (double) part1->sizeItem;
    mi = 0.;
    for(c=0; c<part1->sizeAtom; c++)
        for(d=0; d<part2->sizeAtom; d++)
            mi += p12[c][d]*(log(p12[c][d])-log(p1[c])-log(p2[d]));
    for(c=0; c<part1->sizeAtom; c++)
        free((void*) p12[c]);
    free((void*) p12);
    free((void*) p1);
    free((void*) p2);
    return mi;
}

double normalizedMutualInformationPartition(TypePartition *part1, TypePartition *part2) {
    return 2.*mutualInformationPartition(part1, part2)/(entropyPartition(part1)+entropyPartition(part2));
}

/*compare part1 vs part2*/
double comparePart(TypePartitionCompareIndex type, TypePartition *part1, TypePartition *part2) {
    int i,ii, j,jj, k, ci, cj, common, appar1, appar2, clasmax=0, NbR, NbS, cas, kas, *card1, *card2, **inter;
    double Max, Exp;

    if(part1->sizeItem != part2->sizeItem)
        exitProg(ErrorExec, "Trying to compare two partitions not on the same set");
    switch(type) {
        case correctedRandIndex:
            //  Indice de inter corrigé
        default:
            clasmax = UTILS_MAX(part1->sizeAtom, part2->sizeAtom); // Nb. max de atomes dans le fichier des partitions
            card1 = (int*) malloc(clasmax*sizeof(int));		// Somme Ligne de inter
            card2 = (int*) malloc(clasmax*sizeof(int));		// Somme colonne
            inter = (int**) malloc(clasmax*sizeof(int *));
            assert(card1 != NULL && card2 != NULL && inter != NULL);
            for(i=0; i<clasmax; i++) {
                inter[i] = (int*) malloc(clasmax*sizeof(int));
                assert(inter[i] != NULL);
            }
            //  Initialisation du tableau de contingence des atomes
            for(i=0; i<clasmax; i++) {
                card1[i]=0; card2[i]=0; //cardinal de la atome i dans la partition 1, cardinal de la atome j dans la partition 2
                for(j=0; j<clasmax; j++)
                    inter[i][j]=0; //intersection de la atome i de la partition 1 et de la atome j de la partition 2
            }
            // Table des cardinaux des intersections des atomes
            for(k=0; k<part1->sizeItem; k++) {
                inter[part1->atom[k]][part2->atom[k]]++;
                card1[part1->atom[k]]++;
                card2[part2->atom[k]]++;
            }
            common=0; //nombre de paires communes entre les deux partitions
            for(i=0; i<part1->sizeAtom; i++)
                for(j=0; j<part2->sizeAtom; j++)
                    common += (inter[i][j]*(inter[i][j]-1))/2;
            appar1=0; //nombre de paires appariées dans la partition 1
            appar2=0; //nombre de paires appariées dans la partition 2
            for(i=0; i<part1->sizeAtom; i++)
                appar1 += (card1[i]*(card1[i]-1))/2;
            for(j=0; j<part2->sizeAtom; j++)
                appar2 += (card2[j]*(card2[j]-1))/2;
            free((void*)card1);
            free((void*)card2);
            for(i=0; i<clasmax; i++)
                free((void*)inter[i]);
            free((void*)inter);
//if((1.-((double)((part1->sizeItem*(part1->sizeItem-1))*common-2.*appar1*appar2))/((double)(((part1->sizeItem*(part1->sizeItem-1))*(appar1+appar2))/2.-2.*appar1*appar2)))>1.) {
//printf("size %d, common %d, appar1 %d, appar2 %d\n", part1->sizeItem, common, appar1, appar2);
//printf("%.3lf/%.3lf\n", ((double)(((double)part1->sizeItem*((double)part1->sizeItem-1.))*((double)common)-2.*((double)appar1)*((double)appar2))), ((double)(((((double)part1->sizeItem)*((double)part1->sizeItem-1.))*(((double)appar1)+((double)appar2)))/2.-2.*((double)appar1)*((double)appar2))));
//}
            return ((double)(((double)part1->sizeItem*((double)part1->sizeItem-1.))*((double)common)-2.*((double)appar1)*((double)appar2)))/((double)(((((double)part1->sizeItem)*((double)part1->sizeItem-1.))*(((double)appar1)+((double)appar2)))/2.-2.*((double)appar1)*((double)appar2)));
            //((double)((part1->sizeItem*(part1->sizeItem-1))*common-2.*appar1*appar2))/((double)(((part1->sizeItem*(part1->sizeItem-1))*(appar1+appar2))/2.-2.*appar1*appar2));
    case standardJacquardIndex: // indice de Jaccard
            NbR=0; NbS=0;
            for(i=1; i<part1->sizeItem; i++)
                for(j=0; j<i; j++) {
                    if(part1->atom[i]==part1->atom[j])
                        cas=1;
                    else
                        cas=0;
                    if(part2->atom[i]==part2->atom[j])
                        kas=1;
                    else
                        kas=0;
                    if(cas==1 && kas==1)
                        NbR++;
                    else
                        if(cas==1 || kas==1)
                            NbS++; // Reunies au moins une fois
                }
                return 1.*NbR/(NbS+NbR);
    case standardRandIndex: // Taux de paires simultanement reunies ou separees - Indice de inter
            NbR=0;
            for(i=1; i<part1->sizeItem; i++)
                for(j=0; j<i; j++) {
                    cas = (part1->atom[i]==part1->atom[j])?1:0;
                    kas = (part2->atom[i]==part2->atom[j])?1:0;
                    if(cas==kas)
                        NbR++;
                }
            return 2.*NbR/part1->sizeItem/(part1->sizeItem-1);
    }
}

double comparePartDiff(TypePartitionCompareIndex type, TypePartition *part1, char **name1,  TypePartition *part2, char **name2) {
    TypeLexiTree *dict;
    int i, j, *ok1, *ok2, size = 0, size1 = 0, size2 = 0, *cl1, *cl2;
    TypePartition  tmp1, tmp2;
    double res;

    ok1 = (int*) malloc(part1->sizeItem*sizeof(int));
    ok2 = (int*) malloc(part2->sizeItem*sizeof(int));
    for(i=0; i<part1->sizeItem; i++)
        ok1[i] = -1;
    for(i=0; i<part2->sizeItem; i++)
        ok2[i] = -1;
    dict = newLexiTree();
    for(i=0; i<part1->sizeItem; i++)
        if(addWordLexiTree(name1[i], dict) != i)
            printf("problem\n");
    for(i=0; i<part2->sizeItem; i++) {
        int ind;
        if((ind = findWordLexi(name2[i], dict)) >= 0) {
            ok1[ind] = size;
            ok2[i] = size;
            size++;
        }
    }
    cl1 = (int*) malloc(part1->sizeAtom*sizeof(int));
    cl2 = (int*) malloc(part2->sizeAtom*sizeof(int));
    for(i=0; i<part1->sizeAtom; i++)
        cl1[i] = -1;
    for(i=0; i<part2->sizeAtom; i++)
        cl2[i] = -1;
    for(i=0; i<part1->sizeItem; i++)
        if(ok1[i]>=0 && cl1[part1->atom[i]] < 0)
            cl1[part1->atom[i]] = size1++;
    for(i=0; i<part2->sizeItem; i++)
        if(ok2[i]>=0 && cl2[part2->atom[i]] < 0)
            cl2[part2->atom[i]] = size2++;
    tmp1.sizeItem = size;
    tmp1.sizeAtom = size1;
    tmp1.atom = (int*) malloc(tmp1.sizeItem*sizeof(int));
    for(i=0; i<part1->sizeItem; i++)
        if(ok1[i]>=0)
            tmp1.atom[ok1[i]] = cl1[part1->atom[i]];
    tmp2.sizeItem = size;
    tmp2.sizeAtom = size2;
    tmp2.atom = (int*) malloc(tmp2.sizeItem*sizeof(int));
    for(i=0; i<part2->sizeItem; i++)
        if(ok2[i]>=0)
            tmp2.atom[ok2[i]] = cl2[part2->atom[i]];
/*printf("tmp1\n");
for(i=0; i<part1->sizeItem; i++)
    printf("ok1[%d] = %d -> atom %d -> %d\n", i, ok1[i], part1->atom[i], cl1[part1->atom[i]]);
fprintPartition(stdout, &tmp1, NULL);
printf("\ntmp2\n");
for(i=0; i<part2->sizeItem; i++)
    printf("ok2[%d] = %d -> atom %d -> %d\n", i, ok2[i], part2->atom[i], cl2[part2->atom[i]]);
fprintPartition(stdout, &tmp2, NULL);
*/
    res = comparePart(type, &(tmp1), &(tmp2));
    freeLexiTree(dict);
    free((void*)ok1);
    free((void*)cl1);
    free((void*)tmp1.atom);
    free((void*)ok2);
    free((void*)cl2);
    free((void*)tmp2.atom);
    return res;
}


TypePartitionList *getPartitionList(TypePartition *part) {
    TypePartitionList *pl;
    int i, c;

    pl = (TypePartitionList*) malloc(sizeof(TypePartitionList));
    pl->sizeAtom = part->sizeAtom;
    pl->sizeItem = part->sizeItem;
    pl->start = (int*) malloc(pl->sizeAtom*sizeof(int));
    pl->next = (int*) malloc(pl->sizeItem*sizeof(int));
    for(c=0; c<pl->sizeAtom; c++)
        pl->start[c] = -1;
    for(i=0; i<part->sizeItem; i++) {
        pl->next[i] = pl->start[part->atom[i]];
        pl->start[part->atom[i]] = i;
    }
    return pl;
}

TypePartitionCompact *getPartitionCompact(TypePartition *part) {
    int a, ind=0;
    TypePartitionList *pl = getPartitionList(part);
    TypePartitionCompact *pc = (TypePartitionCompact *) malloc(sizeof(TypePartitionCompact));
    pc->sizeAtom = part->sizeAtom;
    pc->start = (int*) malloc((pc->sizeAtom+1)*sizeof(int));
    pc->item = (int*) malloc((part->sizeItem)*sizeof(int));
    ind = 0;
    for(a=0; a<pc->sizeAtom; a++) {
        int i;
        pc->start[a] = ind;
        for(i=pl->start[a]; i!=-1; i=pl->next[i])
            pc->item[ind++] = i;
    }
    pc->start[a] = ind;
    freePartitionList(pl);
    return pc;
}

void freePartitionCompact(TypePartitionCompact *pl) {
    free((void*) pl->start);
    free((void*) pl->item);
    free((void*) pl);
}
void freePartitionList(TypePartitionList *pl) {
    free((void*) pl->start);
    free((void*) pl->next);
    free((void*) pl);
}

#define INC_SIZE_ITEM_BUF 200
#define MAX_NAME_SIZE 200

/*read partition in Line Number format*/
TypePartition readPartitionNumber(FILE *f) {
    char c;
    TypePartition part;
    int sizeBufItem = INC_SIZE_ITEM_BUF, index;
    TypeLexiTree *dict;
    dict = newLexiTree();

    part.sizeAtom = 0;
    part.sizeItem = 0;
    part.atom = (int*) malloc(sizeBufItem*sizeof(int));
    for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
    while(c != EOF) {
        char tmp[MAX_NAME_SIZE+1];
        int i;
        while(c != EOF && !isline(c)) {
            if(c == '\'' || c == '"') {
                c = fgetc(f);
                for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
                    tmp[i] = c;
                    c = fgetc(f);
                }
                if(c == '\'' || c == '"')
                    c = fgetc(f);
                else
                    exitProg(ErrorReading, "Missing closing \" or '...");
            } else {
                for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issep(c); i++) {
                    tmp[i] = c;
                    c = fgetc(f);
                }
            }
            if(i == MAX_NAME_SIZE)
                exitProg(ErrorExec, "Name too much long...");
            if(i>0) {
                tmp[i++] = '\0';
                index = atoi(tmp);
				if(index>part.sizeItem)
                    part.sizeItem = index;
                if(index>=sizeBufItem) {
                    sizeBufItem += INC_SIZE_ITEM_BUF;
                    part.atom = (int*) realloc((void*)part.atom, sizeBufItem*sizeof(int));
                }
				for(; c != EOF && issep(c); c = fgetc(f));
				if(c == '\'' || c == '"') {
					c = fgetc(f);
					for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
						tmp[i] = c;
						c = fgetc(f);
					}
					if(c == '\'' || c == '"')
						c = fgetc(f);
					else
						exitProg(ErrorReading, "Missing closing \" or '...");
				} else {
					for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issepline(c); i++) {
						tmp[i] = c;
						c = fgetc(f);
					}
				}
				if(i == MAX_NAME_SIZE)
					exitProg(ErrorExec, "Name too much long...");
				if(i>0) {
					tmp[i++] = '\0';
					part.atom[index] = atoi(tmp);
					if(part.atom[index]>=part.sizeAtom)
						part.sizeAtom = part.atom[index]+1;
				}
            }
            for(; c != EOF && issep(c); c = fgetc(f));
        }
        for(; c != EOF && issepline(c); c = fgetc(f));
    }
    part.sizeItem++;
    part.atom = (int*) realloc((void*)part.atom, part.sizeItem*sizeof(int));
    return part;
}


/*read partition in Line format*/
TypePartition readPartitionLine(FILE *f, char***name) {
    char c;
    TypePartition part;
    int sizeBufItem = INC_SIZE_ITEM_BUF, index;
    TypeLexiTree *dict;
    dict = newLexiTree();

    part.sizeAtom = 0;
    part.sizeItem = 0;
    part.atom = (int*) malloc(sizeBufItem*sizeof(int));
    for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
    while(c != EOF) {
        char tmp[MAX_NAME_SIZE+1];
        int i;
        while(c != EOF && !isline(c)) {
            if(c == '\'' || c == '"') {
                c = fgetc(f);
                for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
                    tmp[i] = c;
                    c = fgetc(f);
                }
                if(c == '\'' || c == '"')
                    c = fgetc(f);
                else
                    exitProg(ErrorReading, "Missing closing \" or '...");
            } else {
                for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issepline(c); i++) {
                    tmp[i] = c;
                    c = fgetc(f);
                }
            }
            if(i == MAX_NAME_SIZE)
                exitProg(ErrorExec, "Name too much long...");
            if(i>0) {
                tmp[i++] = '\0';
//printf("%s\n", tmp);
                if(findWordLexi(tmp, dict) >= 0) {
                    printf("Duplicate identifier %s(%d)\n", tmp, findWordLexi(tmp, dict));
                    fprintLexiTree(stdout, dict);
                    exit(1);
                }
                index = addWordLexiTree(tmp, dict);
                if(index>part.sizeItem)
                    part.sizeItem = index;
                if(index>=sizeBufItem) {
                    sizeBufItem += INC_SIZE_ITEM_BUF;
                    part.atom = (int*) realloc((void*)part.atom, sizeBufItem*sizeof(int));
                }
                part.atom[index] = part.sizeAtom;
            }
            for(; c != EOF && issep(c); c = fgetc(f));
        }
        part.sizeAtom++;
        for(; c != EOF && issepline(c); c = fgetc(f));
    }
    part.sizeItem++;
    part.atom = (int*) realloc((void*)part.atom, part.sizeItem*sizeof(int));
    *name = (char**) malloc(part.sizeItem*sizeof(char*));
    fillLexiTree(*name, dict);
    freeLexiTree(dict);
    return part;
}

/*read partition in Clustnsee format*/
TypePartition readPartition(FILE *f, char***name) {
    char c;
    TypePartition part;
    int sizeBufItem = INC_SIZE_ITEM_BUF, index;
    TypeLexiTree *dict;
    dict = newLexiTree();

    part.sizeAtom = -1;
    part.sizeItem = 0;
    part.atom = (int*) malloc(sizeBufItem*sizeof(int));
    for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
    while(c != EOF) {
        char tmp[MAX_NAME_SIZE+1];
        int i;
        if(c == '\'' || c == '"') {
            c = fgetc(f);
            for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
                tmp[i] = c;
                c = fgetc(f);
            }
            if(c == '\'' || c == '"')
                c = fgetc(f);
            else
                exitProg(ErrorReading, "Missing closing \" or '...");
        } else {
            for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issepline(c); i++) {
                tmp[i] = c;
                c = fgetc(f);
            }
        }
/*		if(i == 0)
            exitProg(ErrorExec, "Empty name...");
*/		if(i == MAX_NAME_SIZE)
            exitProg(ErrorExec, "Name too much long...");
        if(i>0) {
            tmp[i++] = '\0';
//printf("%s\n", tmp);
            if(tmp[0] == '#') {
                for(; c != EOF && !isline(c); c = fgetc(f));
            } else {
                if(strncmp(tmp, CLUSTERID, strlen(CLUSTERID)) == 0) {
//	printf("Class\n");
                    part.sizeAtom++;
                    for(; c != EOF && !isline(c); c = fgetc(f));
                } else {
//	printf("Ident\n");
                    if(findWordLexi(tmp, dict) >= 0) {
                        printf("Duplicate identifier %s(%d)\n", tmp, findWordLexi(tmp, dict));
                        fprintLexiTree(stdout, dict);
                        exit(1);
                    }
                    index = addWordLexiTree(tmp, dict);
                    if(index>part.sizeItem)
                        part.sizeItem = index;
                    if(index>=sizeBufItem) {
                        sizeBufItem += INC_SIZE_ITEM_BUF;
                        part.atom = (int*) realloc((void*)part.atom, sizeBufItem*sizeof(int));
                    }
                    part.atom[index] = part.sizeAtom;
                }
            }
        }
        for(; c != EOF && issepline(c); c = fgetc(f));
    }
    part.sizeItem++;
    part.sizeAtom++;
    part.atom = (int*) realloc((void*)part.atom, part.sizeItem*sizeof(int));
    *name = (char**) malloc(part.sizeItem*sizeof(char*));
    fillLexiTree(*name, dict);
    freeLexiTree(dict);
    return part;
}



double getPartitionDensity(TypePartition *part) {
    int i, j, c;
    double dens = 0.;
    for(i=0; i<part->sizeItem; i++)
        for(j=0; j<i; j++)
            if(part->atom[i] == part->atom[j])
                dens++;
    dens /= (double) (part->sizeItem*(part->sizeItem-1))/2;
    return dens;
}


void keepIntersection(TypePartition *part1, char **name1, TypePartition *part2, char **name2, char ***name) {
    TypeLexiTree *dict;
    int *indexI1, *indexI2, indI, *indexC1, *indexC2, *atom, indC, i, n;

    indexI1 = (int*) malloc(part1->sizeItem*sizeof(int));
    for(i=0; i<part1->sizeItem; i++)
        indexI1[i] = -1;
    indexC1 = (int*) malloc(part1->sizeAtom*sizeof(int));
    for(i=0; i<part1->sizeAtom; i++)
        indexC1[i] = -1;
    indexI2 = (int*) malloc(part2->sizeItem*sizeof(int));
    for(i=0; i<part2->sizeItem; i++)
        indexI2[i] = -1;
    indexC2 = (int*) malloc(part2->sizeAtom*sizeof(int));
    for(i=0; i<part2->sizeAtom; i++)
        indexC2[i] = -1;
    dict = getDictFromTable(name1, part1->sizeItem);
    indI = 0;
    for(i=0; i<part2->sizeItem; i++)
        if((n = findWordLexi(name2[i], dict)) >= 0) {
            indexI1[n] = indI;
            indexI2[i] = indI;
            indI++;
        }
    if(indI == 0) {
        fprintf(stderr, "Empty intersection between two partitions. Exit...\n");
        exit(1);
    }
    *name = (char**) malloc(indI*sizeof(char*));
    for(i=0; i<part1->sizeItem; i++)
        if(indexI1[i]>=0)
            (*name)[indexI1[i]] = strdup(name1[i]);
    atom = part1->atom;
    part1->atom = (int*) malloc(indI*sizeof(int));
    indC = 0;
    for(i=0; i<part1->sizeItem; i++) {
        if(indexI1[i] >= 0) {
    printf("i %d index %d/%d atom %d/%d -> %d/%d\n", i, indexI1[i], indI, atom[i], part1->sizeAtom, indexC1[atom[i]], indC);
            if(indexC1[atom[i]]<0)
                indexC1[atom[i]] = indC++;
    printf("i %d index %d/%d atom %d/%d -> %d/%d\n%s -< %d\n", i, indexI1[i], indI, atom[i], part1->sizeAtom, indexC1[atom[i]], indC, (*name)[indexI1[i]], indexC1[atom[i]]+1);
            part1->atom[indexI1[i]] = indexC1[atom[i]];
        }
    }
    part1->sizeAtom = indC;
    part1->sizeItem = indI;
    free((void*)atom);
    atom = part2->atom;
    part2->atom = (int*) malloc(indI*sizeof(int));
    indC = 0;
    for(i=0; i<part2->sizeItem; i++) {
        if(indexI2[i] >= 0) {
            if(indexC2[atom[i]]<0)
                indexC2[atom[i]] = indC++;
            part2->atom[indexI2[i]] = indexC2[atom[i]];
        }
    }
    part2->sizeAtom = indC;
    part2->sizeItem = indI;
    free((void*)atom);
    freeLexiTree(dict);
    free((void*)indexI1);
    free((void*)indexC1);
    free((void*)indexI2);
    free((void*)indexC2);
}


