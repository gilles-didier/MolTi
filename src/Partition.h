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




#ifndef PartitionF
#define PartitionF

#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum TCP {
    correctedRandIndex=1,
    standardJacquardIndex,
    standardRandIndex,
} TypePartitionCompareIndex;

typedef struct PARTITION {
    int sizeItem, sizeAtom, *atom;
} TypePartition;

typedef struct PARTITION_LIST {
    int sizeItem, sizeAtom, *start, *next;
} TypePartitionList;

typedef struct PARTITION_COMPACT {
    int sizeAtom, *start, *item;
} TypePartitionCompact;

/*compare part1 vs part2*/
double comparePart(TypePartitionCompareIndex typ, TypePartition *part1, TypePartition *part2);
void fprintPartition(FILE *f, TypePartition *part, char **name);
double logProbPart(TypePartition *part, double **distri);
TypePartitionList *getPartitionList(TypePartition *part);
TypePartitionCompact *getPartitionCompact(TypePartition *part);
void freePartitionList(TypePartitionList *pl);
void freePartitionCompact(TypePartitionCompact *pl);
void fprintPartitionList(FILE *f, TypePartition *part, char **name);
void fprintPartitionClustNSee(FILE *f, TypePartition *part, char **name);
void fprintPartitionStandard(FILE *f, TypePartition *part, char **name);
double comparePartDiff(TypePartitionCompareIndex type, TypePartition *part1, char **name1,  TypePartition *part2, char **name2);
TypePartition getPermutedPartition(TypePartition *part);
/*read partition in Clustnsee format*/
TypePartition readPartition(FILE *f, char***name);
int *getClassSize(TypePartition *part);
double entropyPartition(TypePartition *part);
double mutualInformationPartition(TypePartition *part1, TypePartition *part2);
double normalizedMutualInformationPartition(TypePartition *part1, TypePartition *part2);
double getPartitionDensity(TypePartition *part);
TypePartition readPartitionLine(FILE *f, char***name);
void keepIntersection(TypePartition *part1, char **name1, TypePartition *part2, char **name2, char ***name);
void freePartition(TypePartition *part);
/*read partition in Line Number format*/
TypePartition readPartitionNumber(FILE *f);
#ifdef __cplusplus
}
#endif

#endif
