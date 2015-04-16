
CC = gcc

LDFLAGS = -g -L/usr/local/lib -L/usr/lib/ -lm 


ALL = molti-console test bonf

all: molti-console test bonf

molti-console: Compute.c
	$(CC) EdgesComposition.c Louvain.c Partition.c Graph.c  Utils.c Compute.c -o $@ -L/usr/local/lib -L/usr/lib/ -lm -g 

test: Test.c
	$(CC) Louvain.c Partition.c Graph.c  Utils.c ERMG.c RandomGraph.c Test.c -o $@ -L/usr/local/lib -L/usr/lib/ -lm -g

bonf: Bonferroni.c
	$(CC) Partition.c StatAnnotationPartition.c Annotation.c  Utils.c Bonferroni.c -o $@ -L/usr/local/lib -L/usr/lib/ -I/usr/local/include -lm -lgsl -lgslcblas -g 
