
CC = gcc

LDFLAGS = -g -L/usr/local/lib -L/usr/lib/ -lm 


ALL = molti-console test testgl bonf consensus testth testglth

all: molti-console test testgl bonf consensus testth testglth

molti-console: Compute.c
	$(CC) EdgesComposition.c Louvain.c Partition.c Graph.c  Utils.c Compute.c -o $@ -L/usr/local/lib -L/usr/lib/ -lm -O3

consensus: Consensus.c
	$(CC) EdgesComposition.c Louvain.c Partition.c Graph.c  Utils.c Consensus.c -o $@ -L/usr/local/lib -L/usr/lib/ -lm -O3 

test: Test.c
	$(CC) Louvain.c Partition.c Graph.c  Utils.c ERMG.c RandomGraph.c Test.c -o $@ -L/usr/local/lib -L/usr/lib/ -lm -O3

testgl: TestGenLouvain.c
	$(CC) Louvain.c Partition.c Graph.c  Utils.c ERMG.c RandomGraph.c TestGenLouvain.c -o $@ -L/usr/local/lib -L/usr/lib/ -lm -O3
	
testth: TestThread.c
	$(CC)  -D_REENTRANT Louvain.c Partition.c Graph.c  Utils.c ERMG.c RandomGraph.c TestThread.c -o $@  -L/usr/local/lib -L/usr/lib/ -lm -lpthread -O3

testglth: TestThreadGenLouvain.c
	$(CC)  -D_REENTRANT Louvain.c Partition.c Graph.c  Utils.c ERMG.c RandomGraph.c TestThreadGenLouvain.c -o $@  -L/usr/local/lib -L/usr/lib/ -lm -lpthread -O3 

bonf: Bonferroni.c
	$(CC) Partition.c StatAnnotationPartition.c Annotation.c  Utils.c Bonferroni.c -o $@ -L/usr/local/lib -L/usr/lib/ -I/usr/local/include -lm -lgsl -lgslcblas -O3 
