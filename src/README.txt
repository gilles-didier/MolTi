The MolTi software suite allows identifying communities from multiplex networks, and annotated the obtained clusters. MolTi performs both clustering and annotation enrichment tests from multiplex networks, and allow an easy exploration of the results.

--------------------

 'molti-console' takes as input several graph files in the simple "ncol" format (one interaction per line, the 2 interactors being separated by a tab), from which it detect communities by maximizing their Multiplex-modularity. It returns a file containing the community partition. Option '-p' sets the resolution parameter and Option '-s' provides the community partitions of all the individual graphs and that of their sum and union graphs (for comparison purpose).

usage: molti-console [options] <input file1> <input file2> ...

options are
	-o <file name>      set the output file name prefix
	-p <real number>    set the Newman modularity resolution parameter (set to 1 as default)
	-s                  compute partition on each graph individually and on the sum graph
	-h                  display help


--------------------

 'bonf' computes the annotation enrichment with an exact Fisher test and Bonferroni correction for multiple-testing, for all the communities of a given partition. It takes as input an annotation "ontology" file and a partition file (partition file may be in Clust'N'See \citep{Spinelli2013 format/option '-f c' or in flat format, i.e. a line per community/option '-f f'). It returns a file listing all the communities  significantly (i.e. with a corrected exact Fisher test lower than the value provided with option 'threshold') enriched in at least one ontology term (all the ontology terms significantly represented in a class are displayed below this one). The option '-o' allows providing a file containing descriptions of the ontology terms, which will be displayed in the output file.

usage: 'bonf' [options] <annotation file> <partition file> [<output file>]

options are
	-o <ontology descriptions file>     load ontology terms descriptions
	-t <real number>                    set the threshold (set to 0.001 as default)
	-f c or f                           indicate the partition format
	-h                                  display help

--------------------

 'test' first simulates random multiplex networks with g vertices and from 1 to t layers with a balanced community structure of c communities (i simulations for each number of layers). It next detects communities by using aggregation and multiplex-modularity approaches and computes the adjusted Rand index (and the normalized mutual information) between the reference community structure and the detected one. It writes a '.csv' file containing the means and the standard deviations of the adjusted Rand indexes for each method, with the following format:
 Column 1: number of layers
 Column 2 and 3: mean and standard deviation obtained with the multiplex-modularity approach
 Column 4 and 5: mean and standard deviation obtained with the sum-aggregation approach
 Column 6 and 7: mean and standard deviation obtained with the intersection
 Column 6 and 7: mean and standard deviation obtained with the union aggregation approach. 

usage: 'test' <options> <output file name>

Options:
	-g <number>	set the number of vertices of the random graphs (50 as default)
	-t <number>	set the max number of random graphs (10 as default)
	-c <number>	set the number of classes (3 as default)
	-i <number>	set the number of iterations (4 as default)
	-p <prob intra> <prob inter>	add a new pair of probas
	-a  <number>	set the modularity parameter (1 as default)
	-h	display help message

--------------------

'testth' is a multi-threaded version of 'test'. It has an extra-option:
	-m <number>	set the number of simultaneous threads

--------------------

'testgl' is a version of 'test' including a call to GenLouvain. It has an extra-option:
	-l <path>	set the path to GenLouvain matlab script

--------------------

'testglth' is a multi-threaded version of 'testgl' with the two extra-options above.
