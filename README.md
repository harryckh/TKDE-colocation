
Notes
=======================

1. This code was used for the empirical study of the paper
"Fraction-Score: A Generalized Support Measure for Weighted and Maximal Co-location Pattern Mining".
2. This code is developed by Harry Kai-Ho Chan (email: kai-ho@ruc.dk).
3. This code is written in C/C++.
4. This code runs under Unix/Linux.

Usage
=======================

Step 1: Compile the source code
	
	make

Step 2: Run the code

	./a.out

Step 3: Input Files (Format of those files can be found in Appendix A)
By default, we have file "config.txt" for configuration

Step 4: Collect the co-location patterns and running statistics 
[you can ignore this step if you don't want to collect the information of
co-location patterns and running statistics]

the co-location patterns are stored in "result.txt"
which format is explained in Appendix D.

the running statistics are stored in "stat.txt"
which format is explained in Appendix E.


Appendix A. The format of config.txt
============================

	<Algorithm indicator> 
	<Measurement indicator>
	<dimension>
	<# of objects>
	<Location file>
	<# of keywords>
	<Keyword file>
	<IR-tree option>
	<IR-tree file>
	<# of querying times>
	<min_sup>
	<dist_thr>

Explanation of the content in config.txt

	<Algorithm indicator> 
		= 41: Filtering and Verification (Combinatorial-Search)
		= 42: Filtering and Verification (Dia-CoSKQ-Adapt)
		= 43: Filtering and Verification (Optimization-Search)
		= 46: Maximal Pattern Mining
		= 48: SGCT
		= 5: Join-less

	<Measurement indicator>
		= 1: Participation-based
		= 2: Fraction-Score

	<dimension>
		: the number of dimensions of the dataset.

	<# of objects>
		: the number of spatial objects.

	<Location file>
		: the file containing the locations of the spatial objects,
	which format is explained in Appendix II.

	<# of keywords>
		: the total number of all possible keywords contained by the objects.

	<Keyword file>
		: the file containing the keywords of the spatial objects,
	which format is explained in Appendix III.

	<IR-tree option>
		= 0: the IR-tree has been built (which will be built and stored in <IR-tree file>)
		= 1: the IR-tree has been built (which is stored in <IR-tree file file>)

	<IR-tree file>
		: the file for storing a new (or an existing) IR-tree.

	<# of querying times>
		: the number of queries that will be performed for the same setting, 
		the average statistics (i.e., query time) based on which will be used.

	<min_sup>
		: the minimum support threshold for the co-location patterns

	<dist_thr>
		: the distance threshold for defining a neighbor set

	(See file config.txt for example)

Appendix B. The format of "Location file"
============================

	<object ID1>, <1st coordinate>, <2nd coordinate>, ..., <d-th coordinate>
	<object ID2>, <1st coordinate>, <2nd coordinate>, ..., <d-th coordinate>
	<object ID3>, <1st coordinate>, <2nd coordinate>, ..., <d-th coordinate>
	...
	<object IDn>, <1st coordinate>, <2nd coordinate>, ..., <d-th coordinate>
	Note that n = # of objects, d = # of dimensions

(See file in the dataset folder for example)

Appendix C. The format of "Keyword file"
=============================
	<object ID1>, <keyword>, <weight>
	<object ID2>, <keyword>, <weight>
	<object ID3>, <keyword>, <weight>
	...
	<object IDn>, <keyword>, <weight>

(See file in the dataset folder for example)

Appendix D. The format of "result.txt"
=============================

	<# of size 1 patterns>
	<1st feature ID of the 1st size 1 pattern> <support>
	<1st feature ID of the 2nd size 1 pattern> <support>
	...
	<1st feature ID of the k-th size 1 pattern> <support>

	<# of size 2 patterns>
	<1st feature ID of the 1st size 2 pattern> <2nd feature ID of the 1st size 2 pattern> <support>
	<1st feature ID of the 2nd size 2 pattern> <2nd feature ID of the 2nd size 2 pattern> <support>
	...
	<1st feature ID of the k-th size 2 pattern> <2nd feature ID of the k-th size 2 pattern> <support>

	...
	<# of size m patterns>
	<1st feature ID of the 1st size m pattern> <2nd feature ID of the 1st size m pattern> ... <m-th feature ID of the 1st size m pattern><support>
	<1st feature ID of the 2nd size m pattern> <2nd feature ID of the 2nd size m pattern> ... <m-th feature ID of the 2nd size m pattern><support>
	...
	<1st feature ID of the k-th size m pattern> <2nd feature ID of the k-th size m pattern> ... <m-th feature ID of the k-th size m pattern><support>

(See file result.txt for example)

Appendix E. The format of "stat.txt"
=============================

	<the time of building the IR-tree>
	<the average time of performing a query>

	<Memory usage>
	<IR-tree memory usage>

	<Filter 1 count>
	<Filter 2 count>
	<Filter 3 count>
	<Filter 4 count>
	<Verification count>

	<Filter 1 time>
	<Filter 2 time>
	<Filter 3 time>
	<Filter 4 time>
	<Verification time>

(See file stat.txt for example)



