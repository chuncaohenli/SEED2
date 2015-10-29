# SEED2

##What is SEED2

SEED2 is the 2nd version of SEED(https://github.com/baoe/SEED), a software for clustering large sets of Next Generation Sequences (NGS) with hundreds of millions of reads in a time and memory efficient manner. 

##How is SEED2 different from SEED?

SEED is an efficient software for clustering large sets of hundreds of millions of reads, but sometimes the cluster result can be uneven, and the size of the cluster can be too large or too small. 

SEED2 is designed to solve this problem. And the differences are given below.

1. SEED2 needs to specify mismatch and shift for 2 times, compared with SEED
2. SEED2 gets more even cluster results while it cost more time than SEED

##How to use SEED2?

###Input&Output

--input <FILE>	input fastq file
Only FASTQ format is supported in the current version. The sequence length should be between 21 bp and 100 bp with the max variation of 
5 bp.

--output <FILE>	output file

###Options

####Mode
	--default		run in default mode
				use default spaced seed

	--fast   		run in fast mode
use a bigger spaced seed weight to save running time. It is only applicable for sequences longer than 58 bp and may need more memory

	--short  		run in short mode
use a smaller spaced seeds weight for sequences as short as 21 bp. This setting often results in longer compute times

####Sensitivity tuning

	-M1 <int>	max # mismatches in the first seed alignment (0 - 3, default 3);
mismatch is the maximum number of mismatches allowed from the center sequence in each cluster

	-S1 <int>	max # shifts in first seed alignment (0 - 6, default 3)
shift is the maximum number of shifts allowed from the center sequence in each cluster

	-M2 <int>	max # mismatches in first seed alignment (0 - 20, default 10)

	-S2 <int>	max # shifts in first seed alignment (0 - 6, default 3)

	--QV1    	threshold for the base call quality values (QV);
QV1 is the threshold for the base call quality values (QV) that are provided in the FASTQ files as Phred scores. SEED ignores those mismatches where the sum of the Phred scores of the mismatching bases is lower than the specified QV1 threshold value (0 - 2 * 93). The default value for QV1 is 0

	--QV2    	another QV threshold;
It prevents co-clustering of sequences where the sum of all mismatched positions is higher than the threshold value (0 - 6 * 93). The default value for QV2 is 6 * 93

	--reverse	co-cluster sequences in sense and anti-sense orientation (reverse and complement)
	
    --KM1 use [1st mapping strategy][1] for K-means part (default choice)
	--KM2 use [2nd mapping strategy][2] for K-means part

####Others

	--input2 <FILE>	specifies the paired sequences

###Examples

1. Basic usage:

        SEED2 –-input input.fastq –-output output.txt

2. Run in user defined mode and tune the sensitivity manually

        SEED2 –-input input.fastq -–output output.txt –L 30 –W 12 –S1 2 –M1 4 --reverse
