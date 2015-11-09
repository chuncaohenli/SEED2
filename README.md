# SEEDg

[What is SEEDg] (#What is SEEDg)  
[How SEEDg works?] (#How SEEDg works)

[How is SEEDg different from SEED?] (#How is SEEDg different from SEED?)  
[How to use SEEDg?] (#How to use SEEDg?)  

<a name="What is SEEDg"/>
##What is SEEDg

SEEDg is the 2nd version of SEED(https://github.com/baoe/SEED), a software for clustering large sets of Next Generation Sequences (**only human genome is the effective input**) with hundreds of millions of reads in a time and memory efficient manner. 
<a name="How is SEEDg different from SEED?"/>

<a name="How SEEDg works">
##How SEEDg works
##1. A glance of SEEDg

![flowchart](http://1.easybuy1.sinaapp.com/SEEDg/flowchart.PNG)

##2. Algorithm detail
This part introduces the process of SEEDg in detail.
<a name="Step 1"/>
###Step 1-SEED on all sequences
Use original *SEED* algorithm (https://github.com/baoe/SEED) to cluster reads for the first time and get clusters result. The mismatch and shift are small so the rule to judge two clusters are simialr is strict 
![SEED](http://1.easybuy1.sinaapp.com/SEEDg/seedg1.PNG)
###Step 2-SEED on center sequences
1. Extract center sequence from each cluster, and use original *SEED* to cluster all center sequences. 
![extract center sequences](http://1.easybuy1.sinaapp.com/SEEDg/seedg2.PNG)
In this time, the mismatch and shift are enlarged and the rule to judge the similarity is loosen. We hope to optimize the result of first step. 
![SEED on cores](http://1.easybuy1.sinaapp.com/SEEDg/seedg3.PNG)
2. Adjust the cluster result
According to the cluster result of center sequences, we can adjust the cluster result of the [Step 1](# Step 1)

![kmeans](http://1.easybuy1.sinaapp.com/SEEDg/seedg4.PNG)

###Step 3-kmeans on large clusters
Now we successfully use high-efficient algorithm *SEED* to get a rough cluster result, and next step we will find out those clusters which are too large and use *kmeans* algorithm to split them into small sub-clusters precisely.
![kmeans](http://1.easybuy1.sinaapp.com/SEEDg/seedg5.PNG)
####kmeans - 2 strategies
In kmeans algorithm, I need to represent the read with array. And I have two strategies, each has pros and cons.
<a name="1map">
####1st mapping algorithm-KM1
Use a N-D array to represent a N-length read. It will cost much memory and the rule is strict so the result is precise.
eg.
For two reads
    ACGTACGTACGT
    AACGTACGTACG
they will be represented as (A-1,C-2,G-3,T-4)

    [1,2,3,4,1,2,3,4,1,2,3,4]
    [1,1,2,3,4,1,2,3,4,1,2,3]
calculate the cosine similarity,

    SIM = 0.84
<a name="2map">
####2nd mapping algorithm-KM2
Use the percentage of ACGT in a read, a 4-D array to represent a read. In this way, there are less memory needed, and the algorithm can tolerant more differents so the result is less precise but can cover more special situations,like

eg.
For two reads

    ACGTACGTACGT
    AACGTACGTACG
they will be represented as ([ratio of A,ratio of C, ratio of G, ratio of T])

    [0.25,0.25,0.25,0.25]
    [0.33,0.25,0.25,0.17]
calculate the cosine similarity,

    SIM = 0.97

Compared with KM1, KM2 are more likely to get the right result, to classify these two reads into one cluster.





##3. Experiment result
I made some experiments to test the performance of the SEED  and compare the result with SEED. ***Jaccard index*** is used to value the precision.

| Data set| SEED   |  SEEDg  |
| -----   | :---:  | :----:  |
| SRR1698795| 0.9859 |   0.9955     |
| SRR2035183|   0.993411   |   0.997966   |

##How is SEEDg different from SEED?

SEED is an efficient software for clustering large sets of hundreds of millions of reads, but sometimes the cluster result can be uneven, and the size of the cluster can be too large or too small. 

SEEDg is designed to solve this problem. And the differences are given below.

1. SEEDg needs to specify mismatch and shift for 2 times, compared with SEED
2. SEEDg offers more sensitivity tuning options for users
3. SEEDg gets more even cluster results while it cost more time than SEED

<a name="How to use SEEDg?"/>
##How to use SEEDg?

###Input&Output

`--input <FILE>`    input fastq file
Only FASTQ format is supported in the current version. The sequence length should be between 21 bp and 100 bp with the max variation of 
5 bp.

`--output <FILE>` output file

###Options

####Mode
`--default` run in default mode,use default spaced seed

`--fast`   	run in fast mode
use a bigger spaced seed weightto save running time. It is only applicable for sequences longer than 58 bp and may need more memory

`--short`  		run in short mode
use a smaller spaced seeds weight for sequences as short as 21 bp. This setting often results in longer compute times

####Sensitivity tuning

`--M1 `	max # mismatches in the first seed alignment (0 - 3, default 3);
mismatch is the maximum number of mismatches allowed from the center sequence in each cluster

`--S1 `	max # shifts in first seed alignment (0 - 6, default 3)
shift is the maximum number of shifts allowed from the center sequence in each cluster

`--M2 `	max # mismatches in second seed alignment (0 - 20, default 10)

`--S2 `	max # shifts in second seed alignment (0 - 6, default 3)

`--QV1`    	threshold for the base call quality values (QV);
QV1 is the threshold for the base call quality values (QV) that are provided in the FASTQ files as Phred scores. SEED ignores those mismatches where the sum of the Phred scores of the mismatching bases is lower than the specified QV1 threshold value (0 - 2 * 93). The default value for QV1 is 0

`--QV2`    	another QV threshold;
It prevents co-clustering of sequences where the sum of all mismatched positions is higher than the threshold value (0 - 40 * 93). The default value for QV2 is 40 * 93

`--reverse`	co-cluster sequences in sense and anti-sense orientation (reverse and complement)
	
`--KM1` use [1st mapping strategy](#1map) for K-means part (default choice)

`--KM2` use [2nd mapping strategy](#2map) for K-means part

####Others

`--input2 <FILE>`	specifies the paired sequences

###Examples

1. Basic usage:
```
$ SEEDg –-input input.fastq –-output output.txt
```

2. Run in user defined mode and tune the sensitivity manually
```
$ SEEDg –-input input.fastq -–output output.txt –L 30 –W 12 –S1 2 –M1 4 --reverse
```









