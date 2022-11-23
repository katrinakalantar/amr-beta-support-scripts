# amr-beta-support-scripts
Scripts to support analysis of CZ ID AMR beta module. Will add these ad-hoc as needs arise for supporting various projects.


Script #1 (Nov 22) -- 

Quick pass at generating gene x sample matrices with filtering capability (to support GC trainings):

```
python3 compile_amr.py [input file] [min_read] [min_read_breadth] [min_read_depth]
```

* `input_file` is the `combined_amr_result.csv` file generated via bulk downloads modal in CZ ID
* `min_read` is the minimum number of reads
* `min_read_breadth` is the minimum read coverage breadth (values from 0 to 100, integers)
* `min_read_depth` is the minimum read coverage depth (values from 0 to 100, integers)


An example:
```
python3 ./scripts/compile_amr.py ./data/combined_amr_results.csv 5 10 1 my_output
```

This will generate two outputs -- 
* a `gene` x `samples` matrix, where counts are 0 if the gene is not present in the sample, 1 if it is present but doesn't pass the thresholds, and 2 if it is present and passes thresholds.
* a `drug class` x `samples` matrix, where counts correspond to the number of threshold-passing genes of that drug class