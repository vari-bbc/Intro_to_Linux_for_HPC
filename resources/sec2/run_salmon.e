Version Info: This is the most recent version of salmon.
### salmon (selective-alignment-based) v1.10.0
### [ program ] => salmon 
### [ command ] => quant 
### [ threads ] => { 4 }
### [ libType ] => { A }
### [ index ] => { /varidata/research/projects/bbc/versioned_references/2022-10-06_14.25.40_v10/data/hg38_gencode/indexes/salmon/hg38_gencode/ }
### [ mates1 ] => { fastqs/SRR1039520_1.fastq.gz }
### [ mates2 ] => { fastqs/SRR1039520_2.fastq.gz }
### [ output ] => { salmon/SRR1039520 }
### [ validateMappings ] => { }
Logs will be written to salmon/SRR1039520/logs
[2025-10-03 10:28:40.892] [jointLog] [info] setting maxHashResizeThreads to 4
[2025-10-03 10:28:40.892] [jointLog] [info] Fragment incompatibility prior below threshold.  Incompatible fragments will be ignored.
[2025-10-03 10:28:40.892] [jointLog] [info] Usage of --validateMappings implies use of minScoreFraction. Since not explicitly specified, it is being set to 0.65
[2025-10-03 10:28:40.892] [jointLog] [info] Setting consensusSlack to selective-alignment default of 0.35.
[2025-10-03 10:28:40.892] [jointLog] [info] parsing read library format
[2025-10-03 10:28:40.892] [jointLog] [info] There is 1 library.
[2025-10-03 10:28:40.966] [jointLog] [info] Loading pufferfish index
[2025-10-03 10:28:40.982] [jointLog] [info] Loading dense pufferfish index.
-----------------------------------------
| Loading contig table | Time = 10.913 s
-----------------------------------------
size = 37226611
-----------------------------------------
| Loading contig offsets | Time = 271.43 ms
-----------------------------------------
-----------------------------------------
| Loading reference lengths | Time = 22.908 ms
-----------------------------------------
-----------------------------------------
| Loading mphf table | Time = 872.89 ms
-----------------------------------------
size = 3779163173
Number of ones: 37226610
Number of ones per inventory item: 512
Inventory entries filled: 72709
-----------------------------------------
| Loading contig boundaries | Time = 5.8189 s
-----------------------------------------
size = 3779163173
-----------------------------------------
| Loading sequence | Time = 690.69 ms
-----------------------------------------
size = 2662364873
-----------------------------------------
| Loading positions | Time = 5.2853 s
-----------------------------------------
size = 3457937385
-----------------------------------------
| Loading reference sequence | Time = 741.8 ms
-----------------------------------------
-----------------------------------------
| Loading reference accumulative lengths | Time = 26.398 ms
-----------------------------------------
[2025-10-03 10:29:05.677] [jointLog] [info] done
[2025-10-03 10:29:05.698] [jointLog] [info] Index contained 227257 targets
[2025-10-03 10:29:05.752] [jointLog] [info] Number of decoys : 194
[2025-10-03 10:29:05.752] [jointLog] [info] First decoy index : 227024 




[2025-10-03 10:29:06.021] [jointLog] [info] Automatically detected most likely library type as IU

[A[32mprocessed[31m 500000 [32mfragments[0m
hits: 2028752, hits per frag:  4.19191







[2025-10-03 10:29:09.156] [jointLog] [info] Computed 88371 rich equivalence classes for further processing
[2025-10-03 10:29:09.156] [jointLog] [info] Counted 465266 total reads in the equivalence classes 
[2025-10-03 10:29:09.158] [jointLog] [info] Number of mappings discarded because of alignment score : 149760
[2025-10-03 10:29:09.158] [jointLog] [info] Number of fragments entirely discarded because of alignment score : 12331
[2025-10-03 10:29:09.158] [jointLog] [info] Number of fragments discarded because they are best-mapped to decoys : 7197
[2025-10-03 10:29:09.158] [jointLog] [info] Number of fragments discarded because they have only dovetail (discordant) mappings to valid targets : 3713
[2025-10-03 10:29:09.192] [jointLog] [warning] Only 465266 fragments were mapped, but the number of burn-in fragments was set to 5000000.
The effective lengths have been computed using the observed mappings.

[2025-10-03 10:29:09.192] [jointLog] [info] Mapping rate = 93.0532%

[2025-10-03 10:29:09.192] [jointLog] [info] finished quantifyLibrary()
[2025-10-03 10:29:09.193] [jointLog] [info] Starting optimizer
[2025-10-03 10:29:09.311] [jointLog] [info] Marked 0 weighted equivalence classes as degenerate
[2025-10-03 10:29:09.337] [jointLog] [info] iteration = 0 | max rel diff. = 198.919
[2025-10-03 10:29:11.967] [jointLog] [info] iteration = 100 | max rel diff. = 18.9177
[2025-10-03 10:29:14.755] [jointLog] [info] iteration = 200 | max rel diff. = 0.559677
[2025-10-03 10:29:17.994] [jointLog] [info] iteration = 300 | max rel diff. = 2.38007
[2025-10-03 10:29:21.381] [jointLog] [info] iteration = 400 | max rel diff. = 0.152702
[2025-10-03 10:29:24.023] [jointLog] [info] iteration = 500 | max rel diff. = 0.0229273
[2025-10-03 10:29:25.327] [jointLog] [info] iteration = 551 | max rel diff. = 0.00931626
[2025-10-03 10:29:25.369] [jointLog] [info] Finished optimizer
[2025-10-03 10:29:25.369] [jointLog] [info] writing output 

Version Info: This is the most recent version of salmon.
### salmon (selective-alignment-based) v1.10.0
### [ program ] => salmon 
### [ command ] => quant 
### [ threads ] => { 4 }
### [ libType ] => { A }
### [ index ] => { /varidata/research/projects/bbc/versioned_references/2022-10-06_14.25.40_v10/data/hg38_gencode/indexes/salmon/hg38_gencode/ }
### [ mates1 ] => { fastqs/SRR1039521_1.fastq.gz }
### [ mates2 ] => { fastqs/SRR1039521_2.fastq.gz }
### [ output ] => { salmon/SRR1039521 }
### [ validateMappings ] => { }
Logs will be written to salmon/SRR1039521/logs
[2025-10-03 10:29:27.902] [jointLog] [info] setting maxHashResizeThreads to 4
[2025-10-03 10:29:27.902] [jointLog] [info] Fragment incompatibility prior below threshold.  Incompatible fragments will be ignored.
[2025-10-03 10:29:27.902] [jointLog] [info] Usage of --validateMappings implies use of minScoreFraction. Since not explicitly specified, it is being set to 0.65
[2025-10-03 10:29:27.902] [jointLog] [info] Setting consensusSlack to selective-alignment default of 0.35.
[2025-10-03 10:29:27.902] [jointLog] [info] parsing read library format
[2025-10-03 10:29:27.902] [jointLog] [info] There is 1 library.
[2025-10-03 10:29:27.902] [jointLog] [info] Loading pufferfish index
[2025-10-03 10:29:27.902] [jointLog] [info] Loading dense pufferfish index.
-----------------------------------------
| Loading contig table | Time = 10.625 s
-----------------------------------------
size = 37226611
-----------------------------------------
| Loading contig offsets | Time = 88.234 ms
-----------------------------------------
-----------------------------------------
| Loading reference lengths | Time = 512.78 us
-----------------------------------------
-----------------------------------------
| Loading mphf table | Time = 613.23 ms
-----------------------------------------
size = 3779163173
Number of ones: 37226610
Number of ones per inventory item: 512
Inventory entries filled: 72709
-----------------------------------------
| Loading contig boundaries | Time = 6.3706 s
-----------------------------------------
size = 3779163173
-----------------------------------------
| Loading sequence | Time = 547.49 ms
-----------------------------------------
size = 2662364873
-----------------------------------------
| Loading positions | Time = 5.0192 s
-----------------------------------------
size = 3457937385
-----------------------------------------
| Loading reference sequence | Time = 429.89 ms
-----------------------------------------
-----------------------------------------
| Loading reference accumulative lengths | Time = 1.8927 ms
-----------------------------------------
[2025-10-03 10:29:51.600] [jointLog] [info] done
[2025-10-03 10:29:51.620] [jointLog] [info] Index contained 227257 targets
[2025-10-03 10:29:51.678] [jointLog] [info] Number of decoys : 194
[2025-10-03 10:29:51.678] [jointLog] [info] First decoy index : 227024 




[2025-10-03 10:29:51.921] [jointLog] [info] Automatically detected most likely library type as IU

[A[32mprocessed[31m 500000 [32mfragments[0m
hits: 2042234, hits per frag:  4.16263







[2025-10-03 10:29:54.785] [jointLog] [info] Computed 86450 rich equivalence classes for further processing
[2025-10-03 10:29:54.785] [jointLog] [info] Counted 467828 total reads in the equivalence classes 
[2025-10-03 10:29:54.787] [jointLog] [info] Number of mappings discarded because of alignment score : 146500
[2025-10-03 10:29:54.787] [jointLog] [info] Number of fragments entirely discarded because of alignment score : 11702
[2025-10-03 10:29:54.787] [jointLog] [info] Number of fragments discarded because they are best-mapped to decoys : 6322
[2025-10-03 10:29:54.787] [jointLog] [info] Number of fragments discarded because they have only dovetail (discordant) mappings to valid targets : 3878
[2025-10-03 10:29:54.819] [jointLog] [warning] Only 467828 fragments were mapped, but the number of burn-in fragments was set to 5000000.
The effective lengths have been computed using the observed mappings.

[2025-10-03 10:29:54.819] [jointLog] [info] Mapping rate = 93.5656%

[2025-10-03 10:29:54.819] [jointLog] [info] finished quantifyLibrary()
[2025-10-03 10:29:55.149] [jointLog] [info] Starting optimizer
[2025-10-03 10:29:55.282] [jointLog] [info] Marked 0 weighted equivalence classes as degenerate
[2025-10-03 10:29:55.309] [jointLog] [info] iteration = 0 | max rel diff. = 203.028
[2025-10-03 10:29:57.974] [jointLog] [info] iteration = 100 | max rel diff. = 19.6442
[2025-10-03 10:30:00.719] [jointLog] [info] iteration = 200 | max rel diff. = 0.462189
[2025-10-03 10:30:03.319] [jointLog] [info] iteration = 300 | max rel diff. = 1.92755
[2025-10-03 10:30:05.909] [jointLog] [info] iteration = 400 | max rel diff. = 0.0347574
[2025-10-03 10:30:07.538] [jointLog] [info] iteration = 463 | max rel diff. = 0.00959217
[2025-10-03 10:30:07.581] [jointLog] [info] Finished optimizer
[2025-10-03 10:30:07.581] [jointLog] [info] writing output 

