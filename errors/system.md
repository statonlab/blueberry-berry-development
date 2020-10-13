## errors when running LoReAn
the program break at the trinity part. 
```text
salmon: error while loading shared libraries: liblzma.so.0: cannot open shared object file: No such file or directory
```
Search about the error, most answers said 'liblzma.so.0 should be supplied and owned by the xz-libs package.'. So get out of conda environment and try again.

Still getting the same error...but when running salmon independently in the container, it doesn't show an error. 

Try entering container first and then run the command
```bash
# enter a screen
screen
# launch the container
/staton/software/singularity/bin/singularity shell -B /staton/projects/blueberry_fruit_development/ -B /staton/projects/blueberry_fruit_development/Libraries/config:/opt/LoReAn/third_party/software/augustus/config/ lorean_latest.sif
# run the test pipeline
cd /staton/projects/blueberry_fruit_development/Lorean_test2/
lorean \
-rp /staton/projects/blueberry_fruit_development/annotation_20May/1_repeatMask/v1.2_mask/Vdarrowii_genome.v1.2.fasta.masked.out.gff \
-sr /staton/projects/blueberry_fruit_development/annotation_20May/2_RNAmapping/bam/S_10.Aligned.sortedByCoord.out.bam \
-lr /staton/projects/blueberry_fruit_development/raw_data/IsoSeq/all_Isoseq.fasta \
-pr /staton/projects/blueberry_fruit_development/Vcorymbosum/V_corymbosum_Draper_v1.0-proteins.fasta \
-sp darrowii \
/staton/projects/blueberry_fruit_development/Lorean_Sep20/Vdarrowii_genome.v1.2.fasta \
-d \
-t 10 \
--keep_tmp \
--prefix_gene VacDa &
```
Trinity still breaks and still get the salmon shared library error
```txt
We are sorry, commands in file: [FailedCommands] failed.  :-(

Error, cmd: /opt/LoReAn/third_party/software/Trinity/trinity-plugins/BIN/ParaFly -c trinity_GG.cmds -CPU 8 -v -shuffle  died with ret 256 at /opt/LoReAn/third_party/software/Trinity/Trinity line 2745.
	main::process_cmd("/opt/LoReAn/third_party/software/Trinity/trinity-plugins/BIN/"...) called at /opt/LoReAn/third_party/software/Trinity/Trinity line 3414
	main::run_partitioned_cmds("trinity_GG.cmds") called at /opt/LoReAn/third_party/software/Trinity/Trinity line 3372
	main::run_genome_guided_Trinity() called at /opt/LoReAn/third_party/software/Trinity/Trinity line 1367
```
Try to run ParaFly alone to see if it reproduced the error
```bash
cd LoReAn_annotation/run/Trinity/trinity_out_dir
/opt/LoReAn/third_party/software/Trinity/trinity-plugins/BIN/ParaFly -c trinity_GG.cmds -CPU 8 -v -shuffle
```
It runs without error! 
