# TitanCNA-utils
TitanCNA utility scripts

## Instructions for using titanCNA_v1.x.y.R script
**Description**  
R script will run the R component of the TITAN analysis using the TitanCNA R/Bioconductor package.  

**Input files**  
  This script assumes that the necessary input files have been generated.  
  1. GC-corrected, normalized read coverage using the HMMcopy suite  
  2. Tumour allelic read counts at heterozygous SNPs (identifed from the normal sample).

The easiest way to generate these files is by using the downloadable pipeline [KRONOS](https://github.com/MO-BCCRC/titan_workflow).  

**Running the R script**  

1. Clone the git repo and locate the folder containing the Rscript  
  ```
  git clone git@github.com:gavinha/TitanCNA-utils.git
  cd TitanCNA-utils/titan_scripts/  
  ```  
2. Look at the usage of the R script  
  ```
  # from the command line
  > Rscript titanCNA_v1.10.1.R --help
  Usage: Rscript titanCNA_v1.10.1.R [options]


  Options:
        --id=ID
                Sample ID

        --hetFile=HETFILE
                File containing allelic read counts at HET sites. (Required)

        --cnFile=CNFILE
                File containing normalized coverage as log2 ratios. (Required)

        --outDir=OUTDIR
                Output directory to output the results. (Required)

        --numClusters=NUMCLUSTERS
                Number of clonal clusters. (Default: 1)

        --numCores=NUMCORES
                Number of cores to use. (Default: 1)

        --ploidy_0=PLOIDY_0
                Initial ploidy value; float (Default: 2)

        --estimatePloidy=ESTIMATEPLOIDY
                Estimate ploidy; TRUE or FALSE (Default: TRUE)

        --normal_0=NORMAL_0
                Initial normal contamination (1-purity); float (Default: 0.5)

        --estimateNormal=ESTIMATENORMAL
                Estimate normal contamination method; string {'map', 'fixed'} (Default: map)

        --maxCN=MAXCN
                Maximum number of copies to model; integer (Default: 8)
        ...
  ```

3. Example usage of R script
  ```
  # normalized coverage file: test.cn.txt
  # allelic read count file: test.het.txt
  Rscript titanCNA_v1.10.1.R --id test --hetFile test.het.txt --cnFile test.cn.txt \
    --numClusters 1 --numCores 1 --normal_0 0.5 --ploidy_0 2 \
    --chrs "c(1:22, \"X\")" --estimatePloidy TRUE --outDir ./
  ```
  Additional arguments to consider are the following:  
    These arguments can be used to tune the model based on variance in the read coverage data and data-type (whole-exome sequencing or whole-genome sequencing).
    ```
    --alphaK=ALPHAK
                Hyperparameter on Gaussian variance; for WES, use 1000; for WGS, use 10000; 
                float (Default: 10000)

    --alphaKHigh=ALPHAKHIGH
                Hyperparameter on Gaussian variance for extreme copy number states; 
                for WES, use 1000; for WGS, use 10000; float (Default: 10000)
    ```
