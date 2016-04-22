# eQTL_v6_script

I use these scripts to process GTEx.v.6 data.


## 1. Gene expression data

Do the following two sequentially:

1. expression\_geno\_etissue.py
2. expression\_gene\_normal.py


## 2. Genotype data

In the imputed data directory ("phg000520.v2.GTEx\_MidPoint\_Imputation.genotype-calls-vcf.c1"), we have the VCF-format genotype data. Create the following folder:

1. ./genotype\_imputed
2. ./genotype\_imputed/genotype\_vcf
3. ./genotype\_imputed/genotype\_processed
4. ./genotype\_imputed/genotype\_post\_prune



Next gunzip this file: "./GTEx\_Analysis\_20150112\_OMNI\_2.5M\_5M\_450Indiv\_chr1to22\_genot\_imput\_info04\_maf01\_HWEp1E6\_ConstrVarIDs.vcf.gz".

Then run the two routines, "**genotype\_vcf\_splitter.py**" and "**genotype\_vcf\_parser.py**" under "./genotype\_imputed/" sequentially (in the middle, please gzip this file to save some quota: "./GTEx\_Analysis\_20150112\_OMNI\_2.5M\_5M\_450Indiv\_chr1to22\_genot\_imput\_info04\_maf01\_HWEp1E6\_ConstrVarIDs.vcf"), to perform the following two functions:

1. split the VCF file into header file and sub-VCF file for each chromosome (saved to "./genotype\_imputed/genotype\_vcf/")
2. extract the tped/dosage/snp\_exclusion\_list information from each sub-VCF file for each chromosome (saved to "./genotype\_imputed/genotype\_processed/"), and delete the vcf file after processing the current chromosome (to save some disk quota).

**NOTE:** There are duplicated SNPs in chr1, chr2 and chr12 and potentially others (observed; for example, "12\_48000000\_T\_C\_b37" in chr12). So I need to check and remove duplications in the below processing script before QC and LD. The duplicated SNPs are documented below for GTEx's reference (some of them are only duplicated in names, and they might be different SNPs actually IoI; I removed the second one and left the first appearing one):


| Chromosome        | duplicated SNP           |
| ------------- |:-------------:|
| Chr1      | 1\_105000000\_C\_CATT\_b37 |
| Chr2      | 2\_72000000\_G\_A\_b37, 2\_87000000\_CA\_C\_b37      |
| Chr10 | 10\_12000000\_A\_T\_b37      |
| Chr12 | 12\_48000000\_T\_C\_b37      |
| Chr14 | 14\_81000000\_TC\_T\_b37      |


We can find the tfam file under "../phg000520.v2.GTEx\_MidPoint.genotype-qc.MULTI/Imputation/", and we need to create "chrX.tfam" under "./genotype\_imputed/genotype\_processed/" as further processing requires.

After that, we can use the procedure similar to GTEx.v.4 to process the data (snp QC exclusion, pruning; here is the [link](https://github.com/morrisyoung/eQTL_v4_script#5-the-pipeline-for-genotype-qc-and-ld-pruning)).

**Note:** We need "chrX.tped" and "chrX.tfam" (--tfile) to perform the snp QC exclusion (witn snp\_exclusion\_list), and that will generate bed/fam/bim files that we can do the LD pruning. And we can finally extract the dosage information from "chrX.dosage".

Specifically, do the following (updated from the pipeline for last version data):

"genotype\_qc\_ld.py" and "genotype\_post\_prune.py", and the pipeline is as followed:

```
1. for all the chromosome “X”, do all the following (step #2 - #10):
2. cp chrX.tfam “chrX.tfam”
3. [QC] ./plink --tfile “chrX” --exclude “chrX.maf05.exclusion.snplist.txt” --make-bed
4. [pruning] ./plink --bfile plink --indep-pairwise 50kb 5 0.5 (0.5 may possibly be adjusted into other values)
5. [LD statistics] ./plink --bfile plink --r2 --ld-snp-list plink.prune.out --ld-window-kb 50 --ld-window 99999 --ld-window-r2 0.5 (0.5 can be adjusted into other values)
6. [LD statistics further] python genotype_post_prune.py
7. mv plink.prune.in ../genotype_post_prune/chr“X”.prune.in
8. mv plink.prune.out ../genotype_post_prune/chr“X”.prune.out
9. mv chr“X”.post_prune.txt ../genotype_post_prune/
10. rm plink.*
```

The work is done under "/ifs/scratch/c2b2/ip\_lab/sy2515/GTEx/data.v.6/47024/PhenoGenotypeFiles/RootStudyConsentSet\_phs000424.GTEx.v6.p1.c1.GRU/GenotypeFiles/phg000520.v2.GTEx\_MidPoint\_Imputation.genotype-calls-vcf.c1/genotype\_imputed/genotype\_processed/", and the results will be moved to "/ifs/scratch/c2b2/ip\_lab/sy2515/GTEx/data.v.6/47024/PhenoGenotypeFiles/RootStudyConsentSet\_phs000424.GTEx.v6.p1.c1.GRU/GenotypeFiles/phg000520.v2.GTEx\_MidPoint\_Imputation.genotype-calls-vcf.c1/genotype\_imputed/genotype\_post\_prune/".


Finally, we can extract the left SNPs (un-pruned) from the dosage file, wiht script "genotype\_dosage\_matrix\_qc\_ld.py". This script will work under "/ifs/scratch/c2b2/ip\_lab/sy2515/GTEx/data.v.6/47024/PhenoGenotypeFiles/RootStudyConsentSet\_phs000424.GTEx.v6.p1.c1.GRU/GenotypeFiles/phg000520.v2.GTEx\_MidPoint\_Imputation.genotype-calls-vcf.c1/genotype\_imputed/", read dosage data from "./genotype\_processed/", and pruning information from "./genotype\_post\_prune/", and output the processed dosage matrix data in "./genotype\_450\_dosage\_matrix\_qc/" (organized by chromosomes and individuals). This script will additionally requires following data:

1. chrX.tfam file, "./genotype\_processed/chrX.tfam", as we need to know in dosage file how these individuals are ordered


