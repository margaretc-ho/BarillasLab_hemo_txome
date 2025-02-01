#!/bin/bash
##################################################################### Investigating potential ncRNAs
###############################################################

v7_dir=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/comparison_of_transcriptomes_v7
Agamhemo_strandfilteredtranscripts_gtf=${v7_dir}/ref_guided_transcriptome_with_percent_identity_collapsing.strandfiltered.gtf
### note this is the same as v6 version, we are just adding further analysis
all_tx_gff_read=${v7_dir}/all_transcripts.comprehensive.v6.gffread.fa
reference_genome=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/reference_genome/VectorBase-59_AgambiaePEST_Genome.fasta

# grep ">" ${v7_dir}/notCodAn/notCodan_v7_addtl_proteincodingtx.fa | sed 's/.p.*//g' | tr -d '>' | cut -d ' ' -f 1 > ${v7_dir}/notCodAn/notCodan_v7_addtl_proteincodingtx.tcons.txt

# grep -v -f ${v7_dir}/notCodAn/notCodan_v7_addtl_proteincodingtx.tcons.txt ${v7_dir}/notCodAn/notCodan_transcripts.txt > ${v7_dir}/notCodAn/notCodan_v7_noncoding.tcons.txt

# grep -f ${v7_dir}/notCodAn/notCodan_v7_noncoding.tcons.txt ${v7_dir}/notCodAn/ref_guided_idcollapse.strandfiltered.notCodAn.gtf > ${v7_dir}/notCodAn/notCodan_v7_noncoding.gtf

# module purge
# module load gffread

# gffread -w ${v7_dir}/notCodAn/notCodan_v7_noncoding.gffread.fa -g ${reference_genome} ${v7_dir}/notCodAn/notCodan_v7_noncoding.gtf

##################################################################### Run RNAMMER on potential ncRNAs
###############################################################

module purge
module load Anaconda3/2022.05
source activate infernal

rfam_cm="/hpcdata/bcbb/homc/conda_envs/envs/infernal/databases/Rfam.cm"
rfam_clanin="/hpcdata/bcbb/homc/conda_envs/envs/infernal/databases/Rfam.clanin"
v7_Agam_noncoding_fa=${v7_dir}/notCodAn/notCodan_v7_noncoding.gffread.fa

## Template command
## cmscan --nohmmonly --rfam --cut_ga --fmt 2 --oclan --oskip --clanin -o my.cmscan.out --tblout my.cmscan.tblout Rfam.cm my.fa

###Use the cmscan program to identify all subsequences that match any Rfam family with a score above the gathering cutoff (GA) selected by the Rfam curators:

# mkdir -p ${v7_dir}/notCodAn/noncoding_rnammer

# cmscan --nohmmonly --rfam --cut_ga --fmt 2 --oclan --oskip --clanin ${rfam_clanin} -o ${v7_dir}/notCodAn/noncoding_rnammer/v7_Agam_noncoding_fa.cmscan.out --tblout ${v7_dir}/notCodAn/noncoding_rnammer/v7_Agam_noncoding_fa.cmscan.tblout ${rfam_cm} ${v7_Agam_noncoding_fa}

####tblout is tabular output file that lists info on each hit above the GA cutoff for all families, one line per hit
#the set of hits in the tblout file will no include any hits that overlap in the same sequence between models in the same clan
#this lack of overlapping hits is one potential advantage of using cmscan on command line instead of RFAM website search

#the my.cmscan.out file includes more info such as a alignments of all hits to query models including overlaps that may have been removed in the tblout file

# grep -o "TCONS.*" ${v7_dir}/notCodAn/noncoding_rnammer/v7_Agam_noncoding_fa.cmscan.tblout | tr -s ' ' | cut -d ' ' -f 1 | sort -u > ${v7_dir}/notCodAn/noncoding_rnammer/v7_Agam_noncoding_fa.cmscan.tblout.tcons.txt


##################################################################### Noncoding - Potential lncRNAs
###############################################################

#mkdir -p ${v7_dir}/notCodAn/lncRNAs/

#grep -v -f ${v7_dir}/notCodAn/noncoding_rnammer/v7_Agam_noncoding_fa.cmscan.tblout.tcons.txt ${v7_dir}/notCodAn/notCodan_v7_noncoding.tcons.txt > ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_maybe_lncRNA.tcons.txt

classcodes=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/classcodes/TCONS_to_gffcompare_class_code_relative_to_AgamP4.14.txt

# grep -f ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_maybe_lncRNA.tcons.txt ${classcodes} | cut -f 1,3 > ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_maybe_lncRNA.tcons.classcodes.txt

##################################################################### Retrieve 2 AGAP = genes (add to known ncRNA list)
#################################################################

#### I have seen that there are two rRNA genes in this dataset, let's remove them
#### [homc@ai-hpcn038 lncRNAs]$ cut -f 2 notCodan_v7_maybe_lncRNA.tcons.classcodes.txt | sort | uniq -c | sort -nk 1
####       2 =
####       5 s
####      20 k
####      24 y
####      42 n
####      61 m
####      85 e
####     108 p
####     151 c
####     218 x
####     257 i
####     283 o
####     333 j
####    1814 u

####[homc@ai-hpcn038 lncRNAs]$ grep "=" notCodan_v7_maybe_lncRNA.tcons.classcodes.txt
####TCONS_00050427  =
####TCONS_00050428  =

### This corresponds to 
### TCONS_00050427 AGAP028391 lsu rRNA
### TCONS_00050428 AGAP028393 ssu rRNA

### Moved these to known ncRNA files

##################################################################### Remaining lncRNA
#################################################################

####    1814 u (keep)
####     257 i
####     218 x
#### total = 2289

# grep "u\|i\|x" ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_maybe_lncRNA.tcons.classcodes.txt > ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.tcons.class.u.i.x.txt

#### ambiguous (keep)
####      20 k
####      24 y
####      42 n
####      61 m
####     151 c
####     283 o
#### total = 581

# grep "k\|y\|n\|m\|c\|o" ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_maybe_lncRNA.tcons.classcodes.txt > ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.tcons.class.k.y.n.m.c.o.txt

### likely can drop
####      85 e
####     108 p
####       5 s
#### total = 198

grep "e\|p\|s" ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_maybe_lncRNA.tcons.classcodes.txt > ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.tcons.class.e.p.s.txt

### 333 j class noncoding RNAs

#grep "j" ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_maybe_lncRNA.tcons.classcodes.txt > ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.tcons.class.j.txt

#### 3203 (2289 + 581 + 333) likely lncRNAs 

#cat ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.tcons.class.u.i.x.txt ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.tcons.class.k.y.n.m.c.o.txt ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.tcons.class.j.txt > ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.tcons.class.uixkynmcoj.txt


#cut -f 1  ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.tcons.class.uixkynmcoj.txt >  ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.uixkynmcoj.tconsonly.txt

#grep -f ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.uixkynmcoj.tconsonly.txt ${v7_dir}/notCodAn/ref_guided_idcollapse.strandfiltered.notCodAn.gtf > ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.uixkynmcoj.tconsonly.gtf

# module purge
# module load gffread

# gffread -w ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.uixkynmcoj.gffread.fa -g ${reference_genome} ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.uixkynmcoj.tconsonly.gtf

### Fasta length distribution

module purge
module load samtools

### run samtools faidx
### http://www.htslib.org/doc/faidx.html
### second column has the fasta length

#samtools faidx ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.uixkynmcoj.gffread.fa -o ${v7_dir}/notCodAn/lncRNAs/notCodan_v7_lncRNA.uixkynmcoj.gffread.fa.fai