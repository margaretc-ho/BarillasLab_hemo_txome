#!/bin/bash
##################################################################### Investigating sequence conservation of ncRNAs
###############################################################

v7_dir=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/comparison_of_transcriptomes_v7
Agamhemo_strandfilteredtranscripts_gtf=${v7_dir}/ref_guided_transcriptome_with_percent_identity_collapsing.strandfiltered.gtf
### note this is the same as v6 version, we are just adding further analysis

AgamP4_conservation_dir=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/conservation/AgamP4_conservation_score
AgamP4_conservation_bg=${AgamP4_conservation_dir}/AgamP4_phyloP.bg

### need to convert chrm names from chr2R to AgamP4_2R

#sed 's/chr/AgamP4_/g' ${AgamP4_conservation_bg} > ${AgamP4_conservation_dir}/AgamP4_phyloP.chrnamefix.bg


reference_genome=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/reference_genome/VectorBase-59_AgambiaePEST_Genome.fasta
lncRNA_dir=${v7_dir}/notCodAn/lncRNAs
#${lncRNA_dir}/notCodan_v7_lncRNA.uixkynmco.gffread.fa

lncRNA_gtf=${lncRNA_dir}/notCodan_v7_lncRNA.uixkynmco.tconsonly.gtf


#grep -w "transcript" ${lncRNA_gtf} > ${lncRNA_dir}/notCodan_v7_lncRNA.uixkynmco.tconsonly.transcriptspan.gtf

#module purge
#module load bedops

#gtf2bed < ${lncRNA_gtf} > ${lncRNA_dir}/seqconserv/notCodan_v7_lncRNA.uixkynmco.tconsonly.bed

#gtf2bed < ${lncRNA_dir}/notCodan_v7_lncRNA.uixkynmco.tconsonly.transcriptspan.gtf > ${lncRNA_dir}/seqconserv/notCodan_v7_lncRNA.uixkynmco.tconsonly.transcriptspan.bed

### simple bed files
cut -f 1,2,3 ${lncRNA_dir}/seqconserv/notCodan_v7_lncRNA.uixkynmco.tconsonly.transcriptspan.bed |sort -k1,1 -k2,2n > ${lncRNA_dir}/seqconserv/notCodan_v7_lncRNA.uixkynmco.tconsonly.transcriptspan.simple.bed


### sequence conservation 

module purge
module load bedtools

bedtools map -a ${lncRNA_dir}/seqconserv/notCodan_v7_lncRNA.uixkynmco.tconsonly.transcriptspan.simple.bed -b ${AgamP4_conservation_dir}/AgamP4_phyloP.chrnamefix.bg -c 4 -o mean > ${lncRNA_dir}/seqconserv/test.txt