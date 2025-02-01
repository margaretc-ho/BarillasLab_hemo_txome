#!/bin/bash
#####################################################################
### Add = (Stringtie exact match) class code transcripts from AGAPnotCodAn
#####################################################################

v7_dir=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/comparison_of_transcriptomes_v7
Agamhemo_strandfilteredtranscripts_gtf=${v7_dir}/ref_guided_transcriptome_with_percent_identity_collapsing.strandfiltered.gtf
### note this is the same as v6 version, we are just adding further analysis

reference_CDS_fa=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/reference_genome/VectorBase-59_AgambiaePEST_AnnotatedCDSs.fasta

AGAPnotCodAn_dir=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/comparison_of_transcriptomes_v7/AGAPnotCodAn
AGAPnotCodAn_classcodes=${AGAPnotCodAn_dir}/AGAPnotCodAn_transcripts.classcodes.txt
#AGAPnotCodAn_fa=${AGAPnotCodAn_dir}/AGAPnotCodAn_transcripts.comprehensive.v6.gffread.fa

mkdir -p ${AGAPnotCodAn_dir}/exactclass/
cd ${AGAPnotCodAn_dir}/exactclass
grep "=" ${AGAPnotCodAn_classcodes} | cut -f 1 > AGAPnotCodAn.exactclass.tcons.txt

### convert TCONS to AGAP from Stringtie

############### cannot use this one directly because it doesn't give detailed enough info about AGAP number - need precise AGAP transcript ID -- but we use it later to find which of the 101 = transcripts we are missing
grep -f AGAPnotCodAn.exactclass.tcons.txt ${Agamhemo_strandfilteredtranscripts_gtf} | grep -w "transcript" | cut -f 9 | cut -d ";" -f 1,3 | tr -d " " | tr -d "\"" | sed s'/gene_name//g' | sed s'/transcript_id//g' | sed s'/;/\t/g' > ${AGAPnotCodAn_dir}/exactclass/find_missing/AGAPnotCodAn.exactclass.tconsAGAP.txt
cut -f 2 ${AGAPnotCodAn_dir}/exactclass/find_missing/AGAPnotCodAn.exactclass.tconsAGAP.txt > ${AGAPnotCodAn_dir}/exactclass/find_missing/AGAPnotCodAn.exactclass.AGAPwewant.txt
#####################################################################

##### get TCONS to precise AGAP transcript id (from the cmp_ref) from
##### from v3 GTF of transcriptome (before cd-hit filtering, strand filtering)

# v3_dir=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/archivedversions/v3/comparison_of_transcriptomes_v3
# v3_gtf=${v3_dir}/comparison_of_transcriptomes.combined.gtf
# mkdir -p ${v3_dir}/exactclass

# grep AGAP ${v3_gtf} | grep -w transcript | grep "="  |  grep "cmp_ref" | cut -f 9 | cut -d ";" -f 1,5 | grep -v "contained" | tr -d " " | tr -d "\"" | sed s'/transcript_id//g' | sed s'/cmp_ref//g' | sed s'/;/\t/g' > ${v3_dir}/exactclass/v3_exactclass.tconsAGAP-precise.txt

# grep -f AGAPnotCodAn.exactclass.tcons.txt ${v3_dir}/exactclass/v3_exactclass.tconsAGAP-precise.txt > AGAPnotCodAn.exactclass.tconsAGAP-precise.txt

cut -f 2  AGAPnotCodAn.exactclass.tconsAGAP-precise.txt > AGAPnotCodAn.exactclass.AGAP-precise.txt

### extract CDS fasta using AGAP number
### this gets 98 out of 101
module purge
module load seqtk
seqtk subseq ${reference_CDS_fa} AGAPnotCodAn.exactclass.AGAP-precise.txt > AGAPnotCodAn.exactclass.tconsAGAP-precise.refCDS.fa

### grab missing AGAPnotCodan exact class transcript TCONS
mkdir -p ${AGAPnotCodAn_dir}/exactclass/find_missing
grep ">" AGAPnotCodAn.exactclass.tconsAGAP-precise.refCDS.fa | cut -d "|" -f 1 | tr -d ">"  | tr -d " " | sed 's/-[A-Z]*//g' > ${AGAPnotCodAn_dir}/exactclass/find_missing/AGAPnotCodAn.exactclass.AGAPwehave.txt

grep -v -f ${AGAPnotCodAn_dir}/exactclass/find_missing/AGAPnotCodAn.exactclass.AGAPwehave.txt ${AGAPnotCodAn_dir}/exactclass/find_missing/AGAPnotCodAn.exactclass.AGAPwewant.txt > ${AGAPnotCodAn_dir}/exactclass/find_missing/missingAGAP.txt

grep -f missingAGAP.txt ${reference_CDS_fa}

### Two were in ${reference_CDS_fa} after all, but the AGAP numbers I was looking up were missing the "precise" transcript ending (e.g. -RB).  Did manual VeuPathDB lookup and the rest of the missing turns out to be rRNA, ncRNA
### Manually added the CDS to fasta and the correct fasta header from ${reference_CDS_fa}

missingCDS_fa=${AGAPnotCodAn_dir}/exactclass/find_missing/missing_AGAPnotCodAn.exactclass.AGAP.fa

cat ${AGAPnotCodAn_dir}/exactclass/find_missing/missing_AGAPnotCodAn.exactclass.AGAP.fa AGAPnotCodAn.exactclass.tconsAGAP-precise.refCDS.fa > AGAPnotCodAn.exactclass.tconsAGAP-precise.refCDS_to_add.fa

### Manually edited missingAGAP.txt into ncRNA_AGAP.txt to remove the 2 protein coding transcripts
### and renamed ncRNA_AGAP so I can note the handful of known ncRNAs that we found

### Note that I fixed AGAPnotCodAn.exactclass.tconsAGAP-precise.txt
### to add the two TCONS and precise AGAP numbers of the two missing protein coding genes
### Searched two AGAP genes from missing_AGAPnotCodAn.exactclass.AGAP.fa
### in AGAPnotCodAn.exactclass.tcons.txt and added them to  AGAPnotCodAn.exactclass.tconsAGAP-precise.txt