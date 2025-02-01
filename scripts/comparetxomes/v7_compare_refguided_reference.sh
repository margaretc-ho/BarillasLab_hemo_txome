#!/bin/bash
#####################################################################
### Compare REFGUIDED HEMO vs REFERENCE ###
#####################################################################

module purge
module load bedtools

agam_hemodir=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes
v7_dir=${agam_hemodir}/comparison_of_transcriptomes_v7

# uses the v7 version with removal on transcripts with no strand assigned
# filter for transcript span only
reference_gtf=${agam_hemodir}/reference_genome/VectorBase-59_AgambiaePEST.gtf
#grep -w transcript ${reference_gtf} > ${agam_hemodir}/reference_genome/VectorBase-59_AgambiaePEST.transcriptspanonly.gtf

refguided_gtf=${v7_dir}/ref_guided_transcriptome_with_percent_identity_collapsing.strandfiltered.gtf
#filter for protein coding genes only AND filter for transcript span only

pctx1_tcons_list=${v7_dir}/Agamhemo_codAn_partial/TCONS_of_proteincodinggenes_alltx.fa.txt 
pctx2_tcons_list=${v7_dir}/notCodAn/notCodan_v7_addtl_proteincodingtx.tcons.txt

### combines the two lists of protein coding genes #1 from Codan #2 of additional protein coding genes IDed by exact class and transdecoder
cat ${pctx1_tcons_list} ${pctx2_tcons_list} > ${v7_dir}/TCONS_of_proteincodingtx_CodAn_v7_addtl.tcons.txt
grep -f ${v7_dir}/TCONS_of_proteincodingtx_CodAn_v7_addtl.tcons.txt ${refguided_gtf} | grep -w "transcript" > ${v7_dir}/v7_ref_guided_tx_with_pct_identity_collapsing.proteincoding.txspanonly.gtf

reference_tx_gtf=${agam_hemodir}/reference_genome/VectorBase-59_AgambiaePEST.transcriptspanonly.gtf
refguided_tx_gtf=${v7_dir}/v7_ref_guided_tx_with_pct_identity_collapsing.proteincoding.txspanonly.gtf

printf "####### Starting Counts of Transcripts ####### \n"


wc -l ${reference_tx_gtf}
wc -l ${refguided_tx_gtf}

module load bedops
gtf2bed < ${refguided_tx_gtf} > RG_tx.transcriptspanonly.bed
gtf2bed < ${reference_tx_gtf} > RF_tx.transcriptspanonly.bed

refguided_tx_gtf=refguided_tx.transcriptspanonly.bed
reference_tx_gtf=reference_tx.transcriptspanonly.bed

# ### pairwise combinations of filenames

printf "####### 98pct identity required ####### \n"

for i in *.bed
do
  for j in *.bed
  do
    if [ "$i" \!= "$j" ]
    then
    i_name=$(basename $i | sed 's/_tx.transcriptspanonly.bed//g')
    j_name=$(basename $j | sed 's/_tx.transcriptspanonly.bed//g')
    echo "$i_name notin $j_name"
	bedtools intersect -f 0.98 -v -a $i -b $j | wc -l
	fi
  done
done

printf "####### 98pct identity required ####### \n"

for i in *.bed
do
  for j in *.bed
  do
    if [ "$i" \!= "$j" ]
    then
    i_name=$(basename $i | sed 's/_tx.transcriptspanonly.bed//g')
    j_name=$(basename $j | sed 's/_tx.transcriptspanonly.bed//g')
    echo "$i_name intersect $j_name"
	bedtools intersect -f 0.98 -u -a $i -b $j | wc -l
	fi
  done
done


### pairwise combinations of filenames

printf "####### 50pct identity required ####### \n"

for i in *.bed
do
  for j in *.bed
  do
    if [ "$i" \!= "$j" ]
    then
    i_name=$(basename $i | sed 's/_tx.transcriptspanonly.bed//g')
    j_name=$(basename $j | sed 's/_tx.transcriptspanonly.bed//g')
    echo "$i_name notin $j_name"
	bedtools intersect -f 0.5 -v -a $i -b $j | wc -l
	fi
  done
done

printf "####### 50pct identity required ####### \n"

for i in *.bed
do
  for j in *.bed
  do
    if [ "$i" \!= "$j" ]
    then
    i_name=$(basename $i | sed 's/_tx.transcriptspanonly.bed//g')
    j_name=$(basename $j | sed 's/_tx.transcriptspanonly.bed//g')
    echo "$i_name intersect $j_name"
	bedtools intersect -f 0.5 -u -a $i -b $j | wc -l
	fi
  done
done
