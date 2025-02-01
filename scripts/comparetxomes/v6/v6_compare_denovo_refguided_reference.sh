#!/bin/bash
#####################################################################
### Compare DENOVO REFGUIDED REFERENCE ###
#####################################################################

module purge
module load bedtools

agam_hemodir=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes
v6_dir=${agam_hemodir}/comparison_of_transcriptomes_v6

# filter for transcript span only
denovo_gtf=${agam_hemodir}/denovo/de_novo_transcriptome_with_genome_hits_and_expression_filtering.gtf
#grep -w transcript ${denovo_gtf} > ${agam_hemodir}/denovo/de_novo_transcriptome_with_genome_hits_and_expression_filtering.transcriptspanonly.gtf

# uses the v6 version with removal on transcripts with no strand assigned
# filter for transcript span only
reference_gtf=${agam_hemodir}/reference_genome/VectorBase-59_AgambiaePEST.gtf
#grep -w transcript ${reference_gtf} > ${agam_hemodir}/reference_genome/VectorBase-59_AgambiaePEST.transcriptspanonly.gtf

refguided_gtf=${v6_dir}/ref_guided_transcriptome_with_percent_identity_collapsing.strandfiltered.gtf
#filter for protein coding genes only AND filter for transcript span only
grep -f ${v6_dir}/Agamhemo_codAn_partial/TCONS_of_proteincodinggenes_alltx.fa.txt ${refguided_gtf} | grep -w "transcript" > ${v6_dir}/ref_guided_transcriptome_with_percent_identity_collapsing.proteincoding.transcriptspanonly.gtf

denovo_tx_gtf=${agam_hemodir}/denovo/de_novo_transcriptome_with_genome_hits_and_expression_filtering.transcriptspanonly.gtf
reference_tx_gtf=${agam_hemodir}/reference_genome/VectorBase-59_AgambiaePEST.transcriptspanonly.gtf
refguided_tx_gtf=${v6_dir}/ref_guided_transcriptome_with_percent_identity_collapsing.proteincoding.transcriptspanonly.gtf

printf "####### Starting Counts of Transcripts ####### \n"

wc -l ${denovo_tx_gtf}
wc -l ${reference_tx_gtf}
wc -l ${refguided_tx_gtf}

module load bedops
gtf2bed < ${refguided_tx_gtf} > RG_tx.transcriptspanonly.bed
gtf2bed < ${reference_tx_gtf} > RF_tx.transcriptspanonly.bed
gtf2bed < ${denovo_tx_gtf} > DN_tx.transcriptspanonly.bed

refguided_tx_gtf=refguided_tx.transcriptspanonly.bed
reference_tx_gtf=reference_tx.transcriptspanonly.bed
denovo_tx_bed=denovo_tx.transcriptspanonly.bed

### pairwise combinations of filenames

# printf "####### 98pct identity required ####### \n"

# for i in *.bed
# do
#   for j in *.bed
#   do
#     if [ "$i" \!= "$j" ]
#     then
#     i_name=$(basename $i | sed 's/_tx.transcriptspanonly.bed//g')
#     j_name=$(basename $j | sed 's/_tx.transcriptspanonly.bed//g')
#     echo "$i_name notin $j_name"
# 	bedtools intersect -f 0.98 -v -a $i -b $j | wc -l
# 	fi
#   done
# done

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

printf "####### 98pct identity required ####### \n"

for i in *.bed
do
  for j in *.bed
  do
  	for k in *.bed
  	do
  		if [ "$i" \!= "$j" ] && [ "$j" \!= "$k" ] && [ "$i" \!= "$k" ] 
    	then
    	i_name=$(basename $i | sed 's/_tx.transcriptspanonly.bed//g')
    	j_name=$(basename $j | sed 's/_tx.transcriptspanonly.bed//g')
    	k_name=$(basename $k | sed 's/_tx.transcriptspanonly.bed//g')
    	echo "$i_name intersect $j_name intersect $k_name"
		bedtools intersect -f 0.98 -u -a $i -b $j | bedtools intersect -f 0.98 -u -a - -b $k | wc -l
		fi
	done
  done
done

printf "####### 98pct identity required ####### \n"

for i in *.bed
do
  for j in *.bed
  do
  	for k in *.bed
  	do
  		if [ "$i" \!= "$j" ] && [ "$j" \!= "$k" ] && [ "$i" \!= "$k" ] 
    	then
    	i_name=$(basename $i | sed 's/_tx.transcriptspanonly.bed//g')
    	j_name=$(basename $j | sed 's/_tx.transcriptspanonly.bed//g')
    	k_name=$(basename $k | sed 's/_tx.transcriptspanonly.bed//g')
    	echo "$i_name notin $j_name notin $k_name"
		bedtools intersect -v -f 0.98 -a $i -b $j | bedtools intersect -v -f 0.98 -a - -b $k | wc -l
		fi
	done
  done
done

### pairwise combinations of filenames

# printf "####### 50pct identity required ####### \n"

# for i in *.bed
# do
#   for j in *.bed
#   do
#     if [ "$i" \!= "$j" ]
#     then
#     i_name=$(basename $i | sed 's/_tx.transcriptspanonly.bed//g')
#     j_name=$(basename $j | sed 's/_tx.transcriptspanonly.bed//g')
#     echo "$i_name notin $j_name"
# 	bedtools intersect -f 0.5 -v -a $i -b $j | wc -l
# 	fi
#   done
# done

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

printf "####### 50pct identity required ####### \n"

for i in *.bed
do
  for j in *.bed
  do
  	for k in *.bed
  	do
  		if [ "$i" \!= "$j" ] && [ "$j" \!= "$k" ] && [ "$i" \!= "$k" ] 
    	then
    	i_name=$(basename $i | sed 's/_tx.transcriptspanonly.bed//g')
    	j_name=$(basename $j | sed 's/_tx.transcriptspanonly.bed//g')
    	k_name=$(basename $k | sed 's/_tx.transcriptspanonly.bed//g')
    	echo "$i_name intersect $j_name intersect $k_name"
		bedtools intersect -f 0.5 -u -a $i -b $j | bedtools intersect -f 0.5 -u -a - -b $k | wc -l
		fi
	done
  done
done

printf "####### 50pct identity required ####### \n"

for i in *.bed
do
  for j in *.bed
  do
  	for k in *.bed
  	do
  		if [ "$i" \!= "$j" ] && [ "$j" \!= "$k" ] && [ "$i" \!= "$k" ] 
    	then
    	i_name=$(basename $i | sed 's/_tx.transcriptspanonly.bed//g')
    	j_name=$(basename $j | sed 's/_tx.transcriptspanonly.bed//g')
    	k_name=$(basename $k | sed 's/_tx.transcriptspanonly.bed//g')
    	echo "$i_name notin $j_name notin $k_name"
		bedtools intersect -v -f 0.5 -a $i -b $j | bedtools intersect -v -f 0.5 -a - -b $k | wc -l
		fi
	done
  done
done