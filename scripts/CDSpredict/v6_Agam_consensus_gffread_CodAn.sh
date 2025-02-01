#!/bin/bash
#####################################################################
### AGAM_v6 GFFread and CodAn ###
#####################################################################

### v5 version corrects the mistake I made where I provided only full length transcript data to CodAn and not exon information
### Colton helped me identify the mistake and here, gffread is run on the full GTF information (transcript, exons, UTR) and then resulting gffread.fa is provided to CodAn

### v6 version adds a filter to remove transcripts from the Stringtie GTF that have strand "."
### "StringTie refuses to consider spliced alignments as valid in cases where the splice sites are very inconsistent", between reads, which can be the case with long reads."
### Usually resulting from low read coverage. In summary, its difficult to assign a strand to these transcripts and it is likely they aren't real.

module purge
module load Anaconda3/2022.05
source activate codan_env
export PATH=$PATH:/hpcdata/bcbb/homc/conda_envs/envs/codan_env/CodAn/bin/
#module load gffread

v6_dir=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/comparison_of_transcriptomes_v6
v6_Agamhemo_consensustranscripts_gtf=${v6_dir}/v5_cd-hit_notyetfilteredforstrand//re-quantified_expression_after_percent_identity_filtering/ref_guided_transcriptome_with_percent_identity_collapsing.gtf
alltx_codan_dir=${v6_dir}/Agamhemo_codAn_partial

mkdir -p ${alltx_codan_dir}
# ### GFF read and CodAn ###
reference_genome=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/reference_genome/VectorBase-59_AgambiaePEST_Genome.fasta
partial_model=/hpcdata/bcbb/homc/conda_envs/envs/codan_env/CodAn/models/INV_partial

cd   ${v6_dir}
### filter Stringtie GTF to remove transcripts that have strand "." which means strand is unassigned
### added requirement that codan use the plus strand on the fasta to predict (since we are following the strand convention from the GTF file -> gffread fasta)

# awk '$7!="." {print $0}' ${v6_Agamhemo_consensustranscripts_gtf} >  ref_guided_transcriptome_with_percent_identity_collapsing.strandfiltered.gtf
# gffread -w all_transcripts.comprehensive.v6.gffread.fa -g ${reference_genome} ref_guided_transcriptome_with_percent_identity_collapsing.strandfiltered.gtf
# codan.py -c 3 -s plus -t ${v6_dir}/all_transcripts.comprehensive.v6.gffread.fa -o ${alltx_codan_dir} -m ${partial_model} | tee ${alltx_codan_dir}/codan_partial.log

# # ### CodAn Translate Partial model ###
TranslatePartialpy=/hpcdata/bcbb/homc/conda_envs/envs/codan_env/CodAn/scripts/TranslatePartial.py
python ${TranslatePartialpy} ${alltx_codan_dir}/ORF_sequences.fasta ${alltx_codan_dir}/ORF_sequences.pep

### Minimap CodAn transcripts (all) back to A.gam reference genome ###
# module purge
# module load minimap2
# module load samtools

# minimap2 -ax splice:hq --secondary=no ${reference_genome} ${alltx_codan_dir}/ORF_sequences.fasta > ${alltx_codan_dir}/codan_alltx_partialORF_align_VB59_Agam.primaryalignmentonly.sam
# samtools view -S -b ${alltx_codan_dir}/codan_alltx_partialORF_align_VB59_Agam.primaryalignmentonly.sam | samtools sort -O BAM -o ${alltx_codan_dir}/codan_alltx_partialORF_align_VB59_Agam.primaryalignmentonly.bam

# module load Anaconda3/2022.05
# source activate agat1
# agat_convert_minimap2_bam2gff.pl -i ${alltx_codan_dir}/codan_alltx_partialORF_align_VB59_Agam.primaryalignmentonly.sam -o ${alltx_codan_dir}/codan_alltx_partialORF_align_VB59_Agam.primaryalignmentonly.gff

## record IDs (TCONS) for all protein coding transcripts found by CodAn
#cd ${alltx_codan_dir}
#grep ">" ORF_sequences.fasta | sed 's/>//g' | sort > TCONS_of_proteincodinggenes_alltx.fa.txt
#grep ">" ORF_sequences.pep | sed 's/>//g' | sed 's/frame[0-9]*//g' | sort | uniq -c | sort > TCONS_of_proteincodinggenes_alltx.pep.txt