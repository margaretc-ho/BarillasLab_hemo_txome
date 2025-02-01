#!/bin/bash
#####################################################################
### AGAM_v6 GFFread and CodAn ###
#####################################################################

### v5 version corrects the mistake I made where I provided only full length transcript data to CodAn and not exon information
### Colton helped me identify the mistake and here, gffread is run on the full GTF information (transcript, exons, UTR) and then resulting gffread.fa is provided to CodAn

### v6 version adds a filter to remove transcripts from the Stringtie GTF that have strand "."
### "StringTie refuses to consider spliced alignments as valid in cases where the splice sites are very inconsistent", between reads, which can be the case with long reads."
### Usually resulting from low read coverage. In summary, its difficult to assign a strand to these transcripts and it is likely they aren't real.

### what are the non-codan genes?
### that were not predicted as protein coding by CodAn

v7_dir=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/comparison_of_transcriptomes_v7
Agamhemo_strandfilteredtranscripts_gtf=${v7_dir}/ref_guided_transcriptome_with_percent_identity_collapsing.strandfiltered.gtf
### note this is the same as v6 version, we are just adding further analysis

reference_genome=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/reference_genome/VectorBase-59_AgambiaePEST_Genome.fasta

### subset Stringtie GTF to investigate AGAP yes, not CodAn transcripts

mkdir -p ${v7_dir}/AGAPnotCodAn

#grep -v -f ${v7_dir}/Agamhemo_codAn_partial/TCONS_of_proteincodinggenes_alltx.fa.txt  ${Agamhemo_strandfilteredtranscripts_gtf} | grep AGAP > ${v7_dir}/AGAPnotCodAn/ref_guided_idcollapse.strandfiltered.AGAPnotCodAn.gtf

#grep -w "transcript" ${v7_dir}/AGAPnotCodAn/ref_guided_idcollapse.strandfiltered.AGAPnotCodAn.gtf | cut -f 9 | cut -f 2 -d " " | sed s'/"//g' | sed s'/;//g' > ${v7_dir}/AGAPnotCodAn/AGAPnotCodAn_transcripts.txt

#### What Class Codes?
#### Note we cannot rely on the AGAP data from this because it is old version of transcriptome
classcodes=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/classcodes/TCONS_to_gffcompare_class_code_relative_to_AgamP4.14.txt

#grep -f ${v7_dir}/AGAPnotCodAn/AGAPnotCodAn_transcripts.txt ${classcodes} | cut -f 1,3 > ${v7_dir}/AGAPnotCodAn/AGAPnotCodAn_transcripts.classcodes.txt

# module purge
# module load gffread
# gffread -w ${v7_dir}/AGAPnotCodAn/AGAPnotCodAn_transcripts.comprehensive.v6.gffread.fa -g ${reference_genome} ${v7_dir}/AGAPnotCodAn/ref_guided_idcollapse.strandfiltered.AGAPnotCodAn.gtf

### Trinotate ### rerunning to compare with CodAn
# module purge
# module load transdecoder
# TransDecoder.LongOrfs -S -t ${v7_dir}/AGAPnotCodAn/AGAPnotCodAn_transcripts.comprehensive.v6.gffread.fa --output_dir ${v7_dir}/AGAPnotCodAn/transdecoder

### -S indicates sense strand prediction
# ### this generates longest_orfs.pep which is "longest_orfs.pep   : all ORFs meeting the minimum length criteria, regardless of coding potential."
# ### this is why looking at longest_orfs.pep never seemed to give me good CDS!!!
# ### need to run through TransDecoder.Predict to actually get the CDS

# ### Including homology searches as ORF retention criteria
# #### Running BLASTP 
# module purge
# module load Anaconda3/2022.05 
# source activate trinotate

TRINOTATE_DIR="/hpcdata/bcbb/homc/conda_envs/envs/trinotate/bin"
uniprotdb="${TRINOTATE_DIR}/uniprot_sprot.pep"
pfamhmm="${TRINOTATE_DIR}/Pfam-A.hmm"

#blastp -query ${v7_dir}/AGAPnotCodAn/transdecoder/longest_orfs.pep -db ${uniprotdb} -num_threads 16 -max_target_seqs 1 -outfmt 6 > ${v7_dir}/AGAPnotCodAn/transdecoder/blastp.outfmt6

#### Running PFAM HMMER
#hmmsearch --cpu 16 -E 1e-10 --domtblout ${v7_dir}/AGAPnotCodAn/transdecoder/pfam.domtblout ${pfamhmm} ${v7_dir}/AGAPnotCodAn/transdecoder/longest_orfs.pep > ${v7_dir}/AGAPnotCodAn/transdecoder/pfam.log


# ### Integrating the Blast and Pfam search results into coding region selection
module purge
module load transdecoder

cd ${v7_dir}/AGAPnotCodAn/
TransDecoder.Predict --single_best_only -t ${v7_dir}/AGAPnotCodAn/AGAPnotCodAn_transcripts.comprehensive.v6.gffread.fa --retain_pfam_hits ${v7_dir}/AGAPnotCodAn/transdecoder/pfam.domtblout --retain_blastp_hits ${v7_dir}/AGAPnotCodAn/transdecoder/blastp.outfmt6 --output_dir ${v7_dir}/AGAPnotCodAn/transdecoder/