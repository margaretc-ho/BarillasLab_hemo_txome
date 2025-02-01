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

### subset Stringtie GTF to investigate not AGAP, not CodAn transcripts
### (roughly 2907 of them)

mkdir -p ${v7_dir}/notCodAn

grep -v -f ${v7_dir}/Agamhemo_codAn_partial/TCONS_of_proteincodinggenes_alltx.fa.txt  ${Agamhemo_strandfilteredtranscripts_gtf} > ${v7_dir}/notCodAn/ref_guided_idcollapse.strandfiltered.notCodAn.gtf

grep -w "transcript" ${v7_dir}/notCodAn/ref_guided_idcollapse.strandfiltered.notCodAn.gtf | cut -f 9 | cut -f 2 -d " " | sed s'/"//g' | sed s'/;//g' > ${v7_dir}/notCodAn/notCodan_transcripts.txt