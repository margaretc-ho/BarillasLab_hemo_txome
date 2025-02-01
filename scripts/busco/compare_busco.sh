#!/bin/bash
#####################################################################
### Compare DENOVO REFGUIDED REFERENCE BUSCO ###
#####################################################################

module purge
module load gffread
module load Anaconda3/2022.05
source activate busco

agam_hemodir=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes
v6_dir=${agam_hemodir}/comparison_of_transcriptomes_v6

#this refguided is the transcript GTF filtered to remove transcripts with no strand assigned
refguided_fa=${v6_dir}/all_transcripts.comprehensive.v6.gffread.fa
reference_genome=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/reference_genome/VectorBase-59_AgambiaePEST_Genome.fasta

### convert gffread

denovo_gtf=${agam_hemodir}/denovo/de_novo_transcriptome_with_genome_hits_and_expression_filtering.gtf
#gffread -w ${agam_hemodir}/denovo/de_novo_transcriptome_with_genome_hits_and_expression_filtering.gffread.fa -g ${reference_genome} ${denovo_gtf}

### run BUSCO on arthropod lineage

#busco -i ${reference_genome} -o VB59_AgambiaePEST_BUSCO -l arthropoda -m transcriptome

#busco -i ${agam_hemodir}/denovo/de_novo_transcriptome_with_genome_hits_and_expression_filtering.gffread.fa -o denovo_Agamhemo_BUSCO -l arthropoda -m transcriptome

busco -i ${refguided_fa} -o refguided_Agamhemo_BUSCO -l arthropoda -m transcriptome