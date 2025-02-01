#!/bin/bash
#####################################################################
### Additional protein coding transcripts to add to annotation (v7)
### Concatenate exactAGAPnotCodan, AGAPnotCodan_transdecoder, notAGAPnotCodan_transdecoder together
#####################################################################


v7_dir=/hpcdata/bcbb/homc/Barillaslab/Agam_hemo_transcriptomes/comparison_of_transcriptomes_v7
AGAPnotCodAn_dir=${v7_dir}/AGAPnotCodAn
notAGAPnotCodAn_dir=${v7_dir}/notAGAPnotCodAn
notCodAn_dir=${v7_dir}/notCodAn

# ## concatenate refCDS
# exactAGAPnotCodan_fa=${AGAPnotCodAn_dir}/exactclass/AGAPnotCodAn.exactclass.tconsAGAP-precise.refCDS_to_add.ids.fa

# AGAPnotCodan_transdecoder_fa=${AGAPnotCodAn_dir}/AGAPnotCodAn_transcripts.comprehensive.v6.gffread.fa.transdecoder.filt.cds

# notAGAPnotCodan_transdecoder_fa=${notAGAPnotCodAn_dir}/notAGAPnotCodAn_transcripts.comprehensive.v6.gffread.fa.transdecoder.cds

# cat ${exactAGAPnotCodan_fa} ${AGAPnotCodan_transdecoder_fa} ${notAGAPnotCodan_transdecoder_fa} > ${notCodAn_dir}/notCodan_v7_addtl_proteincodingtx.fa

# exactAGAPnotCodan_pep=${AGAPnotCodAn_dir}/exactclass/AGAPnotCodAn.exactclass.tconsAGAP-precise.refCDS_to_add.ids.pep

# AGAPnotCodan_transdecoder_pep=${AGAPnotCodAn_dir}/AGAPnotCodAn_transcripts.comprehensive.v6.gffread.fa.transdecoder.filt.pep

# notAGAPnotCodan_transdecoder_pep=${notAGAPnotCodAn_dir}/notAGAPnotCodAn_transcripts.comprehensive.v6.gffread.fa.transdecoder.pep

# cat ${exactAGAPnotCodan_pep} ${AGAPnotCodan_transdecoder_pep} ${notAGAPnotCodan_transdecoder_pep} > ${notCodAn_dir}/notCodan_v7_addtl_proteincodingtx.pep

grep -f ${notCodAn_dir}/notCodan_v7_addtl_proteincodingtx.tcons.txt ${notCodAn_dir}/ref_guided_idcollapse.strandfiltered.notCodAn.gtf | grep -w "transcript" | cut -f 9 | tr -d "\";" | sed s'/transcript_id//'g | sed 's/gene_id//g' | sed 's/gene_name//g' > ${notCodAn_dir}/notCodan_v7_addtl_proteincodingtx.tcons-XLOC-AGAP.txt