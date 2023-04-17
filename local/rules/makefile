-include $(PRJ_ROOT)/local/rules/_header.mk
GENCODE_REFENRECE_DATA=$(BIOINFO_REFERENCE_ROOT)/gencode/dataset/$(GENCODE_SPECIES)/$(GENCODE_VERSION)

CONDA_ROOT=/opt/conda
CONDA_VERSION=miniconda3
CONDA_ACTIVATE=set +u; source $(CONDA_ROOT)/$(CONDA_VERSION)/etc/profile.d/conda.sh; conda activate ; conda activate 


#Basic transcripts data:
primary_assembly.annotation.ensg2enst.map: $(GENCODE_REFENRECE_DATA)/primary_assembly.annotation.ensg2enst.map
	ln -s $< $@
basic.annotation.exon.longest_transcript.bed: $(GENCODE_REFENRECE_DATA)/basic.annotation.exon.longest_transcript.bed
	ln -s $< $@

load_transcripts: basic.annotation.exon.longest_transcript.bed primary_assembly.annotation.ensg2enst.map
	load_transcripts_to_db.py basic.annotation.exon.longest_transcript.bed primary_assembly.annotation.ensg2enst.map

#Full transcript data:
transcripts.fa.gz: $(GENCODE_REFENRECE_DATA)/transcripts.fa.gz
	ln -s $< $@

load_full_transcripts: transcripts.fa.gz
	load_full_transcripts_to_db.py transcripts.fa.gz

#Conservation: 
#"Raw" bw with conservation
phastCons30way.bw: $(GENCODE_REFENRECE_DATA)/phastCons30way.bw
	ln -s $< $@
#Conservation with coordinates converted for transcripts
#TODO
h42.conservation.coordinates.converted: hg42.phastCons30way.bw
	echo TODO

#Repeats:
#Input files for coordinate converter
basic.annotation.exon.longest_transcript.db: $(GENCODE_REFENRECE_DATA)/basic.annotation.exon.longest_transcript.db
	cp $< $@
# Bed file containing the repeats
repeat_rmsk.clean_chr.bed.gz: $(GENCODE_REFENRECE_DATA)/repeat_rmsk.clean_chr.bed.gz
	cp $< $@
repeat_rmsk.clean_chr.sorted.bed: repeat_rmsk.clean_chr.bed.gz
	zcat repeat_rmsk.clean_chr.bed.gz | bedtools sort > repeat_rmsk.clean_chr.sorted.bed

repeat_rmsk.clean_chr.sorted.noscore.bed: repeat_rmsk.clean_chr.sorted.bed	
	awk '{OFS="\t"; $$5=0; print}' repeat_rmsk.clean_chr.sorted.bed > repeat_rmsk.clean_chr.sorted.noscore.bed

repeat_rmsk.clean_chr.sorted.noscore.bb: repeat_rmsk.clean_chr.sorted.noscore.bed $(GENCODE_REFENRECE_DATA)/ucsc_hg38.chrom.sizes
	bedToBigBed repeat_rmsk.clean_chr.sorted.noscore.bed $(GENCODE_REFENRECE_DATA)/ucsc_hg38.chrom.sizes repeat_rmsk.clean_chr.sorted.noscore.bb


#Repeats with coordinates converted
repeats_coordinate_converted.csv: basic.annotation.exon.longest_transcript.db repeat_rmsk.clean_chr.sorted.noscore.bb
	$(CONDA_ACTIVATE) /home/mmasera/.conda/envs/repeats_env;\
	convert_coordinates_bed.py basic.annotation.exon.longest_transcript.db repeat_rmsk.clean_chr.sorted.noscore.bb > repeats_coordinate_converted.csv



-include $(PRJ_ROOT)/local/rules/_footer.mk
