-include $(PRJ_ROOT)/local/rules/_header.mk

#your rules here

GENCODE_REFENRECE_DATA=$(BIOINFO_REFERENCE_ROOT)/$(GENCODE_SPECIES)/$(GENCODE_VERSION)

#Conservation: 

#"Raw" bw with conservation
phastCons30way.bw: $(GENCODE_REFENRECE_DATA)/phastCons30way.bw
	ln -s $< $@

#Conservation with coordinates converted for transcripts
#TODO
h42.conservation.coordinates.converted: hg42.phastCons30way.bw
	echo TODO

#Repeats:
#Input files for coordinate converter:
# !! File deve essere copiato perché viene modificato dallo script
basic.annotation.exon.longest_transcript.db: /home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/42/basic.annotation.exon.longest_transcript.db
	cp $< $@
# !!  Non viene modificato quindi può essere copiato o linkato liberamente
repeat_rmsk.bed.gz: /home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/42/repeat_rmsk.txt.gz
	cp $< $@

#With coordinates converted
#load_repeats_in_db.py to load the csv to the database
repeats_converted.csv: basic.annotation.exon.longest_transcript.db repeat_rmsk.bed.gz
	python get_repeats.py basic.annotation.exon.longest_transcript.db repeat_rmsk.bed.gz repeats_converted.csv

-include $(PRJ_ROOT)/local/rules/_footer.mk
