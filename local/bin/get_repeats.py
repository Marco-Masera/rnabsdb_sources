import pandas as pd
import numpy as np
from collections import defaultdict 
import sqlite3
from sys import stdin, stderr, argv
from optparse import OptionParser
from cgat import GTF
import transcriptCoordInterconverter as tci
import sqlite3
import subprocess
import tempfile
import csv

EXONS_DB=''
REPEATS_BED=''

def get_converter(sqlite_cursor, transcript_id, gene_name):
	key=None
	if transcript_id is not None:
		query='SELECT chr, b, e, strand, gene_name, transcript_id, gene_name FROM exons WHERE transcript_id=?'
		key=transcript_id
	else:
		query='SELECT chr, b, e, strand, gene_name, transcript_id, gene_name FROM exons WHERE gene_name=?'
		key=gene_name
	if key is None:
		raise ValueError("transcript_id or gene_name should be not None")
		
	empty=True
	fp = tempfile.NamedTemporaryFile(mode="w+t")
	for row in sqlite_cursor.execute(query, (key,)):
		fp.write('%s\tOLILAB\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; exon_number "0"; exon_id "NA"; transcript_type "NA"; gene_name "%s"; transcript_name "NA";\n' % row)
		empty=False
	fp.seek(0)
	if empty:
		raise Exception("Transcript not found (%s)" % transcript_id)

	transcript = list(GTF.transcript_iterator(GTF.iterator(fp)))[0]
	return (tci.TranscriptCoordInterconverter(transcript), fp)



def genomic_intervalsToTranscript(exons_db, intersect_bed=None, transcript_id=None, gene_name=None):
	conn = sqlite3.connect(exons_db)
	cur = conn.cursor()

	converters={}
	out = []
	
	bedtools_proc = None

	if intersect_bed is not None:
		converter, exons_gtf_fp = get_converter(cur,transcript_id, gene_name)
		converters[transcript_id]=converter
		cmd='bedtools intersect -a {exons} -b {intersect_bed} -sorted -wb'.format(exons=exons_gtf_fp.name, intersect_bed=intersect_bed)
		bedtools_proc = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="ascii")
		input_file=iter(bedtools_proc.stdout.readline,'')
	else:
		input_file=stdin

	for line in input_file:
		line = line.strip().split("\t")
		if  intersect_bed:
			chrom = line[0]
			b = int(line[3])-1 #the input of bedtools intersect is a gtf
			e = int(line[4])
		else:
			chrom = line[0]
			b = int(line[1])
			e = int(line[2])
			transcript_id = line[3]

		converter = converters.get(transcript_id,None)
		if converter is None:
			converter, exons_gtf_fp =get_converter(cur,transcript_id)
			exons_gtf_fp.close()
			converters[transcript_id]=converter

		try:    #TODO workaround #178638453
			transcript_coords = converter.genome_interval2transcript((b,e))
		except ValueError:
			transcript_coords = converter.genome_interval2transcript((b,e-1))

		yield((line,transcript_coords[0],transcript_coords[1]))

	bedtools_proc.communicate()
	if intersect_bed and bedtools_proc.returncode!=0:
		raise Exception("bedtools exit with {code}, error:\n{e} ".format(code=bedtools_proc.returncode, e=bedtools_proc.stderr.read()))


def get_repeat(transcript_id, exons_db, repeats_bed):
	TEclass=defaultdict(list)
	conn = sqlite3.connect(exons_db)
	for line,b,e in genomic_intervalsToTranscript(exons_db, repeats_bed, transcript_id=transcript_id):
		te_class=line[-1]
		TEclass[te_class].append((b,e))
	return TEclass

def save(data, export_path):
	with open('export_path', 'a') as f:
		write = csv.writer(f)
		write.writerows(data)

def load(export_path):
	conn = sqlite3.connect(EXONS_DB)
	cur = conn.cursor()
	sqlite3_cursor = cur.execute('SELECT DISTINCT transcript_id FROM exons WHERE true')
	transcript_ids = sqlite3_cursor.fetchall()
	i=0
	for transcript_id in transcript_ids:
		i += 1
		print(f"Iteration {i} over {len(transcript_ids)}")
		to_exp = []
		transcript_id = transcript_id[0]
		result = get_repeat(transcript_id, EXONS_DB, REPEATS_BED)
		for repClass in result.keys():
			for start, end in result[repClass]:
				to_exp.append([transcript_id, repClass, start, end])
		save(to_exp, export_path)
		sql = f"DELETE from exons WHERE transcript_id = '{transcript_id}'"
		cur.execute(sql)
		conn.commit()
	
  
	

def main():
	global EXONS_DB; global REPEATS_BED
	if (len(argv)<4):
		exit("Please provide path to exons_db, repeats_bed and exported csv")
	EXONS_DB = argv[1]; REPEATS_BED = argv[2]
	load(argv[3])

if __name__ == '__main__': main()
