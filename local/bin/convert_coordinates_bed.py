#!/usr/bin/env python
import numpy as np
from collections import defaultdict 
import sqlite3
from sys import stdin, stderr, argv, stdout
from optparse import OptionParser
import genomeToTranscriptMapper as gttm
import subprocess
import tempfile
import csv
import pyBigWig


def genomic_intervalsToTranscript(exons_db_cursor, bb, transcript_id=None, gene_name=None):
	query='SELECT chr, b, e, strand, gene_name, transcript_id, gene_name FROM exons WHERE transcript_id=?'
	key=transcript_id
	repeats = []
	exons = []
	strand = None

	for exon_row in exons_db_cursor.execute(query, (key,)):
		if (strand is None):
			strand = exon_row[3]

		repeats_intersected = bb.entries(exon_row[0], exon_row[1],  exon_row[2])
		if (repeats_intersected is not None):
			repeats.append(repeats_intersected)
			exons.append((exon_row[1],  exon_row[2]))
	
	if (len(exons)==0):
		return []

	converter = gttm.GenomeToTranscriptMapper(exons, strand)

	for repeats_block in repeats:
		for line in repeats_block:
			b = int(line[0])
			e = int(line[1])
			transcript_coords = converter.convert_interval_genome_to_transcript(b,e)
			yield((transcript_coords[0], transcript_coords[1], line[2]))


def save(data):
	write = csv.writer(stdout, delimiter = '\t')
	write.writerows(data)

def load(exons_db, to_convert_bed, transcript_id=None):
	conn = sqlite3.connect(exons_db)
	cur = conn.cursor()

	if (transcript_id is None):
		sqlite3_cursor = cur.execute("SELECT DISTINCT transcript_id FROM exons WHERE transcript_type = 'lncRNA'")
		transcript_ids = sqlite3_cursor.fetchall()
	else:
		transcript_ids = [[transcript_id]]

	bb = pyBigWig.open(to_convert_bed)
	to_exp = []

	for index, transcript_id in enumerate(transcript_ids):
		transcript_id = transcript_id[0]
		if (index % 100 == 0):
			print(f"{index} out of {len(transcript_ids)}", file=stderr)
		for line in genomic_intervalsToTranscript(cur, bb, transcript_id=transcript_id):
			to_exp.append((line[0], line[1], transcript_id, line[2]))
	save(to_exp)
	

def main():
	usage = "%prog EXONS_DB.sqlite TO_CONVERT.bb > CONVERTED.bed"
	parser = OptionParser(usage=usage)
	
	parser.add_option('-i', '--id', type=str, dest='transcript_id', default=None, help='Return value only for this transcript [default: %default]', metavar='TRANSCRIPT_ID')
	options, args = parser.parse_args()
	
	if len(args) != 2:
		exit('Unexpected argument number.')
	
	EXONS_DB = argv[1]; REPEATS_BED = argv[2]
	load(args[0], args[1], options.transcript_id)



if __name__ == '__main__': main()
