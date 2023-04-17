#!/usr/bin/env python
import pandas as pd
import numpy as np
from collections import defaultdict 
import sqlite3
from sys import stdin, stderr, argv
from optparse import OptionParser
import genomeToTranscriptMapper as gttm
import subprocess
import tempfile
import csv
import pyBigWig

class Transcript:
	def __init__(self,b, e, chr, gene_id, transcript_id, strand):
		self.b = b; self.e = e; self.chr = chr; self.gene_id = gene_id 
		self.transcript_id = transcript_id; self.strand = strand

def genomic_intervalsToTranscript(exons_db_cursor, bb, transcript_id=None, gene_name=None):
	query='SELECT chr, b, e, strand, gene_name, transcript_id, gene_name FROM exons WHERE transcript_id=?'
	key=transcript_id
	segments = []
	transcripts = []
	for row in exons_db_cursor.execute(query, (key,)):
		intersected = bb.entries(row[0], row[1],  row[2])
		transcript = Transcript(row[1],  row[2], row[0], row[4], row[5], row[3])
		if (intersected is not None):
			intersected = [ [line[0], line[1], transcript] for line in intersected]
			segments.append(intersected)
		transcripts.append(transcript)
	
	if (len(segments)==0):
		return []

	converter = gttm.GenomeToTranscriptMapper([(t.b, t.e) for t in transcripts], transcripts[0].strand)

	for line_enxon in segments:
		for line in line_enxon:
			
			b = int(line[0])-1
			e = int(line[1])
			try:
				transcript_coords = converter.convert_interval_genome_to_transcript(b,e)
			except ValueError:
				#transcript_coords = converter.convert_interval_genome_to_transcript(b,e-1)
				print(f"Interval {b}-{e} not in exons", file=stderr)
				transcript_coords = []

			yield((line,transcript_coords[0],transcript_coords[1]))


def save(data, export_path):
	with open(export_path, 'a') as f:
		write = csv.writer(f)
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

	for transcript_id in transcript_ids:
		transcript_id = transcript_id[0]
		for line,b,e in genomic_intervalsToTranscript(cur, bb, transcript_id=transcript_id):
			#line="\t".join(line)
			#print(f"{transcript_id}\t{int(b)}\t{int(e)}\t{line}")
			print(line)
	

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
