#!/usr/bin/env python
import numpy as np
from collections import defaultdict 
import json
import sqlite3
from sys import stdin, stderr, argv, stdout
from optparse import OptionParser
import genomeToTranscriptMapper as gttm
import subprocess
import tempfile
import csv
import pyBigWig
import json

def get_signal_for_transcript(exons_cursor, full_signal, transcript_id):
	query='SELECT chr, b, e, strand FROM exons WHERE transcript_id=? ORDER BY b;'
	strand = None
	signal = []
	exons = exons_cursor.execute(query, (transcript_id,)).fetchall()
	if (len(exons)==0):
		return ""
	strand = exons[0][3]
	for exon_row in exons:
		signal.extend(full_signal.values(exon_row[0], exon_row[1], exon_row[2]))

	if(strand=="-"):
		signal.reverse()

	return signal


def convert(exons_db, to_convert_bed, transcript_id = None):
	conn = sqlite3.connect(exons_db)
	cur = conn.cursor()
	if (transcript_id is None):
		sqlite3_cursor = cur.execute("SELECT DISTINCT transcript_id FROM exons WHERE transcript_type = 'lncRNA'")
		transcript_ids = sqlite3_cursor.fetchall()
	else:
		transcript_ids = [[transcript_id]]

	signal = pyBigWig.open(to_convert_bed)

	for index, transcript_id in enumerate(transcript_ids):
		if (index % 100 == 0):
			print(f"{index} out of {len(transcript_ids)}", file=stderr)
		transcript_id = transcript_id[0]
		converted = get_signal_for_transcript(cur, signal, transcript_id)
		print(f">{transcript_id}")
		print(json.dumps(converted))

def main():
	usage = "%prog EXONS_DB.sqlite SIGNAL.bw  > CONVERTED.fasta"
	parser = OptionParser(usage=usage)
	parser.add_option('-i', '--id', type=str, dest='transcript_id', default=None, help='Return value only for this transcript [default: %default]', metavar='TRANSCRIPT_ID')
	options, args = parser.parse_args()
	
	if len(args) != 2:
		exit('Unexpected argument number.')
	
	EXONS_DB = argv[1]; REPEATS_BED = argv[2]
	convert(args[0], args[1], options.transcript_id)



if __name__ == '__main__': main()
