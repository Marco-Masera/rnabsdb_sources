#!/usr/bin/env python

import MySQLdb
import csv
from optparse import OptionParser
from sys import argv

DB_TAB_NAME = 'lncrna_transcript'

def load_data_to_db(transcripts_bed_path, transcript_id_to_gene_id_path):
    #Load transcripts into memory data structure
    transcripts = {}
    with open(transcripts_bed_path, "rt", encoding='ascii') as transcripts_bed:
        read = csv.reader(transcripts_bed, delimiter='\t')
        for row in read: 
            if (row[6]== 'lncRNA'):
                if (row[3] in transcripts):
                    transcripts[row[3]][1] = min(transcripts[row[3]][1], row[1])
                    transcripts[row[3]][2] = max(transcripts[row[3]][2], row[2])
                else:
                    transcripts[row[3]] = row
    #Read and parse transcript_id to gene_id mapping data
    transcript_to_gene = {}
    with open(transcript_id_to_gene_id_path, "rt", encoding='ascii') as mapping_file:
        read = csv.reader(mapping_file, delimiter='\t')
        for row in read:
            transcript_to_gene[row[1]] = row[0]

    # Open database connection
    db = MySQLdb.connect("localhost","rnabs","password","rnabsdb" )
    # prepare a cursor object using cursor() method
    cursor = db.cursor()
    i = 0;
    for key in transcripts.keys():
        row = transcripts[key]  
        sql = f"INSERT INTO {DB_TAB_NAME}(chromosome,start,end,transcript_id,score,strand,transcript_type,gene_name, gene_id) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s);"   
        cursor.execute(sql, (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], transcript_to_gene[row[3]]))

    # Commit your changes in the database
    db.commit()
    # disconnect from server
    db.close()

def main():
    usage = "%prog TRANSCRIPT_LIST.bed TRANSCRIPT_ID_TO_GENE_ID.map"
    parser = OptionParser(usage=usage)
    options, args = parser.parse_args()
    if len(args) != 2:
        exit('Unexpected argument number.')
	
    load_data_to_db(args[0], args[1])


if __name__ == '__main__': main()
