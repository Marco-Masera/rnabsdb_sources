#!/usr/bin/env python

import MySQLdb
import csv
import pyfastx
from optparse import OptionParser
from sys import argv

DB_TAB_NAME = 'lncrna_fulltranscript'

def load_data_to_db(transcripts_fasta_path):
    # Open database connection
    db = MySQLdb.connect("localhost","rnabs","password","rnabsdb" )
    # prepare a cursor object using cursor() method
    cursor = db.cursor()

    #Transcript_ids is a dict with Key: transcript_id => Value: database_pk_for_transcript
    transcript_ids = {}
    sql = "SELECT id, transcript_id FROM lncrna_transcript WHERE true;"
    cursor.execute(sql)
    fetched = cursor.fetchall()
    for row in fetched:
        transcript_ids[row[1]] = row[0]

    fa = pyfastx.Fastx(transcripts_fasta_path)
    for name, seq in fa:
        id_ = name.split("|")[0]
        if (id_ in transcript_ids):
            sql = f"INSERT INTO {DB_TAB_NAME}(full_transcript, transcript_id) VALUES (%s, %s);"   
            cursor.execute(sql, (seq, transcript_ids[id_]))

    # Commit your changes in the database
    db.commit()
    # disconnect from server
    db.close()

def main():
    usage = "%prog TRANSCRIPTS.fa.gz"
    parser = OptionParser(usage=usage)
    options, args = parser.parse_args()
    if len(args) != 1:
        exit('Unexpected argument number.')
	
    load_data_to_db(args[0])


if __name__ == '__main__': main()
