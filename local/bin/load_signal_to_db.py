#!/usr/bin/env python
import MySQLdb
import csv
from optparse import OptionParser
from Bio import SeqIO
from sys import argv

def int_(number):
    return int(number.split(".")[0])

def load_data_to_db(file_path, tab_name):
    db = MySQLdb.connect("localhost","rnabs","password","rnabsdb" )
    cursor = db.cursor()
    s = set()
    for record in SeqIO.parse("phastCons30way.by.transcript_id.fasta","fasta"):
        if (record.id in s):
            continue 
        s.add(record.id)
        transcript_ref_query = f"SELECT id FROM lncrna_transcript WHERE transcript_id = '{record.id}';"
        cursor.execute(transcript_ref_query)
        fetched = cursor.fetchall()
        if (len(fetched)==0):
            continue
        pk = fetched[0][0]
        sql = f"INSERT INTO {tab_name} (transcript_id, signal_formatted) VALUES (%s, %s);"
        cursor.execute(sql, (pk, record.seq))
        
    db.commit()
    db.close()

def main():
    usage = "%prog SIGNAL.fasta TABLE_NAME"
    parser = OptionParser(usage=usage)
    options, args = parser.parse_args()
    if len(args) != 2:
        exit('Unexpected argument number.')
	
    load_data_to_db(args[0], args[1])


if __name__ == '__main__': main()
