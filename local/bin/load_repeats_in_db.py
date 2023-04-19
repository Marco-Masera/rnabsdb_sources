#!/usr/bin/env python
import MySQLdb
import csv
from optparse import OptionParser
from sys import argv

def int_(number):
    return int(number.split(".")[0])

def load_data_to_db(file_path):
    tab_name = 'lncrna_repeat'
    db = MySQLdb.connect("localhost","rnabs","password","rnabsdb" )
    cursor = db.cursor()
    skipped = 0
    with open(file_path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        #Map with transcript_id as key
        repeats = dict()
        for row in csv_reader:
            transcript_id = row[2]
            if (transcript_id in repeats):
                repeats[transcript_id].append(row)
            else:
                repeats[transcript_id] = [row]

        #Insert into db one transcript by one
        for transcript_id in repeats.keys():
            #Need to get the primary key of the corresponding transcript
            #Line: [transcript_id, repClass, start, end]
            transcript_ref_query = f"SELECT id FROM lncrna_transcript WHERE transcript_id = '{transcript_id}';"
            cursor.execute(transcript_ref_query)
            fetched = cursor.fetchall()
            if (len(fetched)==0):
                continue
            pk = fetched[0][0]
            s = set()
            for repeat in repeats[transcript_id]:
                #Avoid repetitions:
                to_str = str(repeat[:3])
                if (to_str in s):
                    continue
                s.add(to_str)
                repeat_info = repeat[3].split(";")
                repeat_info[2] = repeat_info[2].split("\t")[0]
                sql = f"INSERT INTO {tab_name}(start, end, repClass, repFamily, repName, transcript_id) values (%s, %s, %s, %s, %s, %s);"
                cursor.execute(sql, (int_(repeat[0]), int_(repeat[1]), repeat_info[2], repeat_info[1], repeat_info[0], str(pk)))
        
    db.commit()
    db.close()

def main():
    usage = "%prog REPEATS.csv"
    parser = OptionParser(usage=usage)
    options, args = parser.parse_args()
    if len(args) != 1:
        exit('Unexpected argument number.')
	
    load_data_to_db(args[0])


if __name__ == '__main__': main()
