import MySQLdb
import csv
from sys import argv

def load_into_db():
    path = argv[1]
    tab_name = 'lncrna_repeat'
    db = MySQLdb.connect("localhost","rnabs","password","rnabsdb" )
    cursor = db.cursor()
    skipped = 0
    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        last_id = None; last_transcript = None
        d = dict()
        for row in csv_reader:
            d[str(row)] = row

        for key in d:
            row = d[key]
            #Line: [transcript_id, repClass, start, end]
            if (last_id != None and last_transcript == row[0]):
                id_ = last_id
            else:
                transcript_ref_query = f"SELECT id FROM lncrna_transcript WHERE transcript_id = '{row[0]}';"
                cursor.execute(transcript_ref_query)
                fetched = cursor.fetchall()
                if (len(fetched)==0):
                    skipped += 1
                    continue 
                id_ = fetched[0][0]
                last_id = id_; last_transcript = row[0]
            sql = f"INSERT INTO {tab_name}(start, end, repClass, transcript_id) values ({row[2]},{row[3]},'{row[1]}',{id_});"
            cursor.execute(sql)
    
    db.commit()
    db.close()
    print(f"Skipped {skipped} lines")

if __name__== "__main__": load_into_db()