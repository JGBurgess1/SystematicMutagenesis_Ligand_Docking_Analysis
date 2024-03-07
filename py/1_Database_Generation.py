import sqlite3
import os
import csv

# Create a database for all results returned.
# Export as necessary, using SQL filtering.

db_path = 'qcrB_mutations.db'
con = sqlite3.connect(db_path)
cursor = con.cursor()

# table already exists? 

# csv folder dir

for filename in os.listdir(os.getcwd()):
    #print(filename)
    if filename.endswith(".csv"):
        with open(os.path.join(os.getcwd(), filename), "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            next(csv_reader)

            for row in csv_reader: #one row in the file...
                unique_name = row[5][:-9] #'B_GLU_19_LEU'
                parts = unique_name.split("_") #['B', 'GLU', '19', 'LEU']
                old_aa = parts[1] #GLU
                resnum = parts[2] #19
                new_aa = parts[3] #LEU
                docking_score = row[4] #-11.8277
                glide_e_score = row[15] #-112.34

                SQL_statement = f'''
                INSERT INTO mutations VALUES (?, ?, ?, ?, ?, ?)
                '''
                cursor.execute(SQL_statement, (unique_name, old_aa, resnum, new_aa, docking_score, glide_e_score))
                con.commit()
con.close()

                




