import sqlite3
import os
import csv
import pandas as pd


# Create a database for all results returned.
# Export as necessary, using SQL filtering.

db_path = 'qcrB_mutations.db'
con = sqlite3.connect(db_path)
cursor = con.cursor()

sql_statement = '''
CREATE TABLE IF NOT EXISTS scores (
    mutation TEXT,
    res_Num INTEGER,
    q203_glide REAL,
    mq9_glide REAL,
    lpzs_glide REAL,
    mCSM_membrane_stability REAL DEFAULT 0.0
);
'''
cursor.execute(sql_statement)
con.commit()

# table already exists? 

# csv folder dir

df_total = {}
aa_dict = {"ALA":'A', "ARG":'R', "ASN":'N', "ASP":'D', "CYS":'C', "GLU":'E', "GLN":'Q', "GLY":'G', "HIS":'H', "ILE":'I', "LEU":'L', "LYS":'K', "MET":'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}


for filename in os.listdir(os.getcwd()):
    if filename.endswith(".csv"):
        print(filename)
        df_new = pd.read_csv(filename)
         
        string = filename.split("_")
        chain = string[0]
        res_num = int(string[2])
        mutation_short_code = aa_dict[string[1]] + string[2] + aa_dict[string[3]] #LEU14ALA => L14A

        q203_glide_score = df_new.iloc[0,4]
        mq9_glide_score = df_new.iloc[1,4]
        lpzs_glide_score = df_new.iloc[2,4]
        
        SQL_statement = f'''
        INSERT INTO scores VALUES (?, ?, ?, ?, ?, ?)
        '''
        cursor.execute(SQL_statement, (mutation_short_code, res_num, q203_glide_score, mq9_glide_score, lpzs_glide_score, 0.0))
        con.commit()       

""" 
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
                INSERT INTO mutations VALUES (?, ?, ?, ?, ?)
                '''
                cursor.execute(SQL_statement, (unique_name, old_aa, resnum, new_aa, docking_score))
                con.commit()
"""

con.close()

                




