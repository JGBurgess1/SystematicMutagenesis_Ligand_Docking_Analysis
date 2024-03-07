import sqlite3
import os
import csv
import pandas as pd


# Create a database for all results returned.
# Export as necessary, using SQL filtering.

db_path = 'qcrB_mutations.db'
con = sqlite3.connect(db_path)
cursor = con.cursor()

df_new = pd.read_csv('membrane_mCSM_pred.txt')

for index, row in df_new.iterrows():
    update_query = f'''
    UPDATE scores
    SET mCSM_membrane_stability = '{row[3]}' -- some comment
    WHERE mutation = '{row[1]}' -- matching the mutation from the df to the table
    '''
    cursor.execute(update_query)
    con.commit()

con.close()
