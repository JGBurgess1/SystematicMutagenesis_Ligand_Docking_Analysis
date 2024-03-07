import pandas as pd 
import os
import sqlite3
import matplotlib as plt

# connect to db
con = sqlite3.connect('qcrB_mutations.db')
cursor = con.cursor()

df = pd.read_sql_table('scores', con)

print(df.head())

