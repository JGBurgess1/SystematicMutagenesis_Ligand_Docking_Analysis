import sqlite3
import pandas as pd
import matplotlib.pyplot as plt

db_path = 'qcrB_mutations.db'
con = sqlite3.connect(db_path)
cursor = con.cursor()

SQL_query = 'SELECT * from scores ORDER BY sig_q203 ASC'
# print out the values from the table.
df = pd.read_sql(SQL_query, con)
print(df.head(10).to_string(index = False))
print(df.tail(10).to_string(index = False))

con.close()

plt.scatter(df["res_Num"], df["q203_glide"], label = "Q203_Docking_Score")

#plt.scatter(df['ID'], df['SP_Docking_Score'], label="Docking Score")
#plt.axis((-11,-12,14,30))
plt.savefig('plot.png')

#plt.show()
