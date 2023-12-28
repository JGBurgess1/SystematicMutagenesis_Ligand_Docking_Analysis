To use the programs, you need to have loaded a schrodinger maestro environment.

Run programs in order (1_Mu.. 2_Grid.. 3_Dat..)

This code will systematically mutate every possible residue to every possible amino acid, 
locally energy minimize the structure,
then dock a given ligand file to the mutated protein.

The docking scores are exported to a SQLITE3 database, and can be accessed using SQL queries as required.

