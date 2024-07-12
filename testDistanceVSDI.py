from measure_distance import measure_distances, measure_shortest_distance

pdb_file = './Structures/8wea.pdb'
residue1_type = 'ARG'
residue1_number = 240
residue2_type = 'GLU'
residue2_number = 163


print("All distances:")
measure_distances(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number)

print("\nShortest distance:")
measure_shortest_distance(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number)


print("All distances:")
measure_distances(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number)

print("\nShortest distance:")
measure_shortest_distance(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number)
