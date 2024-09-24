import os
import pandas as pd
from stateAnalysis_tools import calculate_area_from_residues_with_chain

#Example usage of the intracellular gate area calculation
if __name__ == "__main__":
    pdb_file = '/Users/adrianahernandezgonzalez/LabNotebook/Sandeep/AF_vanilla/256-512/partial_alphaCav12_unrelaxed_rank_002_alphafold2_ptm_model_3_seed_081.pdb'
    # S6_I: LEU290, SER 294, S6_II: LEU 638, V642, S6_III: V1071, ILE1075, S6_IV: I1401, F1405

    residues_IG = [('LEU', 290), ('SER', 294), ('LEU', 638), ('VAL', 642), ('VAL',1071),('ILE',1075),('ILE',1401),('PHE',1405)]  # Example residues
    chain = 'A'

    #print("\nArea between residues:")
    residue_pairs_distances, area = calculate_area_from_residues_with_chain(pdb_file, residues_IG, chain)

    print("\nResidue pairs and distances:")
    for pair, distance in residue_pairs_distances.items():
        print(f"{pair}: {distance:.2f} Å")

    print(f"\nTotal area: {area:.2f} square Å")

    # Selectivity Filter: E252,E595, E1024, E1353
    residues_SF= [('GLU', 252), ('GLU', 595), ('GLU', 1024), ('GLU', 1353)]  # Example residues
    chain_SF = 'A'
    #print("\nArea between residues:")
    residue_pairs_distances, area = calculate_area_from_residues_with_chain(pdb_file, residues_SF, chain)

    print("\nResidue pairs and distances:")
    for pair, distance in residue_pairs_distances.items():
        print(f"{pair}: {distance:.2f} Å")

    print(f"\nTotal area: {area:.2f} square Å")

    #up
    #VSDI_S3: ASP108 VSDII_S2: ASP439 VSDIII_S2:PHE817, VSDIV_S3 : LEU1155
    residues_Up= [('ASP', 108), ('ASN', 439), ('PHE', 817), ('LEU', 1155)]  # Example residues
    chain_SF = 'A'
    #print("\nArea between residues:")
    residue_pairs_distances, area = calculate_area_from_residues_with_chain(pdb_file, residues_Up, chain)

    print("\nResidue pairs and distances:")
    for pair, distance in residue_pairs_distances.items():
        print(f"{pair}: {distance:.2f} Å")

    print(f"\nTotal area: {area:.2f} square Å")

    #Down
    #VSDI_S2: ILE68  VSDII_S2: TYR465 VSDIII_S2: ILE841, VSDIV_S2: TYR1129
    residues_Down= [('ILE',68), ('TYR', 465), ('ILE', 841), ('TYR', 1129)]  # Example residues
    chain_SF = 'A'
    #print("\nArea between residues:")
    residue_pairs_distances, area = calculate_area_from_residues_with_chain(pdb_file, residues_Down, chain)

    print("\nResidue pairs and distances:")
    for pair, distance in residue_pairs_distances.items():
        print(f"{pair}: {distance:.2f} Å")

    print(f"\nTotal area: {area:.2f} square Å")




