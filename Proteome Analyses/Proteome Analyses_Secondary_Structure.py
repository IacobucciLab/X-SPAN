#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import tarfile
import os
import gzip

def extract_all_proteins(tar_file, output_dir):
    """
    Extract all .pdb.gz files from a tar file and decompress them.

    Args:
        tar_file (str): Path to the tar file containing all the proteins.
        output_dir (str): Directory to save the extracted protein files.
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Open the tar file
    try:
        with tarfile.open(tar_file, 'r') as tar_ref:
            # Extract all .pdb.gz files
            for member in tar_ref.getmembers():
                if member.name.endswith('.pdb.gz'):
                    try:
                        # Extract and decompress .pdb.gz files
                        with tar_ref.extractfile(member) as f_in:
                            with gzip.GzipFile(fileobj=f_in) as gz_in:
                                pdb_filename = os.path.join(output_dir, member.name.replace('.gz', ''))
                                with open(pdb_filename, 'wb') as f_out:
                                    f_out.write(gz_in.read())
                                    print(f"Extracted and saved {pdb_filename}")
                    except Exception as e:
                        print(f"Error processing {member.name}: {e}")
                else:
                    # Ignore other file types
                    print(f"Ignored {member.name}")
    except Exception as e:
        print(f"Error opening tar file {tar_file}: {e}")

# Example usage
tar_file = 'C:/Users/users/Downloads/UP000005640_9606_HUMAN_v4.tar'  # Update with your actual tar file
output_dir = 'C:/Users/users/Downloads/HUMAN_extracted_proteins'  # Output directory to extract proteins into
extract_all_proteins(tar_file, output_dir)


# In[ ]:


#estrazione struttura secondaria da proteoma senza resname 
import os
import pymol
from pymol import cmd

def extract_and_process_pdb_files(input_directory, output_directory, batch_size):
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Initialize PyMOL in non-GUI mode
    pymol.finish_launching(['pymol', '-qc'])  # -q for quiet mode, -c for no GUI

    # Get a list of all PDB files in the directory
    pdb_files = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if f.endswith('.pdb')]

    def process_batch(pdb_files, batch_index):
        # Clear PyMOL session
        cmd.reinitialize()

        # Load and process each PDB file
        for pdb_file in pdb_files:
            cmd.load(pdb_file)

        # Dictionaries to store secondary structure and pLDDT scores
        secondary_structure = {}
        pLDDT_scores = {}

        # Collect secondary structure
        def collect_secondary_structure(model, resi, ss):
            secondary_structure[(model, resi)] = ss

        # Collect pLDDT scores
        def collect_plddt_scores(model, resi, b):
            pLDDT_scores[(model, resi)] = b

        # Define a PyMOL expression to collect data
        cmd.iterate("n. CA", "collect_secondary_structure(model, resi, ss)", space={"collect_secondary_structure": collect_secondary_structure})
        cmd.iterate("n. CA", "collect_plddt_scores(model, resi, b)", space={"collect_plddt_scores": collect_plddt_scores})

        # Group residues by protein and sort them by ID
        residues_by_protein = {}
        for (model, resi) in secondary_structure.keys():
            if model not in residues_by_protein:
                residues_by_protein[model] = []
            residues_by_protein[model].append((resi, secondary_structure[(model, resi)], pLDDT_scores.get((model, resi), "")))

        for model in residues_by_protein:
            residues_by_protein[model].sort(key=lambda x: int(x[0]))  # Sort by residue ID

        # Write secondary structure information and pLDDT scores into a text file
        output_file_path = os.path.join(output_directory, f"output_batch_{batch_index}.txt")
        with open(output_file_path, 'w') as f:
            f.write("PDB_ID\tResidue_ID\tSecondary_Structure\tpLDDT_Score\n")
            for model in sorted(residues_by_protein.keys()):
                for (resi, ss, pLDDT) in residues_by_protein[model]:
                    pdb_id = "{}_{}".format(model, resi)  # Concatenate model name and residue ID
                    f.write("{}\t{}\t{}\t{:.2f}\n".format(pdb_id, resi, ss, float(pLDDT)))

        print(f"Information saved to {output_file_path}")

    # Process the PDB files in batches
    for i in range(0, len(pdb_files), batch_size):
        batch_files = pdb_files[i:i + batch_size]
        process_batch(batch_files, i // batch_size + 1)

    # Quit PyMOL
    pymol.cmd.quit()

# Example usage
input_directory = 'C:/Users/users/Downloads/HUMAN_extracted_proteins'   # Update with your actual input directory
output_directory = 'C:/Users/users/Downloads/HUMAN_extracted_proteins_output'   # Update with your actual output directory
batch_size = 1000  # Number of PDB files to process per batch

extract_and_process_pdb_files(input_directory, output_directory, batch_size)


# In[ ]:


# Frequency analysis script for residue relationships from multiple .txt files
import pandas as pd
import os

# Input and output folders
input_folder = 'C:/Users/users/Desktop/HUMAN_extracted_proteins_output/'
output_folder = 'C:/Users/users/Desktop/AA_rel_HUMAN/'

# Make sure output folder exists
os.makedirs(output_folder, exist_ok=True)

# List all .txt files in the input folder
txt_files = [f for f in os.listdir(input_folder) if f.endswith('.txt')]

# List of the 20 standard amino acids
amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
               'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# Function to calculate residue pair distances for entire protein
def calculate_relationship_whole(df, residue_name, window_size=20):
    relationships = {f'1-{i}': 0 for i in range(2, window_size+1)}
    for protein_id, protein_df in df.groupby('PDB_ID'):
        residue_positions = protein_df[protein_df['Residue_Name'] == residue_name]['Residue_ID'].tolist()
        for i, pos1 in enumerate(residue_positions):
            for pos2 in residue_positions[i+1:]:
                dist = abs(pos2 - pos1) + 1
                if 2 <= dist <= window_size:
                    relationships[f'1-{dist}'] += 1
    return relationships

# Function to extract continuous segments of specific secondary structure
def find_continuous_segments(df, secondary_structure):
    segments = []
    current_segment = []
    for _, row in df.iterrows():
        if row['Secondary_Structure'] == secondary_structure:
            current_segment.append(row)
        else:
            if current_segment:
                segments.append(current_segment)
            current_segment = []
    if current_segment:
        segments.append(current_segment)
    return segments

# Function to calculate residue pair distances within continuous structural motifs
def calculate_relationship_continuous_segments(segments, residue_name, window_size=20):
    relationships = {f'1-{i}': 0 for i in range(2, window_size+1)}
    for segment in segments:
        residue_positions = [row['Residue_ID'] for row in segment if row['Residue_Name'] == residue_name]
        for i, pos1 in enumerate(residue_positions):
            for pos2 in residue_positions[i+1:]:
                dist = abs(pos2 - pos1) + 1
                if 2 <= dist <= window_size:
                    relationships[f'1-{dist}'] += 1
    return relationships

# Analyze all .txt files in the folder
for idx, file_name in enumerate(txt_files, start=1):
    print(f'[{idx}/{len(txt_files)}] Processing file: {file_name}')
    file_path = os.path.join(input_folder, file_name)
    df = pd.read_csv(file_path, delimiter='\t')

    # Extract base protein ID and convert Residue_ID to numeric
    df['PDB_ID'] = df['PDB_ID'].str.split('_').str[0]
    df['Residue_ID'] = pd.to_numeric(df['Residue_ID'], errors='coerce')

    output_file = os.path.join(output_folder, f'Amino_Acid_Relationships_output_{idx}_HUMAN.xlsx')

    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        for aa in amino_acids:
            print(f'  -> Analyzing amino acid: {aa}')

            # Whole protein distance calculation
            relationships_whole = calculate_relationship_whole(df, aa)

            # Distance calculation within motifs H, S, L
            motif_relationships = {}
            for motif_type in ['H', 'S', 'L']:
                segments = find_continuous_segments(df, motif_type)
                motif_relationships[motif_type] = calculate_relationship_continuous_segments(segments, aa)

            # Build dataframe with all results
            relationship_df = pd.DataFrame({
                'Amino Acid Relationship': [f'1-{i}' for i in range(2, 21)],
                'Whole_Protein': [relationships_whole.get(f'1-{i}', 0) for i in range(2, 21)]
            })

            for motif_type, relationships in motif_relationships.items():
                motif_values = [relationships.get(f'1-{i}', 0) for i in range(2, 21)]
                relationship_df[f'Continuous_{motif_type}'] = motif_values

            # Save sheet named by amino acid
            relationship_df.to_excel(writer, sheet_name=aa, index=False)

print('Analysis completed.')


# In[ ]:




