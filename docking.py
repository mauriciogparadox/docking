from vina import Vina
import MDAnalysis as mda
import numpy as np
import time

# Function to convert pdb to pdbqt

def convert_pdb_to_pdbqt(input_pdb, output_pdbqt):
    # Command to convert PDB to PDBQT using OpenBabel
    command = f'obabel -ipdb {input_pdb} -opdbqt -O {output_pdbqt}'

    # Execute the command using subprocess
    try:
        subprocess.run(command, shell=True, check=True)
        print(f'Conversion successful: {input_pdb} -> {output_pdbqt}')
    except subprocess.CalledProcessError as e:
        print(f'Error during conversion: {e}')

# Functions to create box size

def compute_bounding_box(pdb_file):
    # Load the protein structure
    u = mda.Universe(pdb_file)
    
    # Get the positions of all atoms in the protein
    positions = u.atoms.positions
    
    # Compute the minimum and maximum coordinates
    min_coords = np.min(positions, axis=0)
    max_coords = np.max(positions, axis=0)
    
    # Compute the center and size of the bounding box
    center = (min_coords + max_coords) / 2
    size = max_coords - min_coords

def adjust_box_size(center, size, padding=5.0):
    # Add padding to the size to ensure the ligand can fully explore the binding site
    adjusted_size = size + padding
    return center, adjusted_size

# Function to calculate center of mass

def calculate_center_of_mass(pdb_file):
    atoms = [atom for atom in structure.get_atoms()]
    center = [sum(coord) / len(atoms) for coord in zip(*[atom.coord for atom in atoms])]
    return center

# Provide files

receptor_file = input("Provide a receptor file name:")
receptor_pdbqt = f"{receptor_pdb}qt"
ligand_file = input("Provide a ligand file name:")
ligand_pdbqt = f"{ligand_file}qt"

convert_pdb_to_pdbqt(receptor_file, receptor_pdbqt)
convert_pdb_to_pdbqt(ligand_file, ligand_pdbqt)

center, size = compute_bounding_box(receptor_file)

center, adjusted_size = adjust_box_size(center, size, padding=5.0)

print(f"Box Size: {adjusted_size}")

center_receptor = calculate_center_of_mass(receptor_file)

print('Center of receptor:', center_receptor)

# Start vina computing

v = Vina(sf_name='vina', cpu = -1)

v.set_receptor(receptor_pdbqt)

v.compute_vina_maps(center=center_receptor, box_size=adjusted_size)

v.set_ligand_from_file(ligand_pdbqt)

start_time = time.time()

v.dock(exhaustiveness=32, n_poses=20)

end_time = time.time()

elapsed_time = end_time - start_time

print(f"Time taken: {elapsed_time:.2f} seconds")

pose = input("Provide a pose from the list:")

v.write_poses('example_output.pdbqt', n_poses=pose, overwrite=True)

# END
