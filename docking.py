from vina import Vina
import MDAnalysis as mda
import numpy as np 

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
ligand_file = input("Provide a ligand file name:")

center, size = compute_bounding_box(receptor_file)

center, adjusted_size = adjust_box_size(center, size, padding=5.0)

print(f"Box Size: {adjusted_size}")

center_receptor = calculate_center_of_mass(receptor_file)

print('Center of receptor:', center_receptor)

# Start vina computing

v = Vina(sf_name='vina', cpu = -1)

v.set_receptor(receptor_file)

v.compute_vina_maps(center=center_receptor, box_size=adjusted_size)

v.set_ligand_from_file(ligand_file)

v.dock(exhaustiveness=32, n_poses=20)

pose = input("Provide a pose from the list:")

v.write_poses('example_output.pdbqt', n_poses=pose, overwrite=True)

# END