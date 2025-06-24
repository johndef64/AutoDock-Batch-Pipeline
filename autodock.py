#%%
# AutoDock Vina docking script using Meeko for ligand preparation

# Requirements: 
# - AutoDock Vina
# - RDKit
# - Meeko for ligand preparation
# - Python 3.10 or higher
# - Linux or MacOS (Windows support may vary)

# Installation instructions:
# $ conda create -n vina python=3.10
# $ conda activate vina
# $ conda config --env --add channels conda-forge
# $ conda install -c conda-forge numpy swig boost-cpp libboost sphinx sphinx_rtd_theme
# $ y

# install latest meeko
# pip install git+https://github.com/forlilab/Meeko.git

##

import os
from vina import Vina
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from meeko_converter import prepare_ligand_with_meeko, prepare_receptor_with_meeko_cli
from get_pdbs import download_alphafold_simple, download_and_convert_to_pdb
# make a system of arguments to pass to the command line interface of meeko
import argparse


# Define ligands and their corresponding drugs

ligands_drugs = {
    1: {"name": "Acetaminophen", "pubchem_id": 60961},
    2: {"name": "Dopamine", "pubchem_id": 681},
    3: {"name": "Aspirin", "pubchem_id": 2244},
    4: {"name": "Ibuprofen", "pubchem_id": 3672},
    5: {"name": "Caffeine", "pubchem_id": 2519},
    6: {"name": "Atenolol", "pubchem_id": 2249},
    7: {"name": "Loratadine", "pubchem_id": 3954},
    8: {"name": "Omeprazole", "pubchem_id": 4594},
    9: {"name": "Simvastatin", "pubchem_id": 54454},
    10: {"name": "Metformin", "pubchem_id": 4091}
}

drugs_receptor = {
      1:  {'drug': 'Acetaminophen', 'target': 'Cyclooxygenase-1 (COX-1)', 'uniprot_id': 'P23219'},
      2:  {'drug': 'Dopamine', 'target': 'Dopamine receptor D1', 'uniprot_id': 'P21728'},
      3:  {'drug': 'Aspirin', 'target': 'Cyclooxygenase-1 (COX-1)', 'uniprot_id': 'P23219'},
      4:  {'drug': 'Ibuprofen', 'target': 'Cyclooxygenase-1 (COX-1)', 'uniprot_id': 'P23219'},
      5:  {'drug': 'Caffeine', 'target': 'Adenosine receptor A2A', 'uniprot_id': 'P29274'},
      6:  {'drug': 'Atenolol', 'target': 'Beta-1 adrenergic receptor', 'uniprot_id': 'P08588'},
      7:  {'drug': 'Loratadine', 'target': 'Histamine H1 receptor', 'uniprot_id': 'P35367'},
      8:  {'drug': 'Omeprazole', 'target': 'H+/K+ ATPase (proton pump)', 'uniprot_id': 'P20616'},
      9:  {'drug': 'Simvastatin', 'target': 'HMG-CoA reductase', 'uniprot_id': 'P04035'},
      10:  {'drug': 'Metformin', 'target': 'AMP-activated protein kinase (AMPK) alpha-1', 'uniprot_id': 'Q13131'}
    
}


uniprot_id = "P14416"  # Example UniProt ID for a receptor (DRD2)
pubchem_id = 681  # Example PubChem ID for a ligand (Dopamine)

pubchem_id = 709778  # Example PubChem ID for a ligand
uniprot_id = "P07477"  # Example UniProt ID for a receptor

pubchem_id = 709778  # Example PubChem ID for a ligand
uniprot_id = "P00760"  # TRY1_BOVIN


double_id = 3
pubchem_id =  ligands_drugs[ double_id ]["pubchem_id"]
uniprot_id = drugs_receptor[ double_id ]['uniprot_id']  



def parse_arguments():
    """Parse command line arguments for docking"""
    # set a default receptor and ligand
    default_receptor = uniprot_id  # Example UniProt ID for a receptor
    default_ligand = pubchem_id  # Example PubChem ID for a ligand
    parser = argparse.ArgumentParser(description="AutoDock Vina docking script")
    parser.add_argument("--receptor", type=str, default=default_receptor,
                        help="Receptor UniProt IDe")
    parser.add_argument("--ligand", type=str, default=default_ligand,
                        help="Ligand PubChem ID")
    return parser.parse_args()

def is_jupyter_notebook():
    import IPython
    return IPython.get_ipython() is not None
is_notebook = is_jupyter_notebook()

# get the arguments from the command line
# if not ipykernel is running
if not is_notebook:
    args = parse_arguments()
    uniprot_id = args.receptor
    pubchem_id = args.ligand

#%%
class AutoDockVinaRunner:
    def __init__(self, receptor_file, ligand_file):
        self.receptor_file = receptor_file
        self.ligand_file = ligand_file
        self.vina_instance = None
        
    def prepare_receptor(self, receptor_pdbqt):
        """Load and prepare the receptor protein"""
        self.vina_instance = Vina(sf_name='vina')
        self.vina_instance.set_receptor(receptor_pdbqt)
        
    def prepare_ligand(self, ligand_pdbqt):
        """Load and prepare the ligand"""
        # self.vina_instance = Vina(sf_name='vina')
        self.vina_instance.set_ligand_from_file(ligand_pdbqt)
        
    def define_search_space(self, center_x, center_y, center_z, 
                           size_x=20, size_y=20, size_z=20):
        """Define the binding site search space"""
        self.vina_instance.compute_vina_maps(
            center=[center_x, center_y, center_z],
            box_size=[size_x, size_y, size_z]
        )
        
    def run_docking(self, exhaustiveness=8, n_poses=9):
        """Execute the docking simulation"""
        self.vina_instance.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        return self.vina_instance.energies(n_poses=n_poses)
        
    def save_results(self, output_file):
        """Save docking results to file"""
        self.vina_instance.write_poses(output_file, n_poses=9, overwrite=True)




# Ensure directories exist
ligand_folder = "ligand_structures"
receptor_folder = "receptor_structures"
docking_folder = "docking_results"
os.makedirs(ligand_folder, exist_ok=True)
os.makedirs(receptor_folder, exist_ok=True)
os.makedirs(docking_folder, exist_ok=True)

# Download receptor and ligand structures
download_alphafold_simple(uniprot_id, receptor_folder)
download_and_convert_to_pdb(pubchem_id, ligand_folder)

# check if ligand and receptor files exist
ligand_file = os.path.join(ligand_folder, f"ligand_{str(pubchem_id)}.pdb")
receptor_file = os.path.join(receptor_folder, f"AF-{uniprot_id}.pdb")

if not os.path.exists(ligand_file) or not os.path.exists(receptor_file):
    missing_files = True
    raise FileNotFoundError(f"PDB file does not exist.")
else:
    missing_files = False
    print(f"Found ligand file: {ligand_file}")
    print(f"Found receptor file: {receptor_file}")



# Example usage
def main():
    """ Main function to run the docking simulation
    """

    # Prepare ligand and receptor files
    ligand = f"docking.{pubchem_id}.meeko.pdbqt"
    receptor = f"docking.{uniprot_id}.meeko.pdbqt"
    ligand = os.path.join(ligand_folder, ligand)
    receptor = os.path.join(receptor_folder, receptor)

    if not os.path.exists(f"{ligand}"):
        print("Preparing ligand file...")
        prepare_ligand_with_meeko(f"{ligand_file}", ligand)
    if not os.path.exists(f"{receptor}"):
        print("Preparing receptor file...")
        prepare_receptor_with_meeko_cli(f"{receptor_file}", receptor.replace(".pdbqt", ""))


    print(f"Running docking with ligand: {ligand} and receptor: {receptor}")
    # Initialize the docking runner
    docker = AutoDockVinaRunner(receptor, ligand)
    
    # Prepare receptor and ligand
    docker.prepare_receptor(receptor)
    docker.prepare_ligand(ligand)
    
    # Define search space (coordinates of binding site)
    docker.define_search_space(
        center_x=15.0, center_y=10.0, center_z=5.0,
        size_x=20, size_y=20, size_z=20
    )

    # Run docking
    energies = docker.run_docking(exhaustiveness=16, n_poses=5)
    
    # Salva il report usando il print
    report_file = os.path.join(docking_folder, f"docking_report_{uniprot_id}_{pubchem_id}.txt")
    with open(report_file, 'w') as f:
        # Scrivi l'header
        f.write(f"Docking Results for {uniprot_id} - {pubchem_id}\n")
        f.write("=" * 50 + "\n\n")
        
        # Scrivi le energie di binding
        f.write("Binding energies (kcal/mol):\n")
        for i, energy in enumerate(energies):
            f.write(f"Pose {i+1}: {energy[0]:.3f}\n")

    # Save result structures
    docker.save_results(os.path.join(docking_folder, f"docking_results_{uniprot_id}_{pubchem_id}.pdbqt"))
    
    
    print("Docking completed successfully.")
    # print(f"---ENERGIES {energies}---")

    # Print binding energies
    print("Binding energies (kcal/mol):")
    for i, energy in enumerate(energies):
        print(f"Pose {i+1}: {energy[0]:.3f}")

# inclusa la tabella con affinità e RMSD

if __name__ == "__main__":
    if not missing_files:
        print("Ligand and receptor files exist, proceeding with docking...")
        main()
    else:
        print("Ligand or receptor files are missing, please check the download process.")
    

