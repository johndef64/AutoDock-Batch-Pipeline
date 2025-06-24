#%%
from meeko import MoleculePreparation
from rdkit import Chem
from rdkit.Chem import AllChem

def prepare_ligand_with_meeko(input_file, output_file):
    """Prepare ligand using Meeko"""
    # Read molecule
    if input_file.endswith('.mol'):
        mol = Chem.MolFromMolFile(input_file)
    elif input_file.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(input_file)
    else:
        mol = Chem.MolFromMolFile(input_file)  # fallback
    
    if mol is None:
        raise ValueError(f"Could not read molecule from {input_file}")
    
    # Add explicit hydrogens
    mol_with_h = Chem.AddHs(mol)
    
    # Generate 3D coordinates if not present
    if mol_with_h.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol_with_h)
        AllChem.MMFFOptimizeMolecule(mol_with_h)
    
    # Prepare with Meeko
    preparator = MoleculePreparation()
    preparator.prepare(mol_with_h)
    
    # Write PDBQT
    with open(output_file, 'w') as f:
        f.write(preparator.write_pdbqt_string())


if __name__ == "__main__":
    from get_pdbs import download_and_convert_to_pdb

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
    pubchem_id =  ligands_drugs[ 9 ]["pubchem_id"] # Example PubChem ID for a ligand (Dopamine)
    ligand_folder = "ligand_structures"
    os.makedirs(ligand_folder, exist_ok=True)
    download_and_convert_to_pdb(pubchem_id, ligand_folder)

    ligand_file = os.path.join(ligand_folder, f"ligand_{str(pubchem_id)}.pdb")

    prepare_ligand_with_meeko(ligand_file, os.path.join(ligand_folder, f"docking.{pubchem_id}.meeko.pdbqt"))

drugs_receptor = {
    'receptors': [
        {'drug': 'Acetaminophen', 'target': 'Cyclooxygenase-1 (COX-1)', 'uniprot_id': 'P23219'},
        {'drug': 'Dopamine', 'target': 'Dopamine receptor D1', 'uniprot_id': 'P21728'},
        {'drug': 'Aspirin', 'target': 'Cyclooxygenase-1 (COX-1)', 'uniprot_id': 'P23219'},
        {'drug': 'Ibuprofen', 'target': 'Cyclooxygenase-1 (COX-1)', 'uniprot_id': 'P23219'},
        {'drug': 'Caffeine', 'target': 'Adenosine receptor A2A', 'uniprot_id': 'P29274'},
        {'drug': 'Atenolol', 'target': 'Beta-1 adrenergic receptor', 'uniprot_id': 'P08588'},
        {'drug': 'Loratadine', 'target': 'Histamine H1 receptor', 'uniprot_id': 'P35367'},
        {'drug': 'Omeprazole', 'target': 'H+/K+ ATPase (proton pump)', 'uniprot_id': 'P20616'},
        {'drug': 'Simvastatin', 'target': 'HMG-CoA reductase', 'uniprot_id': 'P04035'},
        {'drug': 'Metformin', 'target': 'AMP-activated protein kinase (AMPK) alpha-1', 'uniprot_id': 'Q13131'}
    ]
}



# %%
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
import os
from rdkit import Chem


def prepare_receptor_with_meeko(input_file, output_file):
    """
    Prepare receptor using Meeko Python API
    
    Args:
        input_file (str): Input PDB file path
        output_file (str): Output PDBQT file path
    """
    
    # Check if input file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file {input_file} not found")
    
    try:
        # Read PDB file using RDKit
        mol = Chem.MolFromPDBFile(input_file, removeHs=False)
        if mol is None:
            raise ValueError("Could not read PDB file with RDKit")
        
        # Prepare molecule using Meeko
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        
        # Write PDBQT file
        writer = PDBQTWriterLegacy()
        pdbqt_string = writer.write_string(preparator.mol)
        
        with open(output_file, 'w') as f:
            f.write(pdbqt_string)
        
        print(f"Receptor prepared successfully with Meeko: {output_file}")
        return True
        
    except Exception as e:
        print(f"Error preparing receptor with Meeko: {e}")
        return False


# Example usage
# prepare_receptor_with_meeko("receptor.pdb", "receptormeeeeeeko.pdbqt")

#%%
# Example usage
# prepare_receptor_with_meeko("receptor.pdb", "receptormeeeeeeko.pdbqt")

# Or with more control
# prepare_receptor_with_autodock(
#     "receptor.pdb", 
#     "receptor.pdbqt",
#     remove_waters=True,
#     add_hydrogens=True
# )

# %%


import subprocess
import os


def prepare_receptor_with_meeko_cli(input_file, output_basename):
    """
    Prepare receptor using Meeko command-line interface
    
    Args:
        input_file (str): Input PDB file path
        output_basename (str): Output basename for PDBQT file
    """
    
    # Check if input file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file {input_file} not found")
    
    try:
        # Prepare command
        cmd = [
            "mk_prepare_receptor.py",
            "--read_pdb", input_file,
            "-o", output_basename,
            "-p"  # Generate receptor PDBQT file
        ]
        
        # Run command
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        output_file = f"{output_basename}"
        print(f"Receptor prepared successfully: {output_file}")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"Error running mk_prepare_receptor.py: {e}")
        print(f"stderr: {e.stderr}")
        return False
    except Exception as e:
        print(f"Error preparing receptor: {e}")
        return False


# Example usage
# folder = "protein_structures"
# prepare_receptor_with_meeko_cli(os.path.join(folder, "AF-P07477-alphafold.pdb"), f"docking.P07477.meeko.pdbqt")
