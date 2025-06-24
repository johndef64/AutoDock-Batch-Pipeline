#%%
# from bioservices import UniProt

# # Initialize UniProt service
# u = UniProt(verbose=False)

# # Search for trypsin in humans
# query = "trypsin AND organism_id:9606"  # 9606 is the taxonomy ID for Homo sapiens
# results = u.search(query, frmt="tsv", columns="id,protein_names,organism_name")

# print(results)

# %%

import requests

def get_uniprot_ids_from_name(protein_name, organism_id):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    query = f"{protein_name} AND organism_id:{organism_id}"
    params = {
        "query": query,
        "format": "tsv",
        "fields": "accession,protein_name,organism_name"
    }
    
    response = requests.get(base_url, params=params)
    
    if response.status_code == 200:
        return response.text
    else:
        return f"Error: HTTP {response.status_code} - {response.text}"

# Get trypsin UniProt IDs for humans
results = get_uniprot_ids_from_name("trypsin", 9606)
import pandas as pd
from io import StringIO
pd.read_table(StringIO(results))


#%%
import urllib.parse
import urllib.request
import os
import json

def get_pdb_or_alphafold(uniprot_id, output_dir=".", prefer_experimental=True):
    """
    Scarica file PDB sperimentali o strutture AlphaFold per un UniProt ID.
    
    Args:
        uniprot_id (str): L'ID UniProt
        output_dir (str): Directory dove salvare i file
        prefer_experimental (bool): Se True, prova prima i PDB sperimentali
    
    Returns:
        dict: Informazioni sui file scaricati
    """
    results = {
        'experimental_pdbs': [],
        'alphafold_structure': None,
        'uniprot_id': uniprot_id
    }
    
    # Crea la directory di output se non esiste
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    if prefer_experimental:
        # Step 1: Prova a scaricare PDB sperimentali
        experimental_files = download_experimental_pdb(uniprot_id, output_dir)
        results['experimental_pdbs'] = experimental_files
        
        if experimental_files:
            print(f"Trovati {len(experimental_files)} file PDB sperimentali per {uniprot_id}")
            return results
        else:
            print(f"Nessun PDB sperimentale trovato per {uniprot_id}. Scarico da AlphaFold...")
    
    # Step 2: Scarica da AlphaFold se non ci sono PDB sperimentali
    alphafold_file = download_alphafold_structure(uniprot_id, output_dir)
    results['alphafold_structure'] = alphafold_file
    
    return results

def download_experimental_pdb(uniprot_id, output_dir):
    """Scarica PDB sperimentali da RCSB PDB"""
    try:
        # Ottieni PDB IDs da UniProt
        url = 'https://www.uniprot.org/uploadlists/'
        params = {
            'from': 'ACC+ID',
            'to': 'PDB_ID',
            'format': 'tab',
            'query': uniprot_id
        }
        
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        
        with urllib.request.urlopen(req) as response:
            pdb_ids_raw = response.read().decode('utf-8')

        # Estrai PDB IDs
        lines = pdb_ids_raw.strip().split('\n')
        pdb_ids = []
        
        for line in lines[1:]:  # Salta l'header
            if '\t' in line:
                parts = line.split('\t')
                if len(parts) >= 2 and parts[1]:
                    pdb_ids.append(parts[1])

        if not pdb_ids:
            return []

        # Scarica i file PDB
        downloaded_files = []
        base_url = "https://files.rcsb.org/download/"
        
        for pdb_id in pdb_ids:
            pdb_file_url = f"{base_url}{pdb_id}.pdb"
            output_path = os.path.join(output_dir, f"{pdb_id}_experimental.pdb")
            
            try:
                print(f"Scaricando PDB sperimentale {pdb_id}...")
                urllib.request.urlretrieve(pdb_file_url, output_path)
                downloaded_files.append(output_path)
                print(f"Scaricato con successo {pdb_id}.pdb")
            except Exception as e:
                print(f"Errore nel scaricare {pdb_id}: {e}")

        return downloaded_files
        
    except Exception as e:
        print(f"Errore nel recuperare PDB sperimentali: {e}")
        return []

def download_alphafold_structure(uniprot_id, output_dir, file_format="pdb"):
    """
    Scarica struttura predetta da AlphaFold
    
    Args:
        uniprot_id (str): UniProt ID
        output_dir (str): Directory di output
        file_format (str): "pdb" o "cif"
    """
    try:
        # URL per scaricare da AlphaFold Database
        if file_format.lower() == "pdb":
            alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
            file_extension = "pdb"
        else:
            alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.cif"
            file_extension = "cif"
        
        output_path = os.path.join(output_dir, f"AF-{uniprot_id}-alphafold.{file_extension}")
        
        print(f"Scaricando struttura AlphaFold per {uniprot_id}...")
        urllib.request.urlretrieve(alphafold_url, output_path)
        print(f"Scaricata struttura AlphaFold: {output_path}")
        
        return output_path
        
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"Struttura AlphaFold non disponibile per {uniprot_id}")
        else:
            print(f"Errore HTTP nel scaricare da AlphaFold: {e}")
        return None
    except Exception as e:
        print(f"Errore nel scaricare da AlphaFold: {e}")
        return None

def convert_cif_to_pdb_with_biopython(cif_file, output_pdb_file):
    """
    Converte file CIF in PDB usando Biopython (richiede: pip install biopython)
    """
    try:
        from Bio.PDB import MMCIFParser, PDBIO
        
        parser = MMCIFParser()
        structure = parser.get_structure('structure', cif_file)
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb_file)
        
        print(f"Convertito {cif_file} in {output_pdb_file}")
        return output_pdb_file
        
    except ImportError:
        print("Biopython non installato. Installa con: pip install biopython")
        return None
    except Exception as e:
        print(f"Errore nella conversione CIF->PDB: {e}")
        return None

# Esempio di utilizzo
if __name__ == "__main__":
    # Esempi con diversi UniProt IDs
    test_ids = [
        "P07477",  # Emoglobina (dovrebbe avere PDB sperimentali)
        "P07478",  # Proteina che potrebbe avere solo AlphaFold
    ]
    
    output_directory = "receptor_structures"
    
    for uniprot_id in test_ids:
        print(f"\n=== Processando {uniprot_id} ===")
        
        results = get_pdb_or_alphafold(uniprot_id, output_directory)
        
        print(f"\nRisultati per {uniprot_id}:")
        if results['experimental_pdbs']:
            print(f"  PDB sperimentali: {len(results['experimental_pdbs'])} file")
            for file_path in results['experimental_pdbs']:
                print(f"    - {file_path}")
        
        if results['alphafold_structure']:
            print(f"  Struttura AlphaFold: {results['alphafold_structure']}")
        
        if not results['experimental_pdbs'] and not results['alphafold_structure']:
            print(f"  Nessuna struttura trovata per {uniprot_id}")

#%%

import urllib.request
import os

def download_alphafold_simple(uniprot_id, output_dir=".", file_format="pdb"):
    """
    Scarica struttura AlphaFold per un UniProt ID
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # URL AlphaFold Database
    if file_format.lower() == "pdb":
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        extension = "pdb"
    else:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.cif"
        extension = "cif"
    
    output_path = os.path.join(output_dir, f"AF-{uniprot_id}.{extension}")
    
    try:
        print(f"Scaricando {uniprot_id} da AlphaFold...")
        urllib.request.urlretrieve(url, output_path)
        print(f"Scaricato: {output_path}")
        return output_path
    except Exception as e:
        print(f"Errore: {e}")
        return None

# Uso
if __name__ == "__main__":
    uniprot_id = "P69905"
    file_path = download_alphafold_simple(uniprot_id)

#%%


import urllib.request
import urllib.parse
import os
import json
import time

def download_ligand_structure(pubchem_id, output_dir=".", file_format="sdf", id_type="cid"):
    """
    Scarica la struttura 3D di un ligando da PubChem
    
    Args:
        pubchem_id (str/int): PubChem ID (CID, SID, o nome)
        output_dir (str): Directory dove salvare il file
        file_format (str): Formato del file ("sdf", "mol", "pdb", "xyz")
        id_type (str): Tipo di ID ("cid", "sid", "name")
    
    Returns:
        str: Path del file scaricato o None se errore
    """
    # Crea directory se non esiste
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # URL base PubChem REST API
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    
    # Costruisci URL per il download
    if file_format.lower() == "pdb":
        # Per formato PDB, usa il servizio 3D
        url = f"{base_url}/compound/{id_type}/{pubchem_id}/record/PDB/?record_type=3d&response_type=save"
        file_extension = "pdb"
    else:
        # Per altri formati
        url = f"{base_url}/compound/{id_type}/{pubchem_id}/record/{file_format.upper()}/?record_type=3d&response_type=save"
        file_extension = file_format.lower()
    
    output_path = os.path.join(output_dir, f"pubchem_{pubchem_id}.{file_extension}")
    
    try:
        print(f"Scaricando ligando PubChem ID {pubchem_id} in formato {file_format.upper()}...")
        
        # Scarica il file
        urllib.request.urlretrieve(url, output_path)
        
        # Verifica che il file non sia vuoto
        if os.path.getsize(output_path) == 0:
            os.remove(output_path)
            print(f"Errore: File vuoto per PubChem ID {pubchem_id}")
            return None
        
        print(f"Scaricato con successo: {output_path}")
        return output_path
        
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"Struttura 3D non disponibile per PubChem ID {pubchem_id}")
        else:
            print(f"Errore HTTP {e.code}: {e}")
        return None
    except Exception as e:
        print(f"Errore nel download: {e}")
        return None


def get_pubchem_info(pubchem_id, id_type="cid"):
    """
    Ottiene informazioni base su un composto da PubChem
    """
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{id_type}/{pubchem_id}/property/MolecularFormula,MolecularWeight,IUPACName/JSON"
        
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode('utf-8'))
        
        if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
            props = data['PropertyTable']['Properties'][0]
            return {
                'cid': props.get('CID'),
                'molecular_formula': props.get('MolecularFormula'),
                'molecular_weight': props.get('MolecularWeight'),
                'iupac_name': props.get('IUPACName', 'N/A')
            }
    except Exception as e:
        print(f"Errore nel recuperare informazioni: {e}")
    
    return None


def search_compound_by_name(compound_name):
    """
    Cerca un composto per nome e restituisce il CID
    """
    try:
        # Codifica il nome per l'URL
        encoded_name = urllib.parse.quote(compound_name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/cids/JSON"
        
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode('utf-8'))
        
        if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
            cids = data['IdentifierList']['CID']
            return cids[0] if cids else None
    except Exception as e:
        print(f"Errore nella ricerca: {e}")
    
    return None


def download_multiple_formats(pubchem_id, output_dir=".", formats=["sdf", "pdb", "mol"]):
    """
    Scarica un ligando in più formati
    """
    results = {}
    
    for fmt in formats:
        print(f"\n--- Scaricando in formato {fmt.upper()} ---")
        file_path = download_ligand_structure(pubchem_id, output_dir, fmt)
        results[fmt] = file_path
        
        # Pausa breve per non sovraccaricare il server
        time.sleep(1)
    
    return results


def convert_name_to_structure(compound_name, output_dir=".", file_format="sdf"):
    """
    Cerca un composto per nome e scarica la sua struttura
    """
    print(f"Cercando '{compound_name}' in PubChem...")
    
    cid = search_compound_by_name(compound_name)
    if not cid:
        print(f"Composto '{compound_name}' non trovato")
        return None
    
    print(f"Trovato CID: {cid}")
    
    # Ottieni informazioni
    info = get_pubchem_info(cid)
    if info:
        print(f"Formula molecolare: {info['molecular_formula']}")
        print(f"Peso molecolare: {info['molecular_weight']}")
        print(f"Nome IUPAC: {info['iupac_name'][:100]}...")
    
    # Scarica struttura
    return download_ligand_structure(cid, output_dir, file_format)



# Esempi di utilizzo
if __name__ == "__main__":
    output_directory = "ligand_structures"
    
    # Esempio 1: Scarica per CID
    print("=== Esempio 1: Aspirina (CID: 2244) ===")
    aspirin_cid = 2244
    
    # Ottieni informazioni
    info = get_pubchem_info(aspirin_cid)
    if info:
        print(f"Informazioni Aspirina:")
        print(f"  Formula: {info['molecular_formula']}")
        print(f"  Peso molecolare: {info['molecular_weight']}")
    
    # Scarica in formato SDF
    sdf_file = download_ligand_structure(aspirin_cid, output_directory, "sdf")
    
    # Scarica in formato PDB
    pdb_file = download_ligand_structure(aspirin_cid, output_directory, "pdb")
    
    print("\n=== Esempio 2: Caffeina per nome ===")
    caffeine_file = convert_name_to_structure("caffeine", output_directory, "sdf")
    
    print("\n=== Esempio 3: Download multipli formati ===")
    # Adenosina (CID: 60961)
    adenosine_files = download_multiple_formats(60961, output_directory, ["sdf", "pdb", "mol"])
    
    print(f"\nFile scaricati per adenosina:")
    for fmt, path in adenosine_files.items():
        if path:
            print(f"  {fmt.upper()}: {path}")
        else:
            print(f"  {fmt.upper()}: Download fallito")



#%%



import urllib.request
import os
from rdkit import Chem
from rdkit.Chem import AllChem

def download_and_convert_to_pdb(pubchem_cid, output_dir=".", add_hydrogens=True):
    """
    Scarica SDF da PubChem e converte in formato PDB
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Step 1: Scarica SDF
    sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/record/SDF/?record_type=3d&response_type=save"
    sdf_path = os.path.join(output_dir, f"temp_{pubchem_cid}.sdf")
    pdb_path = os.path.join(output_dir, f"ligand_{pubchem_cid}.pdb")
    
    try:
        print(f"Scaricando SDF per CID {pubchem_cid}...")
        urllib.request.urlretrieve(sdf_url, sdf_path)
        
        # Step 2: Converti SDF in PDB
        print(f"Convertendo in PDB...")
        mol = Chem.SDMolSupplier(sdf_path)[0]
        
        if mol is None:
            print("Errore: Molecola non valida nell'SDF")
            return None
        
        # Aggiungi idrogeni se richiesto
        if add_hydrogens:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
        
        # Scrivi in formato PDB
        pdb_block = Chem.MolToPDBBlock(mol)
        
        with open(pdb_path, 'w') as f:
            f.write(pdb_block)
        
        # Rimuovi file temporaneo
        os.remove(sdf_path)
        
        print(f"Convertito con successo: {pdb_path}")
        return pdb_path
        
    except Exception as e:
        print(f"Errore nella conversione: {e}")
        if os.path.exists(sdf_path):
            os.remove(sdf_path)
        return None

def download_ligand_pdb_rdkit(pubchem_cid, output_dir=".", compound_name="LIG"):
    """
    Versione migliorata con nome del composto personalizzabile
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Scarica SDF
    sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/record/SDF/?record_type=3d&response_type=save"
    sdf_path = os.path.join(output_dir, f"temp_{pubchem_cid}.sdf")
    pdb_path = os.path.join(output_dir, f"{compound_name}_{pubchem_cid}.pdb")
    
    try:
        urllib.request.urlretrieve(sdf_url, sdf_path)
        
        # Leggi e processa la molecola
        mol = Chem.SDMolSupplier(sdf_path)[0]
        if mol is None:
            raise ValueError("Molecola non valida")
        
        # Aggiungi idrogeni e ottimizza
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        
        # Converti in PDB con nome residuo personalizzato
        pdb_block = Chem.MolToPDBBlock(mol)
        
        # Modifica il nome del residuo nel PDB
        lines = pdb_block.split('\n')
        modified_lines = []
        
        for line in lines:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                # Sostituisci il nome del residuo (colonne 17-20)
                modified_line = line[:17] + compound_name.ljust(3)[:3] + line[20:]
                modified_lines.append(modified_line)
            else:
                modified_lines.append(line)
        
        # Scrivi il file PDB modificato
        with open(pdb_path, 'w') as f:
            f.write('\n'.join(modified_lines))
        
        os.remove(sdf_path)
        print(f"Creato PDB: {pdb_path}")
        return pdb_path
        
    except Exception as e:
        print(f"Errore: {e}")
        if os.path.exists(sdf_path):
            os.remove(sdf_path)
        return None


#%%

import urllib.request
import os
import re

def sdf_to_pdb_simple(pubchem_cid, output_dir=".", residue_name="LIG"):
    """
    Conversione semplice SDF->PDB senza librerie esterne
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Scarica SDF
    sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/record/SDF/?record_type=3d&response_type=save"
    sdf_path = os.path.join(output_dir, f"temp_{pubchem_cid}.sdf")
    pdb_path = os.path.join(output_dir, f"ligand_{pubchem_cid}.pdb")
    
    try:
        urllib.request.urlretrieve(sdf_url, sdf_path)
        
        # Leggi SDF e estrai coordinate
        with open(sdf_path, 'r') as f:
            sdf_content = f.read()
        
        lines = sdf_content.split('\n')
        
        # Trova il blocco delle coordinate (dopo la terza riga)
        if len(lines) < 4:
            raise ValueError("SDF non valido")
        
        # Leggi il numero di atomi dalla quarta riga
        counts_line = lines[3]
        num_atoms = int(counts_line[:3].strip())
        
        # Estrai coordinate degli atomi
        atom_lines = lines[4:4+num_atoms]
        
        # Converti in formato PDB
        pdb_lines = []
        pdb_lines.append("REMARK   Generated from PubChem SDF")
        pdb_lines.append(f"REMARK   PubChem CID: {pubchem_cid}")
        
        for i, line in enumerate(atom_lines, 1):
            parts = line.split()
            if len(parts) >= 4:
                x, y, z, element = parts[0], parts[1], parts[2], parts[3]
                
                # Formato PDB HETATM
                pdb_line = f"HETATM{i:5d}  {element:<3} {residue_name} A   1    {float(x):8.3f}{float(y):8.3f}{float(z):8.3f}  1.00 20.00           {element:>2}"
                pdb_lines.append(pdb_line)
        
        pdb_lines.append("END")
        
        # Scrivi file PDB
        with open(pdb_path, 'w') as f:
            f.write('\n'.join(pdb_lines))
        
        os.remove(sdf_path)
        print(f"Convertito in PDB: {pdb_path}")
        return pdb_path
        
    except Exception as e:
        print(f"Errore nella conversione: {e}")
        if os.path.exists(sdf_path):
            os.remove(sdf_path)
        return None

# Esempi di utilizzo
if __name__ == "__main__":
    output_dir = "ligand_structures"
    pubchem_id = 2244  # Aspirina
    pubchem_id = 709778

    # Metodo 1: Con RDKit (raccomandato)
    print("=== Metodo RDKit ===")
    try:
        pdb_file = download_and_convert_to_pdb(pubchem_id, output_dir)  # Aspirina
    except ImportError:
        print("RDKit non installato. Usa: pip install rdkit")
    
    # Metodo 3: Senza dipendenze
    print("\n=== Metodo semplice ===")
    pdb_file_simple = sdf_to_pdb_simple(pubchem_id, output_dir, "ASP")  # Aspirina

# %%
