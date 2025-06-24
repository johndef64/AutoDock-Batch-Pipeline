# Molecular Docking Pipeline

An automated pipeline for batch molecular docking using AutoDock Vina, enabling docking simulations starting from UniProt IDs for receptors and PubChem IDs for ligands.

## Features

- **Automated docking**: Performs molecular docking using AutoDock Vina
- **Automatic preparation**: Downloads and prepares protein structures from AlphaFold and ligands from PubChem automatically
- **Batch processing**: Supports processing of multiple ligand-receptor pairs
- **Meeko integration**: Uses Meeko for optimal ligand and receptor preparation
- **Predefined database**: Includes 10 common drug-target pairs for quick testing

## System Requirements

- **Operating System**: Linux or macOS (limited Windows support)
- **Python**: ==3.10 
- **Main dependencies**:
  - AutoDock Vina
  - RDKit
  - Meeko
  - NumPy

## Installation

### 1. Create Conda Environment

```bash
conda create -n vina python=3.10
conda activate vina
conda config --env --add channels conda-forge
conda install -c conda-forge numpy swig boost-cpp libboost sphinx sphinx_rtd_theme
```

### 2. Install Meeko

```bash
pip install git+https://github.com/forlilab/Meeko.git
```

### 3. Install Additional Dependencies

```bash
pip install rdkit vina
```

## Project Structure

```
molecular-docking-pipeline/
├── main_docking.py          # Main script
├── meeko_converter.py       # Module for Meeko conversion
├── get_pdbs.py             # Module for structure downloads
├── ligand_structures/       # Directory for ligand structures
├── receptor_structures/     # Directory for receptor structures
└── docking_results/        # Directory for docking results
```

## Usage

### Basic Execution

```bash
python main_docking.py
```

### Execution with Custom Parameters

```bash
python main_docking.py --receptor P23219 --ligand 60961
```

### Available Parameters

- `--receptor`: UniProt ID of target receptor (default: P23219)
- `--ligand`: PubChem ID of ligand (default: 2244)

## Predefined Database

The pipeline includes 10 predefined drug-target pairs:

| ID | Drug | Target | UniProt ID |
|----|------|--------|------------|
| 1 | Acetaminophen | Cyclooxygenase-1 (COX-1) | P23219 |
| 2 | Dopamine | Dopamine receptor D1 | P21728 |
| 3 | Aspirin | Cyclooxygenase-1 (COX-1) | P23219 |
| 4 | Ibuprofen | Cyclooxygenase-1 (COX-1) | P23219 |
| 5 | Caffeine | Adenosine receptor A2A | P29274 |
| 6 | Atenolol | Beta-1 adrenergic receptor | P08588 |
| 7 | Loratadine | Histamine H1 receptor | P35367 |
| 8 | Omeprazole | H+/K+ ATPase (proton pump) | P20616 |
| 9 | Simvastatin | HMG-CoA reductase | P04035 |
| 10 | Metformin | AMPK alpha-1 | Q13131 |

## Workflow

1. **Automatic download**: Downloads protein structures from AlphaFold and ligands from PubChem
2. **Preparation**: Converts structures to required PDBQT formats using Meeko
3. **Configuration**: Defines search space for binding site
4. **Docking**: Executes molecular docking simulation
5. **Results**: Saves best poses and generates detailed reports

## Output

The pipeline generates the following output files:

- **Docked structures**: `docking_results_[UniProt]_[PubChem].pdbqt`
- **Energy report**: `docking_report_[UniProt]_[PubChem].txt`
- **Console log**: Binding energies for each generated pose

### Example Output

```
Binding energies (kcal/mol):
Pose 1: -8.245
Pose 2: -7.891
Pose 3: -7.634
Pose 4: -7.289
Pose 5: -6.957
```

## Advanced Configuration

### Docking Parameters

You can modify docking parameters in the code:

```python
# Search space parameters
docker.define_search_space(
    center_x=15.0, center_y=10.0, center_z=5.0,
    size_x=20, size_y=20, size_z=20
)

# Docking parameters
energies = docker.run_docking(exhaustiveness=16, n_poses=5)
```

### Database Customization

To add new drug-target pairs, modify the `ligands_drugs` and `drugs_receptor` dictionaries in the main code.

## Troubleshooting

### Common Errors

- **FileNotFoundError**: Check internet connection for structure downloads
- **Import Error**: Ensure all dependencies are correctly installed
- **Meeko Error**: Verify Meeko installation from GitHub repository

### Debugging

For detailed debugging, check:
- Existence of downloaded PDB files
- Correct conversion to PDBQT format
- Search space parameters


## License

This project is distributed under the MIT License. See the `LICENSE` file for details.

## Citations

If you use this pipeline in your research, please consider citing:

- AutoDock Vina
- Meeko
- AlphaFold Database
- PubChem Database


## Class Documentation

### AutoDockVinaRunner

Main class for handling molecular docking operations:

**Methods:**
- `prepare_receptor(receptor_pdbqt)`: Loads and prepares the receptor protein
- `prepare_ligand(ligand_pdbqt)`: Loads and prepares the ligand
- `define_search_space(center_x, center_y, center_z, size_x, size_y, size_z)`: Defines the binding site search space
- `run_docking(exhaustiveness, n_poses)`: Executes the docking simulation
- `save_results(output_file)`: Saves docking results to file

### Key Functions

- `parse_arguments()`: Parses command line arguments for receptor and ligand IDs
- `main()`: Main function that orchestrates the entire docking workflow
- `is_jupyter_notebook()`: Detects if running in Jupyter environment for proper argument handling

