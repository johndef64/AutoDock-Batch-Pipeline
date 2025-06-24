#!/usr/bin/env python3
"""
Batch Molecular Docking Script
Runs pairwise docking between lists of receptors and ligands
"""
#%%
import subprocess
import sys
import argparse
from itertools import product
import time
import os

def run_single_docking(receptor_id, ligand_id):
    """
    Run a single docking job using subprocess
    """
    cmd = ["python", "autodock.py", "--receptor", receptor_id, "--ligand", ligand_id]
    
    print(f"🔬 Starting docking: Receptor {receptor_id} vs Ligand {ligand_id}")
    
    try:
        start_time = time.time()
        result = subprocess.run(cmd)#, capture_output=True, text=True, check=True)
        end_time = time.time()
        
        print(f"✅ Completed in {end_time - start_time:.2f} seconds")
        print(f"Output: {result.stdout}")
        
        return {
            'receptor': receptor_id,
            'ligand': ligand_id,
            'status': 'SUCCESS',
            'output': result.stdout,
            'time': end_time - start_time
        }
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed: {e}")
        print(f"Error output: {e.stderr}")
        
        return {
            'receptor': receptor_id,
            'ligand': ligand_id,
            'status': 'FAILED',
            'error': str(e),
            'stderr': e.stderr
        }
run_single_docking("P23219", "3672")


#%%

def run_pairwise_batch(receptors, ligands, mode='pairwise'):
    """
    Run batch docking jobs
    
    Args:
        receptors: List of UniProt IDs
        ligands: List of PubChem IDs
        mode: 'pairwise' (1:1) or 'all_vs_all' (cartesian product)
    """
    results = []
    
    if mode == 'pairwise':
        # Pairwise: receptor[0] vs ligand[0], receptor[1] vs ligand[1], etc.
        if len(receptors) != len(ligands):
            print(f"⚠️  Warning: Receptor list ({len(receptors)}) and ligand list ({len(ligands)}) have different lengths")
            min_len = min(len(receptors), len(ligands))
            receptors = receptors[:min_len]
            ligands = ligands[:min_len]
            print(f"Using first {min_len} pairs")
        
        pairs = list(zip(receptors, ligands))
        
    elif mode == 'all_vs_all':
        # All vs All: every receptor against every ligand
        pairs = list(product(receptors, ligands))
    
    else:
        raise ValueError("Mode must be 'pairwise' or 'all_vs_all'")
    
    print(f"🚀 Starting {len(pairs)} docking jobs in {mode} mode")
    print("=" * 60)
    
    for i, (receptor, ligand) in enumerate(pairs, 1):
        print(f"\n[{i}/{len(pairs)}] Processing pair:")
        result = run_single_docking(receptor, ligand)
        results.append(result)
        print("-" * 40)
    
    return results

def save_results_summary(results, output_file='batch_results.txt'):
    """
    Save a summary of all results to file
    """
    with open(output_file, 'w') as f:
        f.write("BATCH MOLECULAR DOCKING RESULTS\n")
        f.write("=" * 50 + "\n\n")
        
        successful = [r for r in results if r['status'] == 'SUCCESS']
        failed = [r for r in results if r['status'] == 'FAILED']
        
        f.write(f"Total jobs: {len(results)}\n")
        f.write(f"Successful: {len(successful)}\n")
        f.write(f"Failed: {len(failed)}\n\n")
        
        if successful:
            f.write("SUCCESSFUL DOCKINGS:\n")
            f.write("-" * 30 + "\n")
            for result in successful:
                f.write(f"Receptor: {result['receptor']}, Ligand: {result['ligand']}\n")
                f.write(f"Time: {result['time']:.2f}s\n")
                f.write(f"Output: {result['output'][:200]}...\n\n")
        
        if failed:
            f.write("FAILED DOCKINGS:\n")
            f.write("-" * 30 + "\n")
            for result in failed:
                f.write(f"Receptor: {result['receptor']}, Ligand: {result['ligand']}\n")
                f.write(f"Error: {result['error']}\n\n")
    
    print(f"📄 Results summary saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Batch Molecular Docking')
    parser.add_argument('--receptors', nargs='+', required=True,
                       help='List of UniProt IDs for receptors')
    parser.add_argument('--ligands', nargs='+', required=True,
                       help='List of PubChem IDs for ligands')
    parser.add_argument('--mode', choices=['pairwise', 'all_vs_all'], 
                       default='pairwise',
                       help='Docking mode: pairwise (1:1) or all_vs_all (cartesian)')
    parser.add_argument('--output', default='batch_results.txt',
                       help='Output file for results summary')
    
    args = parser.parse_args()
    
    print("🧬 BATCH MOLECULAR DOCKING")
    print(f"Receptors: {args.receptors}")
    print(f"Ligands: {args.ligands}")
    print(f"Mode: {args.mode}")
    
    # Check if autodock.py exists
    if not os.path.exists('autodock.py'):
        print("❌ Error: autodock.py not found in current directory")
        sys.exit(1)
    
    # Run batch docking
    results = run_pairwise_batch(args.receptors, args.ligands, args.mode)
    
    # Save results
    save_results_summary(results, args.output)
    
    # Print final summary
    successful = len([r for r in results if r['status'] == 'SUCCESS'])
    failed = len([r for r in results if r['status'] == 'FAILED'])
    
    print("\n" + "=" * 60)
    print("🎯 BATCH COMPLETED")
    print(f"✅ Successful: {successful}")
    print(f"❌ Failed: {failed}")
    print(f"📊 Success rate: {successful/len(results)*100:.1f}%")

if __name__ == "__main__":
    # Example usage if run without arguments
    if len(sys.argv) == 1:
        print("Example usage:")
        print("python batch_docking.py --receptors P23219 P21728 P29274 --ligands 60961 2244 12345")
        print("python batch_docking.py --receptors P23219 P21728 --ligands 60961 2244 --mode all_vs_all")
        sys.exit(1)
    
    main()
