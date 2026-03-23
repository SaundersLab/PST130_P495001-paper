import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re

pdb_dir = 'PDB_files'
pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
pdb_files.sort()
names = [f[:-4] for f in pdb_files]

n = len(pdb_files)
rmsd_matrix_pruned = np.zeros((n, n))
rmsd_matrix_all = np.zeros((n, n))
n_atoms_matrix_pruned = np.zeros((n, n), dtype=int)
n_atoms_matrix_all = np.zeros((n, n), dtype=int)

for i in range(n):
    for j in range(i+1, n):
        pdb1 = os.path.join(pdb_dir, pdb_files[i])
        pdb2 = os.path.join(pdb_dir, pdb_files[j])
        
        print(f"Aligning {names[i]} vs {names[j]}...")
        
        commands = f"""
open {pdb1}
open {pdb2}
matchmaker #2 to #1
exit
"""
        
        cmd_file = f'temp_cmd_{i}_{j}.cxc'
        with open(cmd_file, 'w') as f:
            f.write(commands)
        
        # Run ChimeraX in nogui mode
        result = subprocess.run(
            ['/Applications/ChimeraX-1.10.app/Contents/bin/ChimeraX', '--nogui', '--cmd', f'runscript {cmd_file}'],
            capture_output=True,
            text=True,
            timeout=60
        )
        
        output = result.stdout + result.stderr
        
        # Pattern: "RMSD between 60 pruned atom pairs is 1.209 angstroms; (across all 108 pairs: 9.171)"
        full_match = re.search(
            r'RMSD between (\d+) pruned atom pairs is ([\d.]+) angstroms.*?across all (\d+) pairs: ([\d.]+)',
            output
        )
        
        if full_match:
            n_pruned = int(full_match.group(1))
            rmsd_pruned = float(full_match.group(2))
            n_all = int(full_match.group(3))
            rmsd_all = float(full_match.group(4))
            
            rmsd_matrix_pruned[i, j] = rmsd_pruned
            rmsd_matrix_pruned[j, i] = rmsd_pruned
            rmsd_matrix_all[i, j] = rmsd_all
            rmsd_matrix_all[j, i] = rmsd_all
            n_atoms_matrix_pruned[i, j] = n_pruned
            n_atoms_matrix_pruned[j, i] = n_pruned
            n_atoms_matrix_all[i, j] = n_all
            n_atoms_matrix_all[j, i] = n_all
            
            print(f"  RMSD (pruned): {rmsd_pruned:.2f} Å ({n_pruned} atoms)")
            print(f"  RMSD (all):    {rmsd_all:.2f} Å ({n_all} atoms)")
        
        if os.path.exists(cmd_file):
            os.remove(cmd_file)

np.fill_diagonal(rmsd_matrix_pruned, 0)
np.fill_diagonal(rmsd_matrix_all, 0)

rmsd_matrix_combined = np.tril(rmsd_matrix_pruned) + np.triu(rmsd_matrix_all, k=1)

# np.savetxt('rmsd_matrix_chimerax_pruned.txt', rmsd_matrix_pruned, fmt='%.2f')
# np.savetxt('rmsd_matrix_chimerax_all.txt', rmsd_matrix_all, fmt='%.2f')
# np.savetxt('rmsd_matrix_chimerax_combined.txt', rmsd_matrix_combined, fmt='%.2f')
# np.savetxt('n_atoms_matrix_chimerax_pruned.txt', n_atoms_matrix_pruned, fmt='%d')
# np.savetxt('n_atoms_matrix_chimerax_all.txt', n_atoms_matrix_all, fmt='%d')

fig, ax = plt.subplots(figsize=(12, 10))

mask_lower = np.triu(np.ones_like(rmsd_matrix_combined, dtype=bool), k=1)
mask_upper = np.tril(np.ones_like(rmsd_matrix_combined, dtype=bool))

cmap_spectral = plt.cm.Spectral_r
cmap_cubehelix = plt.cm.cubehelix_r

vmin_pruned = np.nanmin(rmsd_matrix_pruned)
vmax_pruned = np.nanmax(rmsd_matrix_pruned)
vmin_all = np.nanmin(rmsd_matrix_all)
vmax_all = np.nanmax(rmsd_matrix_all)

sns.heatmap(rmsd_matrix_pruned, mask=mask_lower, 
            cmap=cmap_spectral, annot=True, fmt='.1f', 
            square=True, cbar=False, ax=ax,
            vmin=vmin_pruned, vmax=vmax_pruned,
            xticklabels=names, yticklabels=names,
            linewidths=0.5, linecolor='gray')

sns.heatmap(rmsd_matrix_all, mask=mask_upper,
            cmap=cmap_cubehelix, annot=True, fmt='.2g',
            square=True, cbar=False, ax=ax,
            vmin=vmin_all, vmax=vmax_all,
            xticklabels=names, yticklabels=names,
            linewidths=0.5, linecolor='gray')

from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize

cbar_ax1 = fig.add_axes([0.92, 0.55, 0.02, 0.35])
norm = Normalize(vmin=vmin_pruned, vmax=vmax_pruned)
cb1 = ColorbarBase(cbar_ax1, cmap=cmap_spectral, norm=norm, boundaries=np.linspace(vmin_pruned, vmax_pruned, 0.2))
cb1.set_label('RMSD (Å) - Pruned', rotation=270, labelpad=20)

cbar_ax2 = fig.add_axes([0.92, 0.1, 0.02, 0.35])
norm = Normalize(vmin=vmin_all, vmax=vmax_all)
cb2 = ColorbarBase(cbar_ax2, cmap=cmap_cubehelix, norm=norm, boundaries=np.linspace(vmin_all, vmax_all, 5))
cb2.set_label('RMSD (Å) - All pairs', rotation=270, labelpad=20)

plt.sca(ax)
plt.title('Structural Similarity (ChimeraX Matchmaker)\nLower: Pruned atoms | Upper: All atom pairs', 
          pad=20)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout(rect=[0, 0, 0.9, 1])
mpl.rcParams.update({'text.usetex': False, "svg.fonttype": 'none'})
plt.savefig('structural_similarity_chimerax_split.svg', bbox_inches='tight')
plt.savefig('structural_similarity_chimerax_split.png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 8))
sns.heatmap(rmsd_matrix_pruned, xticklabels=names, yticklabels=names, 
            cmap='cubehelix_r', annot=True, fmt='.2g', square=True,
            cbar_kws={'label': 'RMSD (Å)'})
plt.title('Structural Similarity - Pruned Atoms')
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig('structural_similarity_chimerax_pruned.svg')
plt.close()

plt.figure(figsize=(10, 8))
sns.heatmap(rmsd_matrix_all, xticklabels=names, yticklabels=names, 
            cmap='Spectral_r', annot=True, fmt='.2g', square=True,
            cbar_kws={'label': 'RMSD (Å)'})
plt.title('Structural Similarity - All Atom Pairs')
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig('structural_similarity_chimerax_all.svg')
plt.close()
