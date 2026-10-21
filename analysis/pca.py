"""
===============================================================================
Title:      PCA
Outline:    Perform principal component analysis on the protein dynamics.
Author:     Alejandro Sánchez Cano
Date:       21/10/2026
Time:       1 min
===============================================================================
"""

# Third-party modules
import numpy as np
import pandas as pd
from tqdm import tqdm
import MDAnalysis as mda
from MDAnalysis.analysis import align
from sklearn.decomposition import PCA

# Custom modules
from misc import path
from misc.logger import logger
logger.info('Imported modules successfully.')

# Variables
EVERY_NTH = 1000

# Collector
data = {
    'Coordinates': [],
    'Frame': [],
    'System': [],
    'Repetition': [],
}

# Iterate over simulations
repetitions = [4, 6]
systems = ['WT', 'ST58AN']
for reps, system in zip(repetitions, systems):
    for repetition in tqdm(range(reps), desc=f'Processing {system}'):

        # Reconstruct paths
        in_dir = path.SIMULATION / 'FUL-SOC1' / system / str(repetition)
        xtc = in_dir / 'aligned.xtc'
        pdb = in_dir / 'top.pdb'
        
        # Create universes
        reference = mda.Universe(pdb)
        universe = mda.Universe(pdb, xtc)

        # Iterate over trajectory
        trajctory = tqdm(
            universe.trajectory[::EVERY_NTH],
            desc='Iterating over trajectory'
        )
        for frame in trajctory:

            # Align to reference
            align.alignto(
                universe,
                reference,
                select='protein and name CA',
            )

            # Collect coordinates
            atoms = universe.select_atoms('protein and backbone')
            coordinates = atoms.positions.flatten()
            data['Coordinates'].append(coordinates)
            data['Frame'].append(frame.frame)
            data['System'].append(system)
            data['Repetition'].append(repetition)

# PCA
logger.info('Performing PCA...')
X = np.stack(data['Coordinates'])
pca = PCA(n_components=2)
Xpca = pca.fit_transform(X)
logger.info(f'Explained variance PC1: {pca.explained_variance_ratio_[0]*100:.2f}%')
logger.info(f'Explained variance PC2: {pca.explained_variance_ratio_[1]*100:.2f}%')

# Save PCA results
df = pd.DataFrame(Xpca, columns=['PC1', 'PC2'])
df['System'] = data['System']
df['Repetition'] = data['Repetition']
df['Frame'] = data['Frame']
df.to_csv(path.DATA / 'pca.csv', index=False)




# Plot
#fig, ax = plt.subplots(figsize=(8, 6))
#for name in set(labels):
#    mask = labels == name
#    if name.startswith('WT'):
#        ax.scatter(
#            X_pca[mask, 0],
#            X_pca[mask, 1],
#            label=name,
#            alpha=0.7,
#            marker='o',
#        )
#    elif name.startswith('ST58AN'):
#        ax.scatter(
#            X_pca[mask, 0],
#            X_pca[mask, 1],
#            label=name,
#            alpha=0.7,
#            marker='s',
#        )
#
#ax.set_xlabel(
#    f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)"
#)
#ax.set_ylabel(
#    f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)"
#)
#
#plt.title('PCA of Protein Dynamics')
#plt.legend()
#plt.savefig('a.png')
