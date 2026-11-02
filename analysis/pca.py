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
PC = 2

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
        in_dir = path.SIMULATION / 'FUL-AG' / system / str(repetition)
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
                select=(
                    '(segid A and resid 75:145 and name CA) or '
                    '(segid B and resid 70:135 and name CA)'
                )
            )

            # Collect coordinates
            backbone = universe.select_atoms('protein and backbone')
            coordinates = backbone.positions.flatten()
            data['Coordinates'].append(coordinates)
            data['Frame'].append(frame.frame)
            data['System'].append(system)
            data['Repetition'].append(repetition)

# PCA
logger.info('Performing PCA...')
X = np.stack(data['Coordinates'])
pca = PCA(n_components=PC)
X_pca = pca.fit_transform(X)
logger.info(f'Explained variance PC1: {pca.explained_variance_ratio_[0]*100:.2f}%')
logger.info(f'Explained variance PC2: {pca.explained_variance_ratio_[1]*100:.2f}%')
np.save(path.PCA / 'explained_variance.npy', pca.explained_variance_ratio_)

# Save reference structure
reference = pca.mean_.reshape(-1, 3)
backbone.positions = reference
backbone.write(path.PCA / 'reference.pdb')
# Principal component motion
for idx in range(PC):
    # Get eigenvector
    eigenvector = pca.components_[idx].reshape(-1, 3)
    logger.info(f'PC{idx+1} eigenvector shape: {eigenvector.shape}')
    # Sample along the principal component
    alphas = np.linspace(
        np.percentile(X_pca[:, idx], 5),
        np.percentile(X_pca[:, idx], 95),
        100
    )
    # Save motions
    out_path = path.PCA / f'PC{idx+1}.xtc'
    with mda.Writer(str(out_path), n_atoms=backbone.n_atoms) as W:
        for alpha in tqdm(alphas, desc=f'Generating PC{idx+1} motion'):
            backbone.positions = reference + alpha * eigenvector
            W.write(backbone)

# Save PCA results
df = pd.DataFrame(X_pca, columns=['PC1', 'PC2'])
df['System'] = data['System']
df['Repetition'] = data['Repetition']
df['Frame'] = data['Frame']
df.to_csv(path.PCA / 'pca.csv', index=False)