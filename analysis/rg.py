"""
===============================================================================
Title:      Rg
Outline:    Calculate the radius of gyration of the protein in the MD 
            simulations.
Author:     Alejandro Sánchez Cano
Date:       21/10/2026
Time:       13 min for 10 simulations
===============================================================================
"""

# Third-party modules
import pandas as pd
from tqdm import tqdm
import MDAnalysis as mda

# Custom modules
from misc import path
from misc.logger import logger
logger.info('Imported modules successfully.')

# Iterate over simulations
repetitions = [4, 6]
systems = ['WT', 'ST58AN']
for reps, system in zip(repetitions, systems):
    for repetition in tqdm(range(reps), desc=f'Processing {system}'):

        # Reconstruct paths
        in_dir = path.SIMULATION / 'FUL-AG' / system / str(repetition)
        xtc = in_dir / 'aligned.xtc'
        pdb = in_dir / 'top.pdb'

        # Create universe
        universe = mda.Universe(pdb, xtc)
        protein = universe.select_atoms("protein")

        # Radius of gyration
        data = {
            'Frame': [],
            'Time': [],
            'Rg': [],
        }
        for idx, ts in enumerate(tqdm(universe.trajectory)):
            data['Rg'].append(protein.radius_of_gyration())
            data['Time'].append(ts.time)
            data['Frame'].append(idx)
        rg = pd.DataFrame.from_dict(data)

        # Save data
        out_file = in_dir / 'rg.csv'
        rg.to_csv(out_file, index=False)