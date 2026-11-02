"""
===============================================================================
Title:      RMSD
Outline:    Calculate the RMSD of the protein in the MD simulations.
Author:     Alejandro Sánchez Cano
Date:       21/10/2026
Time:       13 min for 10 simulations
===============================================================================
"""

# Third-party modules
import pandas as pd
from tqdm import tqdm
import MDAnalysis as mda
from MDAnalysis.analysis import rms

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

        # RMSD
        analysis = rms.RMSD(
            universe,
            universe,
            select='protein',
            ref_frame=0,
        )
        analysis.run(verbose=True)
        rmsd = pd.DataFrame(
            analysis.results.rmsd,
            columns=['Frame', 'Time', 'RMSD'],
        )

        # Save data
        out_file = in_dir / 'rmsd.csv'
        rmsd.to_csv(out_file, index=False)