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
        calphas = universe.select_atoms("protein and name CA")

        # RMSF
        analysis = rms.RMSF(calphas).run(verbose=True)
        rmsf = pd.DataFrame({
            'Chain': calphas.chainIDs,
            'Residue': calphas.resids,
            'RMSF': analysis.results.rmsf,
        })

        # Save data
        out_file = in_dir / 'rmsf.csv'
        rmsf.to_csv(out_file, index=False)