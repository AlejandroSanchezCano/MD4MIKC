"""
===============================================================================
Title:      Prepare simulations
Outline:    Prepare the MD simulations for further analysis by removing PBC
            artifacts and aligning the trajectories.
Author:     Alejandro Sánchez Cano
Date:       19/10/2026
Time:       6 h / simulation
===============================================================================
"""

# Built-in modules
import warnings
from pathlib import Path

# Third-party modules
from tqdm import tqdm
import MDAnalysis as mda
from MDAnalysis import transformations, analysis

# Custom modules
from misc import path
from misc.logger import logger
logger.info('Imported modules successfully.')

# Filter warnings
warnings.filterwarnings(
    'ignore',
    message=r'PDB file is missing resid information.*',
    category=UserWarning,
    module="MDAnalysis.topology.PDBParser"
)
warnings.filterwarnings(
    'ignore',
    message=r"Found no information for attr: 'formalcharges'.*",
    category=UserWarning,
    module="MDAnalysis.coordinates.PDB"
)

# Iterate over simulations
repetitions = [4, 6]
systems = ['WT', 'ST58AN']
for reps, system in zip(repetitions, systems):
    for repetition in tqdm(range(reps), desc=f'Processing {system}'):

        # Reconstruct paths
        xtc = (
            path.SIMULATION / 'FUL-SOC1' / system / 
            str(repetition) / f'{system}.xtc'
        )
        pdb = (
            path.SIMULATION / 'FUL-SOC1' / system / 
            str(repetition) / f'{system}.pdb'
        )
        
        # Select protein + ions atoms
        universe = mda.Universe(pdb)
        selection = universe.select_atoms('protein or resname K CL')

        # Save topology
        out_file = pdb.with_name('top.pdb')
        selection.write(out_file)

        #  Create universe
        universe = mda.Merge(selection)
        universe.load_new(xtc)

        # Logging
        natoms = universe.atoms.n_atoms
        nframes = len(universe.trajectory)
        time_step = universe.trajectory.dt
        logger.info(f'Processing {system} repetition {repetition}')
        logger.info(f'Number of atoms: {natoms}')
        logger.info(f'Number of frames: {nframes}')
        logger.info(f'Time step: {time_step}')

        # Remove PBC artifacts
        logger.info('Removing PBC artifacts...')
        universe.guess_TopologyAttrs(
            to_guess=["bonds"],
            vdwradii={"K": 2.75, "Cl": 1.75},
        )
        protein = universe.select_atoms('protein')
        non_protein = universe.select_atoms('not protein')
        transformations = [
            # When more than half of the atoms of a molecule cross the periodic
            # boundary in the visualization, the molecule jumps to the other
            # side of the box. This is a problem especially for multi-chain 
            # proteins. NoJump() ensures continuity between frames without 
            # jumps
            transformations.NoJump(),
            # Unwrapping ensures that a bonded molecule is kept together even
            # if it crosses the periodic boundary
            transformations.unwrap(protein), 
            # After unwrapping, the molecule is whole, but it may jut ouside of
            # the box. Centering the molecule in the box ensures that the
            # molecule is always inside the box
            transformations.center_in_box(protein, wrap=False), # to the center of the box
            # Shits the ions back into the box according to PBC
            transformations.wrap(non_protein)
        ]
        universe.trajectory.add_transformations(*transformations)

        # Align to reference
        logger.info('Aligning trajectory to first frame...')
        universe.trajectory[0]
        reference = mda.Merge(universe.atoms)
        reference.load_new(universe.atoms.positions.copy()[None, :, :])
        out_file = str(xtc.with_name('aligned.xtc'))
        analysis.align.AlignTraj(
            universe,
            reference,
            select='protein and backbone',
            filename=out_file,
            in_memory=False,
        ).run(verbose=True)
