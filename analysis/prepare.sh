#!/bin/bash
# Preparation of the simulation data for further analysis.
# - Remove solvent
# - Remove PBC artifacts
# - Center the system
# - Align the system to the backbone
# Time: 1 h

repetitions=(4 6)
systems=(WT ST58AN)
path="/home/asanchez/chonky/MD4MIKC/simulation/FUL-AG"

for idx in "${!systems[@]}"; do
    reps=${repetitions[$idx]}
    system=${systems[$idx]}

    for ((repetition=0; repetition<reps; repetition++)); do

        # Logging
        echo "  System $system"
        echo "  Repetition $repetition"
        echo "  Path $path/$system/$repetition"

        # Remove solvent
        printf "non-Water\n" | gmx trjconv \
            -s "$path/$system/$repetition/$system.pdb" \
            -f "$path/$system/$repetition/$system.pdb" \
            -o "$path/$system/$repetition/top.pdb"
        
        # Fix topology. The command to remove solvents transforms the chlorine 
        # atoms from "Cl" to "l", which is not recognized by MDAnalysis.
        sed -i 's/ l$/CL/' "$path/$system/$repetition/top.pdb"

        # Keep system whole
        printf "System\n" | gmx trjconv \
            -s "$path/$system/$repetition/top.pdb" \
            -f "$path/$system/$repetition/$system.xtc" \
            -o "$path/$system/$repetition/whole.xtc" \
            -pbc whole

        # Remove jumps
        printf "System\n" | gmx trjconv \
            -s "$path/$system/$repetition/top.pdb" \
            -f "$path/$system/$repetition/whole.xtc" \
            -o "$path/$system/$repetition/nojump.xtc" \
            -pbc nojump

        # Center system
        printf "Protein\nSystem\n" | gmx trjconv \
            -s "$path/$system/$repetition/top.pdb" \
            -f "$path/$system/$repetition/nojump.xtc" \
            -o "$path/$system/$repetition/centered.xtc" \
            -pbc atom \
            -center \
            -ur compact

        # Align backbone
        printf "Backbone\nSystem\n" | gmx trjconv \
            -s "$path/$system/$repetition/top.pdb" \
            -f "$path/$system/$repetition/centered.xtc" \
            -o "$path/$system/$repetition/aligned.xtc" \
            -fit rot+trans
       
       # Simulation start
       printf "System\n" | gmx trjconv \
           -s "$path/$system/$repetition/top.pdb" \
           -f "$path/$system/$repetition/aligned.xtc" \
           -o "$path/$system/$repetition/start.pdb" \
           -dump 0
      
        # Simulation end
        printf "System\n" | gmx trjconv \
            -s "$path/$system/$repetition/top.pdb" \
            -f "$path/$system/$repetition/aligned.xtc" \
            -o "$path/$system/$repetition/end.pdb" \
            -dump 0 \
            -b 2 \
            -tu us

        # Remove intermediate files
        rm "$path/$system/$repetition/whole.xtc"
        rm "$path/$system/$repetition/nojump.xtc"
        rm "$path/$system/$repetition/centered.xtc"

    done
done