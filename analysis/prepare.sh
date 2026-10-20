#!/bin/bash

repetitions=(4 6)
systems=(WT ST58AN)
path="/home/asanchez/chonky/MD4MIKC/simulation/FUL-SOC1"

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

        # Remove intermediate files
        rm "$path/$system/$repetition/whole.xtc"
        rm "$path/$system/$repetition/nojump.xtc"
        rm "$path/$system/$repetition/centered.xtc"

    done
done