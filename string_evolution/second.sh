#!/bin/bash

#---Here, I have total 27 number of points, loop should run from 2 to P-1 points
#---The new_initial and new_final Python scripts are same as new.py except they also print the initil and final string image CV values

for AT in {2..26}; do 

    if [ "$AT" -eq 2 ]; then
        cp new_initial.py "./md$AT/new.py"
    elif [ "$AT" -eq 26 ]; then
        cp new_final.py "./md$AT/new.py"
    else
        cp new.py "./md$AT/new.py"
    fi

    cd "./md$AT" || exit 1  # safely move into directory

    # Remove lines starting with '#' in COLVAR
    sed -i '/^#/d' COLVAR

    # Run the Python script
    python new.py
    
    echo $AT

    cd ..
    
done

mv ../phi-psi-ite-2 ./
