#!/bin/bash

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
