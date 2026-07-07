# TASM
This repository provides the scripts and workflows required to perform the Multidimensional Temperature-Accelerated String Method (TASM). The method is used to identify multiple transition pathways connecting two known metastable states, compute the corresponding minimum free-energy paths (MFEPs), estimate the associated free-energy landscapes, and calculate transition rates using Kramers theory.

The detailed steps are given in the TASM_workflow PDF document, which follows the steps given in the figure below.

<img width="6031" height="2741" alt="figure-1" src="https://github.com/user-attachments/assets/1081fddb-d8dc-4ff8-ad11-2d71b985a3e0" />


# Package requirement
```
pip3 install numpy
pip3 install os
pip3 install argparse
pip3 install pymbar
pip3 install scipy
```
# Folders
The string_evolution folder has the necessary scripts and code to perform multidimensional TASM

Check the MBAR repository (https://github.com/avijeet0901/MBAR) to compute the weights used in the free-energy calculations.

The post_processing folder contains codes to compute string convergence, exchange probability, free energies from MBAR weights, diffusivity, and effective Kramers rate from the free energy profile. 

The example folder contains the necessary input files to simulate two model systems: alanine-dipeptide and phenylalanine-dipeptide in vacuum.


# Citations
Please cite the following papers if you are using the repository or any segment of this repository:

```
@article{kulshrestha2022finite,
  title={Finite temperature string method with umbrella sampling using path collective variables: application to secondary structure change in a protein},
  author={Kulshrestha, Avijeet and Punnathanam, Sudeep N and Ayappa, K Ganapathy},
  journal={Soft Matter},
  volume={18},
  number={39},
  pages={7593--7603},
  year={2022},
  publisher={Royal Society of Chemistry}
}
@article{kulshrestha2022finite,
  title={Multidimensional Temperature-Accelerated String Method: Sampling Multiple Transition Pathways and Slow Variables},
  author={Kulshrestha, Avijeet},
  journal={ChemRxiv},
  year={2026},
  doi = {10.26434/chemrxiv.15005676/v1}
}

```
# Collaboration
If you encounter any issues or have questions about the workflow, please open a GitHub Issue. For discussions about adapting the workflow to your system of interest or potential collaborations, feel free to contact me at avijeetkulshrestha@gmail.com. I'd be happy to help.

Best Regards,

Avijeet Kulshrestha.


