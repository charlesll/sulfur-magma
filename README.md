# sulfur-magma

Model for prediction of sulfur redox state in magmatic liquids

Copyrights (2005-2021) Roberto Moretti

# Citation

If you use this code in your work, please cite the following papers:

Moretti, R., & Ottonello, G. (2003). A polymeric approach to the sulfide capacity of silicate slags and melts. Metallurgical and Materials Transactions B, 34(4), 399–410. [doi.org/10.1007/s11663-003-0066-1](https://doi.org/10.1007/s11663-003-0066-1)

Moretti, R., & Ottonello, G. (2005). Solubility and speciation of sulfur in silicate melts: The Conjugated Toop-Samis-Flood-Grjotheim (CTSFG) model. Geochimica et Cosmochimica Acta, 69(4), 801–823. [https://doi.org/10.1016/j.gca.2004.09.006](https://doi.org/10.1016/j.gca.2004.09.006)

# Licence

MIT licence, see Licence file.

# Contributors

Roberto Moretti (IPGP), moretti@ipgp.fr

Charles Le Losq (IPGP), lelosq@ipgp.fr

# Dependencies

A working fortran compiler. We suggest using gfortran, tested on Mac and Linux. It works well with this software!

# How to use

Download the repository, and use the provided example input file. It first requires compilation of the FORTRAN source, then running the compilated software.

## Compilation

To create the program, with gfortran, in the terminal on Linux or MacOS:

`$ gfortran ctsfg6.for -o ctsfg6.o`

## Running the software

The software takes an input file, INPUT.txt, which contains the compositions of interest.

It returns an output file, OUTPUT.txt

Run in the terminal, after compilation, run the command:

`$ ./ctsfg6.o`

## How to use the input file INPUT.txt

This is a comma-separated text file. Please follow the example file (provided) and be careful with integer/float types of variables

First line: integer, number of compositions. The name is not used, just kept for personal purpose.

Following lines > contain the following information:

SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, MnO, MgO, CaO, Na2O, K2O, P2O5, H2O, S_tot_ppm, so2, T(C), Pbar, xossi_fo2, fs2, index_author, kflag, wmol, kir

SiO2 to H2O: floats, composition in wt%	  

S_tot_ppm: float, ppm S

so2: float, not used

T(C): float, temperature in C

Pbar: float, pressure in bar    

xossi_fo2 : float, oxygen fugacity (log)

fs2: float, sulfur fugacity (log), initial (will be adjusted)

index_author: integer, not used, author index for development purpose   

kflag : integer, 0 > use fO2 (xossi) in input, if 2 > use Fe2O3/FeO in input for calculating fO2 (bemol to correct > recalcule quantity of O2 for Fe3+)

wmol : float, guess H2Omol in input, in wt%

kir : something to do or not to do H2O calc

## OUTPUT FILES

To write.