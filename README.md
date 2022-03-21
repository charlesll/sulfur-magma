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

First line: integer, number of compositions (up to 300). The name is not used, just kept for personal purpose.

Following lines > contain the following information:

SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, MnO, MgO, CaO, Na2O, K2O, P2O5, H2O, S_tot_ppm, so2, T(C), Pbar, xossi_fo2, fs2, index_author, kflag, wmol, kir

SiO2 to H2O: floats, composition in wt%	  

S_tot_ppm: float, ppm S

so2: float, not used

T(C): float, temperature in C

Pbar: float, pressure in bar    

xossi_fo2 : float, oxygen fugacity (log)

fs2: float, sulfur (S2) fugacity (log), initial

index_author: integer, not used, author index for development purpose   

kflag : integer, 0 > use fO2 (xossi) in input, if 2 > use Fe2O3/FeO in input for calculating fO2 (bemol to correct > recalcule quantity of O2 for Fe3+)

wmol : float, guess H2Omol in input, in wt%.

kir : something to do or not to do H2O calc

Please set wmol and kir at 0.01 and 0, respectively. The option is still under development and the feature does not need further exploration for the moment.

For your own convenience you can of course add a note after the last input variable (kir) at  the end of each input string

## OUTPUT FILES

Four output files are produced: ctsfg6.inp, ctsfg6.S2, ctsfg6.so4, ctsfg6.jet

Cstfg6.inp returns the composition normalized to 100% for SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, MnO, MgO, CaO, Na2O, K2O, P2O5, H2O.
Note that sulfur is excluded from normalization.

Ctsfg6.s2 returns the following information:

TK : temperature in Kelvin

Pbars : pressure in bar

logfO2in : input logfO2 (xossi_fO2 in the input file)

logfO2cal : output logfO2 (must be same as logfO2in) if kflag was set to 0 in the input file)

logfS2 : input logfS2

logCs(o) : initial sulfide capacity based on input data

logCs(c) : calculated sulfide capacity (the one you need) after last iteration

logKS2 : the equilibrium constant of the reaction 1/2S2 + O=  <=> 1/2O2 + S=   (Please remember that sulfide capacity is derived from thius quantity)


Ctsfg6.so4  returns the following information:

TK : temperature in Kelvin

Pbars : pressure in bar

logCs(o) : logCs(o) : initial sulfate capacity based on input data

logCs(c) : calculated sulfate capacity (the one you need) after last iteration

LogKSO4 : the equilibrium constant of the reaction 1/2S2 + O= + 3/2O2 <=> SO4=   (Please remember that sulfate capacity is derived from this quantity)

logfO2in : input logfO2 (xossi_fO2 in the input file)

logfO2cal : output logfO2 (must be same as logfO2in if kflag was set to 0 in the input file)

logfS2 : input logfS2

Ctsfg6.jet  returns the following information:

T(K) : temperature in Kelvin

P(bar): pressure in bar

logfO2in : input logfO2 (xossi_fO2 in the input file)

logfO2cal : output logfO2 (must be same as logfO2in if kflag was set to 0 in the input file)

logfS2 : input logfS2

aossi : the activity of O= (see Moretti and Ottonello, 2003 Metall. Mat. Trans. and Moretti and Ottonello, 2005 GCA)

TotAni : the amount of anions

TotCat : the amount of cations

nO= : the amount of free oxygens (O=)

nO- : the amount of non-bridging oxygens (O-)

nO0 : the amount of bridging oxygens (O0)

NBO/T : the amount of non-bridging oxygens over the amount of tetrahedral units

Kpol : the polymerization constant (see Moretti and Ottonello, 2003, 2005 and references therein)

Stot(obs,wt%) : the amount of S (as wt%) in input

Stot(calc,wt%) :  the returned amount of S (as wt%) for input conditions

S_as_S2-(wt%) : the amount of S as sulfide to the calculated amount of bulk sulfur

S_as_S6+(wt%) :  the amount of S as sulfate to the ca lculated amount of bulk sulfur

S6+/tot : the sulfur redox ration (varies between 0 and 1)

log(KSO4/KS2) : the equilibrium constant for reaction S= + 2O2 <=> SO4=

Redox : the iron redox ratio given as log(FeII/FeIII) (used for cross-check of internal features)

Redoz : the iron redox ratio given as log(FeII/FeIII) (same as redox if kflag is set to 0 in input)
         *** USE REDOZ IN ALL CASES, FOR KFLAG EITHER EQUAL TO 0 OR 2 ***

actFe2+ : the ion activity of ion Fe2+

cost_FeO : the equilibrium constant for reaction FeO <=> Fe2+ + O=

kflag : the redox calculation option (0 or 2; see below)


## RUNNING THE PROGRAM

Just write your input. Add as many lines as you wish per run. The number of composition must be indicate in the rist row of input.txt

Type "ctsfg6" (or use any other name you have given to the .exe file while compiling)

Use of "kflag" in input.txt: this can be either 0 or 2 (no other options). If kflag=0 the input logfO2 values will be used, and FeO and Fe2O3 will be calculated accordingly.
If kflag=2 the input FeO and Fe2O3 amounts (i.e. the input FeII/FeIII) will be used and logfO2 calculated accordingly.
A redox value S6+/Stot will be calculated accordingly, which is independent of the amount of sulfur.

Note that the program does not compute fS2 from given "S_tot_ppm" value in input. So, if you want to match that value, you have to change manually "logfS2" in input

As outlined in reference papers (e.g. Moretti and Ottonello, 2005) fO2-fS2 pairs will give you unique pairs of dissolved S and S redox ratio for each P-T-composition condition.

After running you can open ourput files in spreadsheets, plot output variables and play with them.

## REMARKS

CTSFG6 is not intended for inversion of measured S6+/Stot ratios and dissolved S amounts to retrieve fO2 and fS2 values. Rather, it calculates the sulfur redox state and dissolved amount from assigned fO2 and fS2 pairs at given T, P and composition. Therefore, if the user wants to retrieve fO2 from S6+/Stot measurements, he/she must manually change the input logfO2 value until the calculated S6+/Stot ratio that appears in the terminal ouput matches the analytical (EMP or XANES) value. The use of three decimals for logfO2 is encouraged.

Then, to find fS2 from dissolved S-amounts, please manually change the input logfS2 value until the calculated amount of dissolved sulfur that appears in the terminal ouput matches the analytical value for the given redox conditions. Again, the use of three decimals for logfS2 is suggested for ppm-level precision.
