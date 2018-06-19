# computing the power spectrum

## goal
The aim of this code is to compute the power spectrum of galaxies from a mock galaxy catalogue in Fourier space and get the correlation function from the power spectrum, as is explained e.g. in [appendix B in this paper (arXiv:1507.01948)][1].

It is written in Fortran and makes use of  OpenMP and FFTW3. A plotting script is also included, which needs python3 and matplotlib.


##  what does this directory contain?
There are two directories:

 - `./code` contains all the used codes
 - `./output_00067` contains the  mock galaxy data and other necessary files.


## How to run the code

`./code/eval_galaxies.sh` is a bash script that compiles `./code/eval_galaxies.f03`, runs it, and then plots the results by calling the `./code/plot_fortran_correlation.py` script.

The scripts are set up to run from this directory (or more precisely, from the directory that contains the `output_00067/` directory). Furthermore, the script expects at least one argument: The number of the output directory. In the case provided here, run it with 

    $ ./code/eval_galaxies.sh 67

An optional cmdline argument is the number of cells to use per dimension to divide the simulation box in. Just pass a second integer as a second cmdline arg, like so:

    $ ./code/eval_galaxies.sh 67 128

to use a box of $128^{3}$ number of cells.


  [1]: https://arxiv.org/pdf/1507.01948.pdf

