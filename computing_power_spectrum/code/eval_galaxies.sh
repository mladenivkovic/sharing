#!/bin/bash


#==========================================================================
# This scripts compiles and calls the fortran eval_galaxies.f03
# program correctly.
# Assumes you are in the parent directory of output_XXXXX files.
# usage:
# eval_galaxies.sh <output-nr> <ncells>
#   output-nr:  output number of snapshot to work with
#   ncells:     optional; number of cells to divide simulation box in
#               default: 256
#==========================================================================

# set default value for nc
nc=256


if [ $# -lt 3 ]; then
    case $1 in
        ''|*[!0-9]*) echo argument 1 "(" $1 ")" is not an integer?; exit;;
    esac

    output=`printf "output_%05d" "$1"`
    if [ -d "$output" ]; then
        echo "working for directory" $output
    else
        echo "didnt find directory" $output
        exit
    fi

    if [ $# -gt 1 ]; then
        case $2 in
            ''|*[!0-9]*) echo argument 2 "(" $2 ")" is not an integer?; exit;;
        esac
        nc=$2
    fi

else
    echo "I need one argument: the output number to work with"
    echo "optional second argument: ncells"
    exit
fi





#------------------------
# compile f03 program
#------------------------

workdir=$PWD

gfortran $PWD/code/eval_galaxies.f03 -o $PWD/eval_galaxies.o -g -fbacktrace -fopenmp -I/usr/local/include -lfftw3 -lfftw3_omp -Wall

if [ $? -ne 0 ]; then
    echo COMPILATION FAILED. EXITING.
    cd $workdir
    exit
fi








#--------------------------
# call fortran program
#--------------------------

$PWD/eval_galaxies.o $output $nc

if [ $? -ne 0 ]; then
    echo GET_CORRELATION.F03 FAILED. EXITING
    exit
fi



#----------------------------
# Plot Correlations
#----------------------------

$PWD/code/plot_fortran_correlation.py $output
if [ $? -ne 0 ]; then
    echo PLOT_FORTRAN_CORRELATION.PY FAILED. EXITING
    exit
fi

#--------------------------------------
# Show the plot using eye-of-gnome
#--------------------------------------
eog $output/correlation.png

# plot_fortran_galaxies.py
# eog galaxy_density_projection_plot.png
