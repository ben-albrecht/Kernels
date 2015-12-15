#!/bin/bash

source config.sh

echo "n = ${n}"
echo "iterations = ${iterations}"

( cd ../OPENMP/Stencil

if [ ! -f stencil ]; then
    make clean &> /dev/null
    make stencil
fi
echo "Running OPENMP C Stencil"

./stencil ${ppn} ${n} ${iterations}

)

( cd ../CHAPEL/Stencil

if [ ! -f stencil-datablock ]; then
    make clean &> /dev/null
    make stencil-datablock
fi
echo "Running block distributed data-parallel Chapel Stencil"
./stencil-datablock --n=${n} --iterations=${iterations}
)
