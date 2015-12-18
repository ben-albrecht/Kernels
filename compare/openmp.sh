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

./stencil ${ppn} ${iterations} ${n}

)

( cd ../CHAPEL/Stencil

if [ ! -f stencil-data ]; then
    make clean &> /dev/null
    make stencil-data
fi
echo "Running block distributed data-parallel Chapel Stencil"
./stencil-data  --iterations=${iterations} --n=${n}
)
