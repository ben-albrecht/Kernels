#!/bin/bash

source config.sh

echo "n = ${n}"
echo "iterations = ${iterations}"

( cd ../SERIAL/Stencil

if [ ! -f stencil ]; then
    make clean &> /dev/null
    make stencil
fi
echo "Running serial C Stencil"
./stencil ${n} ${iterations}

)

( cd ../CHAPEL/Stencil

if [ ! -f stencil-serial ]; then
    make clean &> /dev/null
    make stencil-serial
fi
echo "Running serial Chapel Stencil"
./stencil-serial --n=${n} --iterations=${iterations}
)
