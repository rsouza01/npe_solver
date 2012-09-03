#!/bin/bash

clear
scons

if [ "$?" -ne "0" ]; then
  echo "Compilação mal sucedida."
  exit 1
fi

./build/npe_solver
