#!/bin/bash
#Automate compilation steps for bloodflow module
cd build/
cmake ..
cmake --build .
cd ../scripts
pip install .
cd ..