#!/bin/bash
#Automate compilation steps for flow module
cd build/
cmake ..
cmake --build .
cd ../scripts
pip install .
cd ..