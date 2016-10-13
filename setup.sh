#!/bin/bash
# Script finishes installation after cloning from Github.

# Run install script:
echo "Installing julia packages"
julia install_packages.jl

echo "Installing seaborn for plotting"
pip install --user seaborn

echo "checking that heuristic is installed..."
# Check if TOPTW needs to be downloaded:
res=$(ls -A ./ILS/TOPTW)
if [ res == "" ]; then
	echo "Cloning Heuristic"
	cd ILS/TOPTW
	git clone https://github.com/stefantj/orienteering_heuristics.git
	mv orienteering_heuristics TOPTW
fi
echo "Finished with setup"
