#!/bin/bash

homepath=/home/e924p726/My-Code-Collection/

rm -ri bin/
mkdir bin/


# Sets up the symlinks
ln -s $homepath/CP2K_Util/wham_generation/setup_wham.py bin/


chmod 777 bin/*
