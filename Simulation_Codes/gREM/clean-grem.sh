#!/bin/bash

echo "Cleaning up files"
rm *_dump_*
rm dump_analyz.o*
mkdir final_data
mv TandS_STWHAM.out Teff.dat area c2 thick lambda_funcs.dat histfrac_STWHAM.out dump_singlesnapshot.lmpstraj final_data/
mkdir histograms 
mv histogram_* histograms/
echo "Cleanup complete"
