#!/bin/bash
cd /var/www/nsf-dev/assets
matlab -nodisplay -nosplash -nodesktop -r "filename='$1';Res=$2;rockType=$3;tensorN=$4;run('ChristoffelPlotDEV.m');exit;"
