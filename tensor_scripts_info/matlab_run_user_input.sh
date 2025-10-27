#!/bin/bash
cd /var/www/nsf-dev/tensor_scripts_info
matlab -nodisplay -nosplash -nodesktop -r "tensorn1=$1;tensorn2=$2;run('VectorProjectionDeconvolutionDEV.m');exit;"
