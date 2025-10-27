#!/bin/bash
cd /var/www/nsf-dev/tensor_scripts_info
rm HexagonalComponentsTest.txt
rm IsotropicComponentsTest.txt
rm OrthorhombicComponentsTest.txt
rm BreakdownTest.txt
matlab -nodisplay  -nosplash -nodesktop -r "try; VectorProjectionDeconvolutionDEV $1; catch; end; quit;"
