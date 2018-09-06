# sap-taslp-noise-covariance-matrix-estimation-2018
Supporting code to reproduce results of "Noise Covariance Matrix Estimation for Rotating Microphone Arrays".

This code is published as-is in the hope that it will be useful in understanding and reproducing the published results. It is not a toolbox and is not documented as such.

It has been tested on macOS 10.12 with Matlab2017b.

In Matlab do

`run_experiments`

and then

`plot_results`


To follow the overall procedure have a look at `scr_20180326_01_main_loop`. The main functions of interest are likely to be 
`fcn_20170626_01_triple_sph_harm_integral`,
`fcn_20170721_01_get_hrtf_as_SH` and `fcn_20180327_01_ewls_ncm_est`.



Alastair H Moore, September 2018
