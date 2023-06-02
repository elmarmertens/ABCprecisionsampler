# ABCprecisionsampler: “Precision-based sampling for state space models without measurement error“

This readme file describes the set of replication files for “Precision-based sampling for state space models without measurement error.“[^draft]

[^draft]: [latest draft](https://drive.google.com/file/d/13Hef9bNI3N6PElvCkDbEeyc_E1XUX9Cf/view?usp=share_link), with [supplementary appendix](https://drive.google.com/file/d/1s4qXSRI96sAy36FEXvBavCugDXHxlTZB/view?usp=share_link)


The project is work in progress, and all results are to be considered preliminary.  The materials provided do not necessarily reflect the views of the Deutsche Bundesbank, or the Eurosystem.

## Author

Elmar Mertens (Deutsche Bundesbank) [^em] 

[^em]: Corresponding author: [em@elmarmertens.com](mailto:em@elmarmertens.com)

## Overview

These replication files provide code to apply the precision-based sampler for state spaces that have not measurement error as described in my paper. The code in this repository allows to replicate the following two applications described in the paper:

- Common trend model with VAR(p) dynamics for the gap variables. [goPrecisonsamplerCommonTrendCycle.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/goPrecisonsamplerCommonTrendCycle) simulates data for a multivariate common trend model with VAR(p) cycle and applies the precision-based sampler as well as the Durbin-Koopmans (DK) sampler, and collects execution times for various model configurations described in the paper. [tabulateTrendVAR.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/tabulateTrendVAR.m) tabulates the results.
- Multivariate trend model with VAR(p) dynamics for the gap variables. [goPrecisonsamplerTrendCycle.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/goPrecisonsamplerTrendCycle) simulates data for a multivariate trend model with VAR(p) cycle and applies the precision-based sampler as well as the Durbin-Koopmans (DK) sampler, and collects execution times for various model configurations described in the paper. [tabulateTrendVAR.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/tabulateTrendVAR.m) tabulates the results.
- VAR(p) model with missing observations. [goPrecisonsamplerVARmissingvalues.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/goPrecisonsamplerVARmissingvalues.m) simulates data for a VAR(p) and applies the precision-based sampler as well as the Durbin-Koopmans (DK) sampler, and collects execution times for various model configurations described in the paper. [tabulateVARmissingvalues.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/tabulateVARmissingvalues.m) tabulates the results.
- Estimation of a common trend model as in [Mertens (2016)](https://doi.org/10.1162/REST_a_00549): 
-- Data has been downloaded from FRED and the Philadelphia Fed's website.  [getFREDdata.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/getFREDdata.m)
constructs the data set and produces INFTRMSRV.csv as input file for use by the estimation routines.  [tabulateFREDdata.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/tabulateFREDdata.m) tabulates the variable names, availability etc.
-- [goCommonTrendinflationPS.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/goCommonTrendinflationPS.m) and [goCommonTrendinflationDK.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/goCommonTrendinflationDK.m) perform estimation using either precision-based (PS) or Durbin-Koopman (DK) methods for the state space.
-- results are ploted by [plotCommontrendInflation.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/plotCommontrendInflation.m) and [reportMCMCtocs.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/reportMCMCtocs.m)
In addition, [goPrecisonsamplerTrendVARlike.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/goPrecisonsamplerTrendVARlike.m) illustrates the use of likelihood calculations with the conventional Kalman filter and precision-based methods.
