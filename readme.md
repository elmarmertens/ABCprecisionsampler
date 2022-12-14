# ABCprecisionsampler: “Precision-based sampling for state space models without measurement error“

This readme file describes the set of replication files for “Precision-based sampling for state space models without measurement error.“[^draft]

[^draft]: [latest draft](https://drive.google.com/file/d/1bzRInkpglMolYdZ2sZeJhs_gbJgPtBXU/view?usp=share_link), with [supplementary appendix](https://drive.google.com/file/d/1ISLJHl2r5Fm9-xG0-ncEVteZKk1DktPY/view?usp=share_link)


The project is work in progress, and all results are to be considered preliminary.  The materials provided do not necessarily reflect the views of the Deutsche Bundesbank, or the Eurosystem.

## Author

Elmar Mertens (Deutsche Bundesbank) [^em] 

[^em]: Corresponding author: [em@elmarmertens.com](mailto:em@elmarmertens.com)

## Overview

These replication files provide code to apply the precision-based sampler for state spaces that have not measurement error as described in my paper. The code in this repository allows to replicate the following two applications described in the paper:

- Common-Trend model with VAR(p) dynamics for the gap variables. [goPrecisonsamplerTrendVAR.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/goPrecisonsamplerTrendVAR.m) simulates data for a VAR(p) and applies the precision-based sampler as well as the Durbin-Koopmans (DK) sampler, and collects execution times for various model configurations described in the paper. [tabulateTrendVAR.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/tabulateTrendVAR.m) tabulates the results.
- VAR(p) model with missing observations. [goPrecisonsamplerVARmissingvalues.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/goPrecisonsamplerVARmissingvalues.m) simulates data for a VAR(p) and applies the precision-based sampler as well as the Durbin-Koopmans (DK) sampler, and collects execution times for various model configurations described in the paper. [tabulateVARmissingvalues.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/tabulateVARmissingvalues.m) tabulates the results.

In addition, [goPrecisonsamplerTrendVARlike.m](https://github.com/elmarmertens/kendallcloverCode/blob/kendallcloverCode/goPrecisonsamplerTrendVARlike.m) illustrates the use of likelihood calculations with the conventional Kalman filter and precision-based methods.