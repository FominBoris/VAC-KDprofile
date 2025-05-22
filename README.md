# VAC-KDprofile
K-coefficients from temperature profile on Venus in far-IR

ATMOSPHERE.INP - read only the first line, name of the atmospheric file.

`gfortrtan main.f90 VAC.f90`

**output file in the root folder (structure):**
km, ln (P), absorption coeff for 1st k-term, ...

Applicability: 60-100 km. The cutoff will be made later (to keep physical relevant figures)

