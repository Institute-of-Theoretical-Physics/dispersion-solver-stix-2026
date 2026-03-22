========================================================================
   PW_STIX: Stix dispersion solver for cold plasma waves
========================================================================

[Authors]      Yasuhito Narita (1), Uwe Motschmann (1), and Tohru Hada (2)
[Contact]     y.narita@tu-braunschweig.de
[Affiliation] (1) Institute of Theoretical Physics, Technical University of Braunschweig, Braunschweig, Germany.
              (2) Interdisciplinary Graduate School of Engineering Sciences, Kyushu University, Kasuga City, Fukuoka, Japan
[License]     MIT License
[DOI]         TBD

------------------------------------------------------------------------
1. OVERVIEW
------------------------------------------------------------------------
This repository provides Python implementation of the dispersion
relation solver for cold plasma waves using the Stix parameter method,
which is the widely-used, standard numerical approach with
dielectric tensor elements S, D, P, R, and L. See Stix, T. H., Waves in Plasmas, American Institute of Physics, New York, 1992.

------------------------------------------------------------------------
2. REQUIREMENTS
------------------------------------------------------------------------
- Python 3.x
- NumPy
- Matplotlib (optional, for plotting)

------------------------------------------------------------------------
3. INSTALLATION
------------------------------------------------------------------------
The code pw_stix.py can run on the terminal.

------------------------------------------------------------------------
4. BASIC USAGE
------------------------------------------------------------------------
Step 1. Set three parameters: (1) the ion-to-electron mass ratio,
(2) the ratio of electron cyclotron frequency to electron plasma
frequency, and (3) propagation angle (in degrees) to the
mean magnetic field.

Step 2. Set the range of frequencies in units of electron
plasma frequency. Set the number of frequenct array elements.

Step 3. Set the name of output file (NumPy binary file, *.npz)
and the name of pdf file (graphical output).

Step 4. Run the code on the terminal by typing, e.g.,

\>\> python3 pw_stix.py


------------------------------------------------------------------------
5. CITATION
------------------------------------------------------------------------
If you use this code in your research, please cite it as:

Yasuhito Narita, Uwe Motschmann, and Tohru Hada, 2026, Stix dispersion solver
(version 1.0.0). Zenodo. doi TBD.

------------------------------------------------------------------------
6. LICENSE
------------------------------------------------------------------------
This project is licensed under the MIT License - see the LICENSE file
for details.


