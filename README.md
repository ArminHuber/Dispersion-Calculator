# Dispersion Calculator v3.3 [![DOI](https://zenodo.org/badge/629601478.svg)](https://doi.org/10.5281/zenodo.23983280)
Calculate guided wave dispersion diagrams for multilayered isotropic plates, rods, and pipes as well as for multilayered anisotropic plates.

## Description
The [MATLAB](https://www.mathworks.com/products/matlab.html)-based Dispersion Calculator (DC) is an interactive software for calculating the phase velocity, energy velocity, and attenuation dispersion as well as mode shapes of guided waves in multilayered isotropic plates, rods, and pipes as well as in multilayered anisotropic plates. Fluid-loading and viscoelasticity can be considered. Polar dispersion diagrams can be calculated for anisotropic plates. DC uses the Rayleigh-Lamb equations for single layer isotropic plates, the corresponding equations for isotropic cylindrical waveguides, and the stiffness matrix method (SMM) combined with the transfer matrix method (TMM) for multilayered plates, rods, and pipes (see literature below). DC is continuously improved and validated by using [DISPERSE](https://www.imperial.ac.uk/non-destructive-evaluation/publications-products-and-services/disperse/). DC was first released in 2018, and is used worldwide today. For more information, read [DispersionCalculator_Description](https://github.com/ArminHuber/Dispersion-Calculator/blob/main/DispersionCalculator_Description.pdf) and [DispersionCalculator_Manual](https://github.com/ArminHuber/Dispersion-Calculator/blob/main/DC_MATLAB_Code/DispersionCalculator_Manual.pdf).

## Download and installation
There are two ways of how to use DC. In the `Releases` section on the right, click on `Dispersion Calculator v3.3`. Here you can 
* Download `DC_v33_Installer.exe` for Windows to install DC as a stand-alone.
* Download `Source code` to use the code in MATLAB.

### Install DC as a stand-alone application on Windows (no MATLAB required)
* Download `DC_v33_Installer.exe`.
* Execute `DC_v33_Installer.exe`.

In case the MATLAB Runtime is not downloaded and installed automatically:
* Download [MATLAB Runtime R2026a (26.1)](https://www.mathworks.com/products/compiler/matlab-runtime.html).
* Install the runtime.
* Execute `DC_v33_Installer.exe`.

### Use the code (MATLAB R2026a or later and Curve Fitting Toolbox required)
* Download `Source code`.
* Copy the `DC_MATLAB_Code` folder to your MATLAB working directory.
* In MATLAB, add the `DC_MATLAB_Code` folder to the MATLAB path (right-click on the folder &rarr; Add to Path &rarr; Selected Folders and Subfolders).
* Open and run `DispersionCalculator.m`.

## Author
Dr. Armin Huber, armin.huber@dlr.de

## Version history
See [ChangeLog](https://github.com/ArminHuber/Dispersion-Calculator/blob/main/ChangeLog.txt).

## License
This project is licensed under the MIT License - see [LICENSE](https://github.com/ArminHuber/Dispersion-Calculator/blob/main/LICENSE) for details.

## Thanks to
* Prof. Dr. Michael Lowe, Imperial College London, London, UK 
* Prof. Dr. Michel Castaings, University of Bordeaux, Bordeaux, France
* Prof. Dr. Stanislav Rokhlin, Ohio State University, Columbus, OH, USA
* Prof. Dr. Marc Deschamps, University of Bordeaux, Bordeaux, France  
* Dr. Eric Ducasse, University of Bordeaux, Bordeaux, France   
* Prof. Dr. Victor Giurgiutiu, University of South Carolina, Columbia, SC, USA
* Prof. Dr. Markus Sause, University of Augsburg, Augsburg, Germany

## Dedication
Dispersion Calculator is dedicated to<br><br>&emsp;&emsp;&emsp;&ensp;&ensp;Sir Patrick Stewart

## Literature

### Books
* S. I. Rokhlin, D. E. Chimenti and P. B. Nagy, *Physical Ultrasonics of Composites* (Oxford University Press, Oxford, 2011).
* J. L. Rose, *Ultrasonic Waves in Solid Media* (Cambridge University Press, Cambridge, 1999).
* A. H. Nayfeh, *Wave Propagation in Layered Anisotropic Media with Applications to Composites* (North-Holland, Amsterdam, 1995).
* B. A. Auld, *Acoustic Fields and Waves in Solids* (Krieger, Malabar, 1990).
* J. D. Achenbach, *Wave Propagation in Elastic Solids* (North-Holland, Amsterdam, 1973). 

### Journal articles
* S. Guo et al., "Guided waves propagation in arbitrarily stacked composite laminates: Between-layers incompatibility issue resolution using hybrid matrix strategy," [Compos. Struct.](https://doi.org/10.1016/j.compstruct.2023.117360) **322**, (2023).
* A. M. A. Huber, "Classification of solutions for guided waves in fluid-loaded viscoelastic composites with large numbers of layers," [J. Acoust. Soc. Am.](https://doi.org/10.1121/10.0020584) **154**(2), 1073–1094 (2023).
* A. M. A. Huber and M. G. R. Sause, "Classification of solutions for guided waves in anisotropic composites with large numbers of layers," [J. Acoust. Soc. Am.](https://doi.org/10.1121/1.5082299) **144**(6), 3236-3251 (2018).
* V. G. A. Kamal and V. Giurgiutiu, "Stiffness transfer matrix method (STMM) for stable dispersion curves solution in anisotropic composites," Proc. SPIE **9064** (2014).
* X. Zhao and J. L. Rose, "Guided circumferential shear horizontal waves in an isotropic hollow cylinder," [J. Acoust. Soc. Am.]( https://doi.org/10.1121/1.1691037) **115**(5), 1912-1916 (2004).
* S. I. Rokhlin and L. Wang, "Stable recursive algorithm for elastic wave propagation in layered anisotropic media: Stiffness matrix method," [J. Acoust. Soc. Am.](https://doi.org/10.1121/1.1497365) **112**(3), 822-834 (2002).
* L. Wang and S. I. Rokhlin, "Stable reformulation of transfer matrix method for wave propagation in layered anisotropic media," [Ultrasonics](https://doi.org/10.1016/S0041-624X(01)00082-8) **39**, 413-424 (2001).
* G. Liu and J. Qu, "Guided Circumferential Waves in a Circular Annulus," [J. Appl. Mech.](https://doi.org/10.1115/1.2789071) **65**(2), 424-430 (1998).
* A. H. Nayfeh, "The general problem of elastic wave propagation in multilayered anisotropic media," [J. Acoust. Soc. Am.](https://doi.org/10.1121/1.400988) **89**(4), 1521–1531 (1991).
* A. H. Nayfeh, "The propagation of horizontally polarized shear waves in multilayered anisotropic media," [J. Acoust. Soc. Am.](https://doi.org/10.1121/1.398580) **86**(5), 2007–2012 (1989).
* R. Kumar, "Flexural Vibrations of Fluid-Filled Circular Cylindrical Shells," Acta Acust. **24**, 137–146 (1971).
* Y.-H. Pao and R. D. Mindlin, "Dispersion of Flexural Waves in an Elastic, Circular Cylinder," [J. Appl. Mech.](https://doi.org/10.1115/1.3644033) **27**(3), 513-520 (1960).
* D. Gazis, "Three‐Dimensional Investigation of the Propagation of Waves in Hollow Circular Cylinders. I. Analytical Foundation," [J. Acoust. Soc. Am.]( https://doi.org/10.1121/1.1907753) **31**(5), 568–573 (1959).

### PhD thesis
* A. M. A. Huber, *Numerical Modeling of Guided Waves in Anisotropic Composites with Application to Air-coupled Ultrasonic Inspection* ([University of Augsburg](https://opus.bibliothek.uni-augsburg.de/opus4/frontdoor/index/index/year/2021/docId/82760), Augsburg, 2020).
