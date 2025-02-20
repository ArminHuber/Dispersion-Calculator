# Dispersion Calculator v3.0 [![DOI](https://zenodo.org/badge/629601478.svg)](https://doi.org/10.5281/zenodo.14167633)
Calculate guided wave dispersion diagrams for isotropic plates, rods, and pipes as well as for multilayered anisotropic laminates.

## Description
The [MATLAB](https://www.mathworks.com/products/matlab.html)-based Dispersion Calculator (DC) is an interactive software for calculating the phase velocity, energy velocity, and attenuation dispersion as well as mode shapes of guided waves in isotropic plates, rods, and pipes as well as multilayered anisotropic laminates. Fluid-loading and viscoelasticity can be considered. Polar dispersion diagrams can be calculated for anisotropic specimens. DC uses the Rayleigh-Lamb equations for isotropic plates, the corresponding equations for isotropic cylindrical waveguides, and the stiffness matrix method (SMM) for multilayered anisotropic plates (see literature below). DC is continuously improved and validated by using [DISPERSE](https://www.imperial.ac.uk/non-destructive-evaluation/products-and-services/disperse/). DC was first released in 2018, and is used worldwide today. For more information, read [DispersionCalculator_Description.pdf](https://github.com/ArminHuber/Dispersion-Calculator/blob/main/DispersionCalculator_Description.pdf) and [DispersionCalculator_Manual.pdf](https://github.com/ArminHuber/Dispersion-Calculator/blob/main/DC_MATLAB_Code/DispersionCalculator_Manual.pdf).

## Download and installation
There are two ways of how to use DC. In the `Releases` section on the right, click on `Dispersion Calculator v3.0`. Here you can 
* download `DC_v30_Installer.exe` to install DC as a stand-alone.
* download `Source code` to use the code in MATLAB.

### Install DC as a stand-alone application (no MATLAB required)
* Download `DC_v30_Installer.exe`.
* Execute `DC_v30_Installer.exe`.

In case the MATLAB Runtime is not downloaded and installed automatically:
* Download the [MATLAB Runtime R2024b (24.2)](https://www.mathworks.com/products/compiler/matlab-runtime.html).
* Install the runtime.
* Execute `DC_v30_Installer.exe`.

### Use the code (MATLAB R2024b or later and Curve Fitting Toolbox required)
* Download `Source code`.
* Copy the `DC_MATLAB_Code` folder to your MATLAB working directory.
* In MATLAB, add the `DC_MATLAB_Code` folder to the MATLAB path (right-click on the folder -> Add to Path -> Selected Folders and Subfolders).
* Open and run `DispersionCalculator.m`.

## Author
Dr. Armin Huber, armin.huber@dlr.de

## Version history
See `ChangeLog.txt`.

## License
This project is licensed under the MIT License - see `LICENSE` for details.

## Thanks to
* Prof. Dr. Michael Lowe, Imperial College London, London, UK 
* Prof. Dr. Michel Castaings, University of Bordeaux, Bordeaux, France
* Prof. Dr. Stanislav Rokhlin, Ohio State University, Columbus, OH, USA
* Prof. Dr. Marc Deschamps, University of Bordeaux, Bordeaux, France  
* Dr. Eric Ducasse, University of Bordeaux, Bordeaux, France   
* Prof. Dr. Victor Giurgiutiu, University of South Carolina, Columbia, SC, USA
* Prof. Dr. Markus Sause, University of Augsburg, Augsburg, Germany

## Literature

### Books
* S. I. Rokhlin, D. E. Chimenti and P. B. Nagy, *Physical Ultrasonics of Composites* (Oxford University Press, Oxford, 2011).
* J. L. Rose, *Ultrasonic Waves in Solid Media* (Cambridge University Press, Cambridge, 1999).
* A. H. Nayfeh, *Wave Propagation in Layered Anisotropic Media with Applications to Composites* (North-Holland, Amsterdam, 1995).

### Journal articles
* A. M. A. Huber, "Classification of solutions for guided waves in fluid-loaded viscoelastic composites with large numbers of layers," [J. Acoust. Soc. Am.](https://doi.org/10.1121/10.0020584) **154**(2), 1073–1094 (2023).
* A. M. A. Huber and M. G. R. Sause, "Classification of solutions for guided waves in anisotropic composites with large numbers of layers," [J. Acoust. Soc. Am.](https://doi.org/10.1121/1.5082299) **144**(6), 3236-3251 (2018).
* V. G. A. Kamal and V. Giurgiutiu, "Stiffness transfer matrix method (STMM) for stable dispersion curves solution in anisotropic composites," Proc. SPIE **9064** (2014).
* S. I. Rokhlin and L. Wang, "Stable recursive algorithm for elastic wave propagation in layered anisotropic media: Stiffness matrix method," [J. Acoust. Soc. Am.](https://doi.org/10.1121/1.1497365) **112**(3), 822-834 (2002).
* L. Wang and S. I. Rokhlin, "Stable reformulation of transfer matrix method for wave propagation in layered anisotropic media," [Ultrasonics](https://doi.org/10.1016/S0041-624X(01)00082-8) **39**, 413-424 (2001).

### PhD thesis
* A. M. A. Huber, *Numerical Modeling of Guided Waves in Anisotropic Composites with Application to Air-coupled Ultrasonic Inspection* ([University of Augsburg](https://opus.bibliothek.uni-augsburg.de/opus4/frontdoor/index/index/year/2021/docId/82760), Augsburg, 2020).
