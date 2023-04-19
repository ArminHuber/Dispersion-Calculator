# Dispersion Calculator v2.4

Calculate guided wave dispersion diagrams for flat isotropic plates and multilayered anisotropic laminates.

## Description
The MATLAB-based Dispersion Calculator (DC) is an interactive software for calculating the phase velocity, energy velocity, and attenuation dispersion as well as mode shapes of guided waves in flat isotropic plates and multilayered anisotropic laminates. Fluid-loading and viscoelasticity can be considered. Polar dispersion diagrams can be calculated for anisotropic specimens. DC uses the Rayleigh-Lamb equations and the stiffness matrix method (SMM) developed by S. I. Rokhlin and L. Wang (see literature below). DC is continuously improved and validated by using DISPERSE. DC was first released in 2018, and is used worldwide today. For more information, read [DispersionCalculator_Description.pdf](https://github.com/ArminHuber/Dispersion-Calculator/files/11271206/DispersionCalculator_Description.pdf) and [DispersionCalculator_Manual.pdf](https://github.com/ArminHuber/Dispersion-Calculator/files/11271208/DispersionCalculator_Manual.pdf)

## Download and usage

There are two ways of how to use the DC.

### Installing the DC as a stand-alone application (no MATLAB required)

* Download DC_v[xy]_installer.exe.
* Execute DC_v[xy]_installer.exe.

In case the MATLAB runtime is not downloaded and installed automatically:

* Check in ChangeLog.txt for the MATLAB runtime version.
* Download the appropriate MATLAB runtime at https://www.mathworks.com/products/compiler/matlab-runtime.html
* Install the runtime.
* Execute DC_v[xy]_installer.exe.

### Use the code in MATLAB (MATLAB and Curve Fitting Toolbox required)
* Download the repository as a zip-file (click the green "Code" button, and press "Download ZIP").
* Copy the DC_MATLAB_Code folder to your MATLAB working directory.
* In MATLAB, add the DC_MATLAB_Code folder to the MATLAB path (right-click on the folder -> Add to Path -> Selected Folder and Subfolder).
* Open and run DispersionCalculator.m.

## Author

Dr. Armin Huber, armin.huber@dlr.de

## Version History

See ChangeLog.txt.

## License

This project is licensed under the GNU GPL v3 License - see LICENSE.md for details.

## Thanks to

* Prof. Dr. Michael Lowe 
* Prof. Dr. Michel Castaings 
* Prof. Dr. Stanislav Rokhlin 
* Prof. Dr. Marc Deschamps 
* Dr. Eric Ducasse 
* Prof. Dr. Victor Giurgiutiu
* Prof. Dr. Markus Sause

## Literature

### Books

* J. L. Rose, Ultrasonic Waves in Solid Media. Cambridge University Press, Cambridge, 1999.
* A. H. Nayfeh, Wave Propagation in Layered Anisotropic Media with Applications to Composites. North-Holland, Amsterdam, 1995.
* S. I. Rokhlin, D. E. Chimenti and P. B. Nagy, Physical Ultrasonics of Composites. Oxford University Press, Oxford, 2011.

### Journal articles

* S. I. Rokhlin and L. Wang, Stable recursive algorithm for elastic wave propagation in layered anisotropic media: Stiffness matrix method, J. Acoust. Soc. Am., vol. 112, pp. 822-834, 2002.
* L. Wang and S. I. Rokhlin, Stable reformulation of transfer matrix method for wave propagation in layered anisotropic media, Ultrasonics, vol. 39, pp. 413-424, 2001.
* V. G. A. Kamal and V. Giurgiutiu, Stiffness Transfer Matrix Method (STMM) for Stable Dispersion Curves Solution in Anisotropic Composites, in Proc. SPIE, 2014.
* A. M. A. Huber and M. G. R. Sause, Classification of solutions for guided waves in anisotropic composites with large numbers of layers, J. Acoust. Soc. Am., vol. 144, pp. 3236-3251, 2018.

### PhD thesis

* A. M. A. Huber, Numerical Modeling of Guided Waves in Anisotropic Composites with Application to Air-coupled Ultrasonic Inspection. University of Augsburg, Augsburg, 2020. https://elib.dlr.de/139819
