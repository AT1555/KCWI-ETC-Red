# KCWI-ETC-Red

A Python translation of the IDL red side ETC hosted at https://www2.keck.hawaii.edu/inst/kcwi/etc.html. 

Disclaimer: This is simply a translation of the existing ETC from IDL. The red side of KCWI is still very new, therefore the calculations from this program are not expected to be perfectly accurate.  

# Prerequisites

Python>=3.8

astropy>=5.0

matplotlib>=3.0.0

# Installation

Simply run `git clone https://github.com/AT1555/KCWI-ETC-Red` in a directory of your choice to download the ETC and its necessary data files. The program is simply a python script, so no installation is necessary. 

# Running the Program

To run the ETC, simply call `python [PATH TO DIRECTORY]/KCWI-ETC-Red/KCWI-ETC-Red.py slicer grating grat_wave f_lam_index seeing exposure_time *args` where `slicer`, `grating`, `grat_wave` ,`f_lam_index`, `seeing`, and `exposure_time` are requried positional arguments (described below), and `*args` are optional keywords arguments (also described below):

## Required positional arguments:

These arguments are required, and must be entered in order.

`slicer`: choose `L`, `M`, or `S`.

`grating`: choose `BL`, `BM`, `BH1`, `BH2`, `BH3`, `RL`, `RM1`, `RM2`, `RH1`, `RH2`, `RH3`, or `RH4`.

`grat_wave`: Enter the wavelength (in Angstroms) at which your input flux is measured as a `float`.

`f_lam_index`: Enter the power law index of your desired f_lambda spectrum as a `float`.

`seeing`: Enter the expected seeing as a `float` in arcseconds.

`exposure_time`: Enter you desired total exposure time for your sequence of exposures in seconds as a `float`.

## Optional keyword arguments:

These arguments are (mostly) optional and must be entered as `keyword=value`.

`ccd_bin`: Enter `ccd_bin=1x1` or `ccdbin=2x2` (defaults to `ccd_bin=1x1`).

`mag_AB`: Enter `mag_AB=[mag]` where `[mag]` is a floating point number for the AB magnitude of the input spectrum. This (or `flux`) is a required keyword.

`flux`: enter `flux=[value]` where `[value]` is a floating point number for the flux (or flux density) of the input spectrum at the wavelength `grat_wave`. If `emline_width=` is called as an optional keyword, the flux is interperted as an integrated line flux with units of erg s-1 cm-2 A-1. If `emline_width=` is not called as an optional keyword, the flux is interperted as an flux density with units of erg s-1 cm-2.

`sb`: enter `sb=True` to treat the `flux` or `mag_AB` as a surface brightness in units of arcsec-2 (defaults to `sb=False`).

`emline_width`: enter `emline_width=[float]` where `[float]` is a width in angstroms of an emission line. If called, this keyword will cause the ETC to treat the specified flux above as the integrated line flux of a top-hat profile emissoin line centered at `grat_wave` with a width of `emline_width`.

`Nframes`: enter `Nframes=[int]` where `[int]` is an integer to set the number of frames to be taken over the course of your total exposure time (defaults to 1).

`nas`: Activates `nod-and-shuffle` mode by specifying the desired exposure time of a nod-and-shuffle cycle. As of 2023/08/31, this is not a reccomended mode for KCWI Red and shoudl not be used in this calculator. 

`spatial_bin`: enter `spatial_bin=[x],[y]` where `[x]` and `[y]` are floating point values in arcseconds over which the output data are to be binned and summed (defaults to `spatial_bin=1.,1.`).

`spectral_bin`: enter `spectral_bin=[float]` where `[float]` is a floating point number in angstroms over which to bin for S/N calculation (defaults to `spectral_bin=1.`).

`plot`: enter `plot=False` to disable plotting of S/N results.

`plotout`: enter `plotout=[filename]` where `[filename]` is a filename without an extension in which to save the S/N plot as a PDF.

`printout`: enter `printout=[filename]` where `[filename]` is a filename without an extension in which to save the S/N data as a TXT file

# Example

To run a simulation using the small slicer, a wavelngth of interest of 9210 angstroms, the RH2 grating, a spectral index of 0, 0.7 arcsecond seeing, a 3600 second total integration time, an integrated flux of 1e-16 erg s-1 cm-2 over a 10A width, using 12 frames (each 5 minutes long) and save the resulting SNR data table as 'SNR_TABLE.txt' and the resulting SNR Plot as 'SNR_Plot.pdf', run:

`python KCWI-ETC-Red.py S RM2 9210 0 0.7 3600 Nframes=12 flux=1e-16 emline_width=10 printout=SNR_Table plotout=SNR_Plot`

# Acknowledgments

This program is simply a Python translation of the IDL red side ETC hosted at https://www2.keck.hawaii.edu/inst/kcwi/etc.html. Full credit for the core functionality and underlying calculations belongs to the original authors of the IDL version. 

If this program helps you in your research, please consider thanking Dr. Anthony Taylor (University of Texas at Austin, Department of Astronomy) in your acknowledgments. 
