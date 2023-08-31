# KCWI-ETC-Red

A Python translation of the IDL red side calclulated hosted at https://www2.keck.hawaii.edu/inst/kcwi/etc.html.

# Prerequisites

Python>=3.8

astropy>=5.0

matplotlib>=3.0.0

# Running the Program

To run the ETC, simply call `python [PATH TO FILE]/KCWI-ETC-Red.py slicer grating grat_wave f_lam_index seeing exposure_time *args` where `slicer`, `grating`, `grat_wave` ,`f_lam_index`, `seeing`, and `exposure_time` are requried positional arguments (described below), and `*args` are optional keywords arguments (also described below):

## Required positional arguments:

These arguments are required, and must be entered in order.

`slicer`: choose `L`, `M`, or `S`

`grating`: choose `BL`, `BM`, `BH1`, `BH2`, `BH3`, `RL`, `RM1`, `RM2`, `RH1`, `RH2`, `RH3`, or `RH4`

`grat_wave`: Enter the wavelength (in Angstroms) at which your input flux is measured as a `float`

`f_lam_index`: Enter the power law index of your desired f_lambda spectrum as a `float`

`seeing`: Enter the expected seeing as a `float` in arcseconds

`exposure_time`: Enter you desired total exposure time for your sequence of exposures in seconds as a `float`

## Optional keyword arguments:

These arguments are (mostly) optional and must be entered as `keyword=value`

`ccd_bin`: Enter `ccd_bin=1x1` or `ccdbin=2x2` (defaults to `ccd_bin=1x1`)

`mag_AB`: Enter `mag_AB=[mag]` where `[mag]` is a floating point number for the AB magnitude of the input spectrum. This (or `flux`) is a required keyword.

`flux`: enter `flux=[value]` where `[value]` is a floating point number for the flux (or flux density) of the input spectrum at the wavelength `grat_wave`. If `emline_width=` is called as an optional keyword, the flux is interperted as an integrated line flux with units of erg s-1 cm-2 A-1. If `emline_width=` is not called as an optional keyword, the flux is interperted as an flux density with units of erg s-1 cm-2.

`sb`: enter `sb=True` to treat the `flux` or `mag_AB` as a surface brightness in units of arcsec-2 (defaults to `sb=False`)

`emline_width`: enter `emline_width=[float]` where `[float]` is a width in angstroms of an emission line. If called, this keyword will cause the ETC to treat the specified flux above as the integrated line flux of a top-hat profile emissoin line centered at `grat_wave` with a width of `emline_width`

`Nframes`: enter `Nframes=[int]` where `[int]` is an integer to set the number of frames to be taken over the course of your total exposure time (defaults to 1).

`nas`: Activates `nod-and-shuffle` mode by specifying the desired exposure time of a nod-and-shuffle cycle. As of 2023/08/31, this is not a reccomended mode for KCWI Red and shoudl not be used in this calculator. 

`spatial_bin`: enter `spatial_bin=[x],[y]` where `[x]` and `[y]` are floating point values in arcseconds over which the output data are to be binned and summed (defaults to `spatial_bin=1.,1.`)

`spectral_bin`: enter `spectral_bin=[float]` where `[float]` is a floating point number in angstroms over which to bin for S/N calculation (defaults to `spectral_bin=1.`)

`plot`: enter `plot=False` to disable plotting of S/N results

`plotout`: enter `plotout=[filename]` where `[filename]` is a filename without an extension in which to save the S/N plot as a PDF

`printout`: enter `printout=[filename]` where `[filename]` is a filename without an extension in which to save the S/N data as a TXT file
