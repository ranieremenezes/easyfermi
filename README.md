# easyFermi
The easiest way to analyze Fermi-LAT data.

# Requirements
_easyFermi_ relies on Python 3, _Fermitools_ and _Fermipy_. 

We recommend the user to install Miniconda 3 or Anaconda 3 before proceeding.

To install _Fermitools_ in the terminal with conda, do:

<pre><code>$ conda create --name fermi -c conda-forge -c fermi python=3.9 "fermitools>=2.2.0" healpy gammapy
</code></pre>

Then activate the fermi environment:

<pre><code>$ conda activate fermi
</code></pre>

And simply install _Fermipy_ and _easyFermi_ with pip:

<pre><code>$ pip install fermipy ipython easyFermi numpy==1.22.4
</code></pre>

# Usage

While in the fermi environment, do:

<pre><code>$ ipython
>>> import easyFermi
</code></pre>

![easyFermi main window](/code/images/easyFermiWindow.png "EasyFermi main window")

Once you finish the analysis, all of the results will be saved in .npy, .fits and .pdf (or .png) files.
You can easily access them with numpy or other softwares like e.g. TopCat.

Let's use the source 4FGL J1229.0+0202 as an example, you can access the light curve data from python by opening the .npy file: 

LC = numpy.load('4fgl_j1229.0+0202_lightcurve.npy')

flux = LC[()]['flux']

fluxError = LC[()]['flux_err']

tmin = LC[()]['tmin_mjd']

tmax = LC[()]['tmax_mjd']

Or you can simply open the file 4fgl_j1229.0+0202_lightcurve.fits directly in TopCat.

# Tutorials

You can check _easyFermi_ tutorials on YouTube:

https://www.youtube.com/channel/UCeLCfEoWasUKky6CPNN_opQ

# Fermipy V1.0.1 light curve problem

In the old version of _Fermipy_ (i.e. V1.0.1, Python 3), the users face a "KeyError: 'fit_success'" issue when trying to build the light curves. 

This issue is solved here:
https://github.com/fermiPy/fermipy/issues/368
