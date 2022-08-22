# easyFermi
The easiest way to analyze Fermi-LAT data.

# Requirements and installation
_easyFermi_ relies on Python 3, _Fermitools_ and _Fermipy_. 

We recommend the user to install Miniconda 3 or Anaconda 3 before proceeding.

To install _Fermitools_ in the terminal with conda, do:

<pre><code>$ conda create --name fermi -c conda-forge -c fermi python=3.9 "fermitools>=2.2.0" healpy gammapy
</code></pre>

Then activate the fermi environment:

<pre><code>$ conda activate fermi
</code></pre>

And simply install _Fermipy_ and _easyFermi_ with pip:

<pre><code>$ pip install fermipy ipython easyFermi
</code></pre>

# Uninstall

To uninstall the _Fermitools_, _Fermipy_ and _easyFermi_, do:

<pre><code>$ conda remove --name fermi --all
</code></pre>

# Usage

While in the fermi environment, do:

<pre><code>$ ipython
>>> import easyFermi
</code></pre>

![easyFermi main window](/code/images/easyFermiWindow.png "EasyFermi main window")

Before running any analysis, **please make sure that none of the working directories have spaces in their names!** This will crash the analysis.

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

# Acknowledging _easyFermi_

To acknowledge _easyFermi_ in a publiaction, please cite  [de Menezes, R (2022)](https://ui.adsabs.harvard.edu/abs/2022arXiv220611272D/abstract).

Since _easyFermi_ relies on _Fermipy_, please also cite [Wood et al. 2017](https://ui.adsabs.harvard.edu/abs/2017ICRC...35..824W/abstract).


I would like to thank Clodomir Vianna for helping me with the design of easyFermi, and to Fabio Cafardo, Lucas Costa Campos and Raí Menezes for their help and strong support in this project. A big thanks to Alessandra Azzollini, Douglas Carlos, Kaori Nakashima, Lucas Siconato, Matt Pui, and Romana Grossova, the first users/testers of easyFermi.

# Caveats

### Fermipy V1.0.1 light curve problem

In the old version of _Fermipy_ (i.e. V1.0.1, Python 3), the users face a "KeyError: 'fit_success'" issue when trying to build the light curves. 

This issue is solved [here](https://github.com/fermiPy/fermipy/issues/368).


### Densely populated areas of the sky

When analyzing a target surrounded by too many $\gamma$-ray sources, as typical in the Galactic plane, the users may face the following warning:

 MINUIT USER ERROR.  PARAMETER NUMBER IS        101
,  ALLOWED RANGE IS ONE TO 100

And the likelihood fit will not work. This happens because the MINUIT minimizer can handle only up to 100 parameters when performing the fit. If this problem happens, we recommend the user to play with the "Free source radius" panel in the main window of _easyFermi_. The user can, for instance, check the box "Only norm." and MINUIT will fit only the normalization of the sources in the ROI. When this box is checked, _easyFermi_ will still fit the spectral shape of strong gamma-ray sources (i.e. above 10 sigma) when optimizing the region of interest, but will skip it during the main fit. In the end, the only real loss when checking the "Only norm." box is in the spectral shape of weak (below 10 sigma) gamma-ray sources.

For more tips on this topic, we recommend the user to take a look at the [Goodness of fit](https://www.youtube.com/watch?v=Ny7aA9EBRUs&t=4s&ab_channel=easyFermi) tutorial on YouTube.



