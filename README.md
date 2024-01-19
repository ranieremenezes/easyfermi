<p align="left" width="100%">
 <img width="100%" height="250" src="https://github.com/clodoN1109/easyFermi/assets/104923248/d1a25a66-0fc6-4484-93fa-aaa8717f4276">
 <h1>easyFermi</h1>
</p> 



The easiest way to analyze Fermi-LAT data.
<p align="center" width="100%">
 <img width="60%" height="300" src="https://github.com/clodoN1109/easyFermi/assets/104923248/2c55876f-9983-4164-bd02-de274c9187ed">
</p> 

easyFermi is a GUI solution for speeding up the analysis of Fermi-LAT data by easying the interaction with the Fermi software analysis tools.



# Requirements

- Linux OS / Mac OS
- <a href="https://docs.conda.io/projects/miniconda/en/latest/">Miniconda 3</a> / <a href="https://conda.io/projects/conda/en/latest/user-guide/install/index.html">Anaconda 3</a>


# Installation on Linux OS and Mac OS
Note: for Mac OS with M1 chip, try instead next section.

In the terminal, run the following commands:

- Creating a virtual environment and installing packages:
<pre><code>$ conda create --name fermi --channel conda-forge --channel fermi "python=3.9" "pyqt=5.15.10" "fermitools>=2.2.0" "healpy" "gammapy" "astropy=5.3.3" "fermipy=1.2" "astroquery=0.4.6"</code></pre>

- Activating the virtual environment:

<pre><code>$ conda activate fermi
</code></pre>

- Installing the easyFermi package:

<pre><code>$ pip install "easyFermi"</code></pre>

# Installation on Mac OS with M1 chip

- Creating a virtual environment:
<pre><code>$ conda create --name fermi python=3.9 </code></pre>

- Activating the virtual environment:
<pre><code>$ conda activate fermi</code></pre>

- Installing _pyqt_ package:
<pre><code>$ conda install pyqt </code></pre>

Note: that _pyqt_ is installed before the _Fermitools_ and _Fermipy_, otherwise it will downgrade some matplotlib library and break the installation.

- Installing _Fermitools_ package:
<pre><code>$ conda install -c conda-forge -c fermi "fermitools>=2.2.0" healpy gammapy
</code></pre>

- Installing _Fermipy_ package:
<pre><code>$ pip install fermipy ipython</code></pre>

- Download the script _easyFermi.py_ from GitHub and put it in your working directory.
From this point on, you can follow the **Usage** section.


# Uninstall

In the terminal, run:
<pre><code>$ conda deactivate</code></pre>
<pre><code>$ conda env remove --name fermi</code></pre>

# Usage

While in the fermi environment, do:

<pre><code>$ ipython
>>> import easyFermi
</code></pre>

<p align="center" width="100%">
 <img width="60%" height=400" src="/code/images/easyFermiWindow.png">
 <h1>easyFermi</h1>
</p> 


Before running any analysis, **please make sure that none of the working directories have spaces in their names!** This will crash the analysis.

Once you finish the analysis, all of the results will be saved in .npy, .fits and .pdf (or .png) files.
You can easily access them with numpy or other softwares like e.g. TopCat.

### Light curves

Let's use the source 4FGL J1229.0+0202 as an example, you can access the light curve data from python by opening the .npy file: 

```
LC = numpy.load('4fgl_j1229.0+0202_lightcurve.npy')
flux = LC[()]['flux']
fluxError = LC[()]['flux_err']
tmin = LC[()]['tmin_mjd']
tmax = LC[()]['tmax_mjd']
```

Or you can simply open the file 4fgl_j1229.0+0202_lightcurve.fits directly with TopCat.

### Excess and significance maps

Let's give another example on how to open the fits file of the Excess map:

```
import astropy.io.fits as pyfits
from astropy.wcs import WCS
import matplotlib.pyplot as plt

Excess_data = pyfits.open("Excess_4FGL_j1229.0+0202_pointsource_powerlaw_2.00_residmap.fits")

wcs = WCS(Excess_data[0].header)
significance_map = Excess_data[0].data

plt.figure()
plt.subplot(projection=wcs)
plt.imshow(significance_map,origin="lower",cmap="RdBu_r")
plt.grid()
plt.title("Significance (in $\sigma$)")
plt.colorbar()

wcs2 = WCS(Excess_data[2].header)
Excess_map = Excess_data[2].data

plt.figure()
plt.subplot(projection=wcs2)
plt.imshow(Excess_map,origin="lower",cmap="RdBu_r")
plt.grid()
plt.title("Excess (in counts)")
plt.colorbar()
plt.show()
```
In this way, you can play with the plots as you wish.

# Tutorials

You can check _easyFermi_ tutorials on YouTube:

https://www.youtube.com/channel/UCeLCfEoWasUKky6CPNN_opQ

# Acknowledging _easyFermi_

To acknowledge _easyFermi_ in a publiaction, please cite  [de Menezes, R (2022)](https://ui.adsabs.harvard.edu/abs/2022arXiv220611272D/abstract).

Since _easyFermi_ relies on _Fermipy_, please also cite [Wood et al. 2017](https://ui.adsabs.harvard.edu/abs/2017ICRC...35..824W/abstract).


I would like to thank Clodomir Vianna for helping me with the design of easyFermi and for making the easyFermi logo, and to Fabio Cafardo, Lucas Costa Campos and Raí Menezes for their help and strong support in this project. A big thanks to Alessandra Azzollini, Douglas Carlos, Kaori Nakashima, Lucas Siconato, Matt Pui, and Romana Grossova, the first users/testers of easyFermi.


# Caveats


### Densely populated areas of the sky

When analyzing a target surrounded by too many $\gamma$-ray sources, as typical in the Galactic plane, the users may face the following warning:

 MINUIT USER ERROR.  PARAMETER NUMBER IS        101
,  ALLOWED RANGE IS ONE TO 100

And the likelihood fit will not work. This happens because the MINUIT minimizer can handle only up to 100 parameters when performing the fit. If this problem happens, we recommend the user to play with the "Free source radius" panel in the main window of _easyFermi_. The user can, for instance, check the box "Only norm." and MINUIT will fit only the normalization of the sources in the ROI. When this box is checked, _easyFermi_ will still fit the spectral shape of strong gamma-ray sources (i.e. above 10 sigma) when optimizing the region of interest, but will skip it during the main fit. In the end, the only real loss when checking the "Only norm." box is in the spectral shape of weak (below 10 sigma) gamma-ray sources.

For more tips on this topic, we recommend the user to take a look at the [Goodness of fit](https://www.youtube.com/watch?v=Ny7aA9EBRUs&t=4s&ab_channel=easyFermi) tutorial on YouTube.



 <br>
 <br>
<p align="center" width="100%">
 <img width="20%" height="200" src="https://github.com/clodoN1109/easyFermi/assets/104923248/a5fd6166-4dce-475b-92e6-78cbcbcd36af">
</p> 
