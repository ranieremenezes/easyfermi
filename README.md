# easyFermi
The easiest way to analyze Fermi-LAT data.

# Requirements
_easyFermi_ relies on Python 3, _Fermitools_ and _Fermipy_. 

We recommend the user to install Miniconda 3 or Anaconda 3 before proceeding.

To install _Fermitools_ with conda, do:

<pre><code>$ conda create -n fermi -c conda-forge -c fermi fermitools python=3 clhep=2.4.4.1
</code></pre>

If you have problems with the installation of Fermitools, please take a look on the video tutorial here: https://fermi.gsfc.nasa.gov/ssc/data/analysis/video_tutorials/

Now that you installed the _Fermitools_, you can install _Fermipy_ with the following command:

<pre><code>$ conda install --name fermi -c conda-forge fermipy matplotlib=3.3.2
</code></pre>

Note that _Fermipy_ doesn't work properly with more recent releases of matplotlib. This is why we are installing matplotlib V3.3.2 here.

For more details, check the documentation of Fermipy here: https://fermipy.readthedocs.io/en/latest/install.html


# Installing easyFermi

Now you just have to open the fermi environment with

<pre><code>$ conda activate fermi
</code></pre>

And then install easyFermi with pip:

<pre><code>$ pip install easyFermi
</code></pre>

# Usage

While in the fermi environment, do:

<pre><code>$ python
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

# Fermipy light curve problem

In _Fermipy_ V1.0.1 (Python 3), some users face a "KeyError: 'fit_success'" issue when trying to build the light curves. 

This issue is solved here:
https://github.com/fermiPy/fermipy/issues/368
