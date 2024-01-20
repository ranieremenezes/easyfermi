<p align="left" width="100%">
 <img width="100%" height="250" src="https://github.com/clodoN1109/easyFermi/assets/104923248/d1a25a66-0fc6-4484-93fa-aaa8717f4276">
 <h1>easyFermi</h1>
</p> 



The easiest way to analyze Fermi-LAT data.
<p align="center" width="100%">
 <img width="60%" height="300" src="https://github.com/clodoN1109/easyFermi/assets/104923248/2c55876f-9983-4164-bd02-de274c9187ed">
</p> 

easyFermi is a solution to facilitate Fermi-LAT data analysis by providing an intuitive graphical interface to interact with Fermi software analysis tools.

# Sections
 - Requirements
 - Installation
 - Uninstallation
 - Running
 - Tutorials
 - Warnings
 - Acknowledgements


# Requirements

- Linux OS / Mac OS
- <a href="https://docs.conda.io/projects/miniconda/en/latest/">Miniconda 3</a> / <a href="https://conda.io/projects/conda/en/latest/user-guide/install/index.html">Anaconda 3</a>

# Installation

- Creating a virtual environment:
<pre>$ <code> conda create --name easyfermi "python=3.9" </code></pre>

- Installing easyfermi and dependencies:
<pre>$ <code>conda install --name easyfermi --channel conda-forge "pyqt=5.15.10" </code></pre>
<pre>$ <code>conda install --name easyfermi --channel conda-forge --channel fermi "easyFermi" "python=3.9" "fermitools>=2.2.0" "healpy" "gammapy" "astropy=5.3.3" "fermipy=v1.2" "astroquery=0.4.6"</code></pre>

# Uninstallation

In the terminal, run:
<pre>$ <code>conda deactivate</code></pre>
<pre>$ <code>conda env remove --name easyfermi</code></pre>

# Running

In the terminal, run:

<pre>$ <code>conda activate easyfermi</code></pre>
<pre>$ <code>python -c "import easyFermi"</code></pre>

<p align="center" width="100%">
 <img width="60%" height=400" src="/code/images/easyFermiWindow.png">
</p> 

# Tutorials

Check for tutorials on the easyFermi YouTube channel:

<a href="https://www.youtube.com/channel/UCeLCfEoWasUKky6CPNN_opQ"> 
 
![youtube](https://img.shields.io/badge/YouTube-FF0000?style=for-the-badge&logo=youtube&logoColor=white)

</a>

# Warnings

### Densely populated areas of the sky

When analyzing a target surrounded by too many $\gamma$-ray sources, as typical in the Galactic plane, the users may face the following warning:

 MINUIT USER ERROR.  PARAMETER NUMBER IS        101
,  ALLOWED RANGE IS ONE TO 100

And the likelihood fit will not work. This happens because the MINUIT minimizer can handle only up to 100 parameters when performing the fit. If this problem happens, we recommend the user to play with the "Free source radius" panel in the main window of _easyFermi_. The user can, for instance, check the box "Only norm." and MINUIT will fit only the normalization of the sources in the ROI. When this box is checked, _easyFermi_ will still fit the spectral shape of strong gamma-ray sources (i.e. above 10 sigma) when optimizing the region of interest, but will skip it during the main fit. In the end, the only real loss when checking the "Only norm." box is in the spectral shape of weak (below 10 sigma) gamma-ray sources.

For more tips on this topic, we recommend the user to take a look at the [Goodness of fit](https://www.youtube.com/watch?v=Ny7aA9EBRUs&t=4s&ab_channel=easyFermi) tutorial on YouTube.



# Acknowledgements

To acknowledge _easyFermi_ in a publiaction, please cite  [de Menezes, R (2022)](https://ui.adsabs.harvard.edu/abs/2022arXiv220611272D/abstract).

Since _easyFermi_ relies on _Fermipy_, please also cite [Wood et al. 2017](https://ui.adsabs.harvard.edu/abs/2017ICRC...35..824W/abstract).

I would like to thank Clodomir Vianna for helping me with the design of easyFermi and for making the easyFermi logo, and to Fabio Cafardo, Lucas Costa Campos and Raí Menezes for their help and strong support in this project. A big thanks to Alessandra Azzollini, Douglas Carlos, Kaori Nakashima, Lucas Siconato, Matt Pui, and Romana Grossova, the first users/testers of easyFermi.

 <br>
<p align="center" width="100%">
 <img height="200" src="https://github.com/clodoN1109/easyFermi/assets/104923248/a5fd6166-4dce-475b-92e6-78cbcbcd36af">
</p> 
