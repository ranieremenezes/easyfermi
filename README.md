<p align="left" width="100%">
 <img width="100%" height="250" src="https://github.com/clodoN1109/easyFermi/assets/104923248/d1a25a66-0fc6-4484-93fa-aaa8717f4276">
 <h1>easyFermi</h1>
</p> 

The easiest way to analyze Fermi-LAT data.

<div align="center">

![easyfermi-demo-rect-ezgif com-video-to-gif-converter](https://github.com/ranieremenezes/easyFermi/assets/104923248/6657b52f-9538-40ff-86c0-5b742dd6f2b0)

</div>

easyFermi is a solution to facilitate Fermi-LAT data analysis by providing an intuitive graphical interface to interact with Fermi software analysis tools.

# Requirements

- Linux OS / Mac OS / Windows with <a href="https://learn.microsoft.com/en-us/windows/wsl/install">WSL (Windows Subsystem for Linux)</a>
- <a href="https://docs.conda.io/projects/miniconda/en/latest/">Miniconda 3</a> / <a href="https://conda.io/projects/conda/en/latest/user-guide/install/index.html">Anaconda 3</a>

# Installation instructions

It is recommended to use a virtual environment.

The following instructions assume an installation of `conda` or `mamba` (a faster version of `conda`).
If you don't have it, it is recommended to start with the [Miniforge](https://github.com/conda-forge/miniforge) distribution.

A virtual environment is an important tool for safely installing the dependencies of an application without inadvertently replacing existing versions that may be needed by other applications or programs.

## Users

1. create and activate the virtual environment
   1. `conda create --name easyfermi`
   2. `conda activate easyfermi`

2. install the dependencies
   - `conda install --channel fermi --channel conda-forge "fermitools=2.2.0"`
   - `conda install --channel fermi --channel conda-forge "python=3.9" "fermipy=v1.2" "scipy=1.11.4" "astropy=5.3.3" "pyqt=5.15.10" "astroquery=0.4.6" "healpy=1.16.1" "gammapy=1.1" "psutil=5.9.8" "matplotlib=3.8.2"`

3. installing `easyfermi`:
   - `pip install --no-deps easyfermi`

- (ONLY FOR WINDOWS) Instal the `libgl1` package: `sudo apt-get install libgl1`


## Developers

1. clone this repository
2. create the virtual environment (`mamba env create -f environment.yml`)
3. install `easyfermi` in editable mode (`pip install -e .`)

# Uninstalling

In the terminal, run:
- `conda deactivate`
- `conda env remove --name easyfermi`

# Running

In the terminal, run:

- `conda activate easyfermi`
- `python -c "import easyFermi"`

<p align="center" width="100%">
 <img width="60%" height=400" src="/code/images/easyFermiWindow.png">
</p> 

# Tutorials

Check for tutorials on the easyFermi YouTube channel:

 <a href="https://www.youtube.com/channel/UCeLCfEoWasUKky6CPNN_opQ">

  ![youtube](https://img.shields.io/badge/YouTube-FF0000?style=for-the-badge&logo=youtube&logoColor=white)

 </a> 

- Learn how to analyze Fermi-LAT data in 10 min
  
 [![Watch the video](https://img.youtube.com/vi/Ny7aA9EBRUs/hqdefault.jpg)](https://youtu.be/Ny7aA9EBRUs)
- Goodness of fit and advanced configurations
  
 [![Watch the video](https://img.youtube.com/vi/OPMOsheCId8/hqdefault.jpg)](https://youtu.be/OPMOsheCId8)

# Warnings

- Densely populated areas of the sky:

When analyzing a target surrounded by too many $\gamma$-ray sources, as typical in the Galactic plane, the users may face the following warning:

`MINUIT USER ERROR. PARAMETER NUMBER IS 101, ALLOWED RANGE IS ONE TO 100`

And the likelihood fit will not work. This happens because the MINUIT minimizer can handle only up to 100 parameters when performing the fit. If this problem happens, we recommend the user to play with the "Free source radius" panel in the main window of _easyFermi_. The user can, for instance, check the box "Only norm." and MINUIT will fit only the normalization of the sources in the ROI. When this box is checked, _easyFermi_ will still fit the spectral shape of strong gamma-ray sources (i.e. above 10 sigma) when optimizing the region of interest but will skip it during the main fit. In the end, the only real loss when checking the "Only norm." box is in the spectral shape of weak (below 10 sigma) gamma-ray sources.

For more tips on this topic, we recommend the user take a look at the [Goodness of fit](https://www.youtube.com/watch?v=Ny7aA9EBRUs&t=4s&ab_channel=easyFermi) tutorial on YouTube.

# Acknowledgments

To acknowledge _easyFermi_ in a publication, please cite [de Menezes, R (2022)](https://ui.adsabs.harvard.edu/abs/2022arXiv220611272D/abstract).

Since _easyFermi_ relies on _Fermipy_, please also cite [Wood et al. 2017](https://ui.adsabs.harvard.edu/abs/2017ICRC...35..824W/abstract).

I would like to thank Clodomir Vianna for helping me with the design of easyFermi and for making the easyFermi logo, and Fabio Cafardo, Lucas Costa Campos and Raí Menezes for their help and strong support in this project. A big thanks to Alessandra Azzollini, Douglas Carlos, Kaori Nakashima, Lucas Siconato, Matt Pui, and Romana Grossova, the first users/testers of easyFermi.

 <br>
<p align="center" width="100%">
 <img height="200" src="https://github.com/clodoN1109/easyFermi/assets/104923248/a5fd6166-4dce-475b-92e6-78cbcbcd36af">
</p> 