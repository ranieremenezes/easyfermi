<p align="left" width="100%">
 <img width="100%" height="250" src="https://github.com/clodoN1109/easyFermi/assets/104923248/d1a25a66-0fc6-4484-93fa-aaa8717f4276">
 <h1>easyfermi</h1>
</p>

<a href="https://easyfermi.readthedocs.io/en/latest/">
 
![easyfermiDocsBadge](https://img.shields.io/badge/docs-easyfermi-green?style=for-the-badge&logo=googledocs&logoColor=white&labelColor=gray&color=blue)

</a>

![ci](https://github.com/ranieremenezes/easyfermi/actions/workflows/ci.yml/badge.svg)

The easiest way to analyze Fermi-LAT data.

<div align="center">

![easyfermi-demo-rect-ezgif com-video-to-gif-converter](https://github.com/ranieremenezes/easyFermi/assets/104923248/6657b52f-9538-40ff-86c0-5b742dd6f2b0)

</div>

_easyfermi_ is a solution to facilitate Fermi-LAT data analysis by providing an intuitive graphical interface to interact with the Fermi science tools.

# Requirements

- Linux OS / Mac OS / Windows with <a href="https://learn.microsoft.com/en-us/windows/wsl/install">WSL (Windows Subsystem for Linux)</a>
- [Miniconda 3](https://docs.conda.io/projects/miniconda/en/latest/),
  [Anaconda 3](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Miniforge](https://github.com/conda-forge/miniforge) (recommended) distribution.

# Installation

The following instructions assume an installation of `conda` or `mamba` (i.e. a faster version of `conda`).

### Conda-based installation 

In the terminal, run:
<pre><code>mamba create --name easyfermi -c conda-forge -c fermi python=3.9 "fermitools>=2.2.0" "healpy=1.16.1" "gammapy=1.1" "fermipy=v1.2" "scipy=1.11.4" "astropy=5.3.3" "pyqt=5.15.9" "astroquery=0.4.6" "psutil=5.9.8" "matplotlib=3.8.2" "emcee=3.1.4" "corner=2.2.2"</code></pre>

This will create the virtual environment and install all dependencies. Then activate the environment and install _easyfermi_:
<pre><code>mamba activate easyfermi</code></pre>
<pre><code>pip install easyfermi</code></pre>

- (ONLY FOR WINDOWS) Install the `libgl1` package:
<pre><code>sudo apt-get install libgl1</code></pre>

- If you want, you can set _easyfermi_ as an environmental variable. For instance, if you use a Bash shell environment, you can open the `.bashrc` file in your home and set:
<pre><code>alias easyfermi="miniforge && mamba activate easyfermi && python -c 'import easyfermi'"</code></pre>
substituting _miniforge_ and _mamba_ by e.g. _anaconda_ and _conda_ if needed. This line of command depends on which distribution of Python you installed and how you set up the _mamba/conda_ environment.

### Cloning the repository installation

Another option to install _easyfermi_ is by cloning the GitHub repository:
<pre><code>git clone https://github.com/ranieremenezes/easyfermi.git</code></pre>

Then change to the _easyfermi_ directory and create the environment:
<pre><code>cd easyfermi</code></pre>
<pre><code>mamba env create -f environment.yml</code></pre>
<pre><code>mamba activate easyfermi</code></pre>

#### For users:

Install _easyfermi_:
<pre><code>pip install --no-deps easyfermi</code></pre>

#### For developers:

Install _easyfermi_ in editable mode:
<pre><code>pip install -e .</code></pre>

# Upgrading

You can check your currently installed version of _easyfermi_ with _pip show_:
<pre><code>pip show easyfermi</code></pre>
   
If you have `easyfermi 2.0.X` installed, upgrade your installation to the latest version by running the following command in the _easyfermi_ environment:
<pre><code>pip install easyfermi --upgrade --no-deps</code></pre>

If instead, you have `easyfermi 1.X.X` installed, please install easyfermi V2 following section **Installation**.


# Uninstalling

In the terminal, run:
<pre><code>mamba deactivate</code></pre>
<pre><code>mamba env remove --name easyfermi</code></pre>

# Running

If you defined the variable _easyfermi_ in your shell environment (see **Installing**), simply type the following in the terminal:
<pre><code>easyfermi</code></pre>

Otherwise, type:
<pre><code>mamba activate easyfermi</code></pre>
<pre><code>python -c "import easyfermi"</code></pre>

Substituting _mamba_ by _conda_ if this is the case for you.

The main window of _easyfermi_ looks like this:
<p align="center" width="100%">
 <img width="100%" height="100%" src="/src/easyfermi/resources/images/easyFermiWindow.png">
</p> 

# Tutorials

Check for tutorials on the _easyfermi_ YouTube channel:

 <a href="https://www.youtube.com/channel/UCeLCfEoWasUKky6CPNN_opQ">

  ![youtube](https://img.shields.io/badge/YouTube-FF0000?style=for-the-badge&logo=youtube&logoColor=white)

 </a> 

- Learn how to analyze Fermi-LAT data in 10 min
  
 [![Watch the video](https://img.youtube.com/vi/Ny7aA9EBRUs/hqdefault.jpg)](https://youtu.be/Ny7aA9EBRUs)
- Goodness of fit and advanced configurations
  
 [![Watch the video](https://img.youtube.com/vi/OPMOsheCId8/hqdefault.jpg)](https://youtu.be/OPMOsheCId8)

# Documentation

The documentation of _easyfermi_ can be found [here](https://easyfermi.readthedocs.io/en/latest/index.html).


# Acknowledgements

To acknowledge _easyfermi_ in a publication, please cite  [de Menezes, R (2022)](https://ui.adsabs.harvard.edu/abs/2022arXiv220611272D/abstract).

Since _easyfermi_ relies on _fermipy_, _gammapy_, _astropy_, and _emcee_, please also cite [Wood et al. 2017](https://ui.adsabs.harvard.edu/abs/2017ICRC...35..824W/abstract), [Donath et al. 2023](https://ui.adsabs.harvard.edu/abs/2023A%26A...678A.157D/abstract), [Astropy Collaboration et al. 2018](https://ui.adsabs.harvard.edu/abs/2018AJ....156..123A/abstract), and [Foreman-Mackey et al. 2013](https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract).

The EBL models adopted in _easyfermi_ are from:
 - Franceschini et al. 2008 (http://adsabs.harvard.edu/abs/2008A%26A...487..837F)
 - Finke et al. 2010 (http://adsabs.harvard.edu/abs/2009arXiv0905.1115F)
 - Dominguez et al. 2011 (http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1007.1459)
 - Franceschini & Rodighiero 2017 (https://ui.adsabs.harvard.edu/abs/2017A%26A...603A..34F/abstract)
 - Saldana-Lopez et al. 2021 (https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.5144S/abstract)

If you make use of EBL correction via _easyfermi_ in your publication, please cite the papers corresponding to the adopted EBL models. The EBL data files in this repository were collected from the _gammapy_ repository at https://github.com/gammapy/gammapy-data/tree/main/ebl

I want to thank Clodomir Vianna for helping me design _easyfermi_, for making the _easyfermi_ logo, and for the several hours of discussion about this project. Clodomir is the one responsible for making _easyfermi_ user-friendly. Thanks to Fabio Cafardo, Lucas Costa Campos, and Raí Menezes for their help and strong support in this project. A big thanks to Alessandra Azzollini, Douglas Carlos, Kaori Nakashima, Lucas Siconato, Matt Pui, and Romana Grossova, the first users/testers of _easyfermi_.

 <br>
<p align="center" width="100%">
 <img height="200" src="https://github.com/clodoN1109/easyFermi/assets/104923248/a5fd6166-4dce-475b-92e6-78cbcbcd36af">
</p> 
