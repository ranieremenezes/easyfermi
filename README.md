<p align="left" width="100%">
 <img width="100%" height="250" src="https://github.com/clodoN1109/easyFermi/assets/104923248/d1a25a66-0fc6-4484-93fa-aaa8717f4276">
 <h1>easyFermi</h1>
</p> 

|ci|

.. |ci| image:: https://github.com/ranieremenezes/easyFermi/workflows/CI/badge.svg?branch=main
    :target: https://github.com/ranieremenezes/easyFermi/actions?query=workflow%3ACI+branch%3Amain
    :alt: Test Status

The easiest way to analyze Fermi-LAT data.

<div align="center">

![easyfermi-demo-rect-ezgif com-video-to-gif-converter](https://github.com/ranieremenezes/easyFermi/assets/104923248/6657b52f-9538-40ff-86c0-5b742dd6f2b0)

</div>

easyfermi is a solution to facilitate Fermi-LAT data analysis by providing an intuitive graphical interface to interact with Fermi software analysis tools.

# Requirements

- Linux OS / Mac OS / Windows with <a href="https://learn.microsoft.com/en-us/windows/wsl/install">WSL (Windows Subsystem for Linux)</a>
- [Miniconda 3](https://docs.conda.io/projects/miniconda/en/latest/),
  [Anaconda 3](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Miniforge](https://github.com/conda-forge/miniforge)(recommended) distribution.

# Installation

It is recommended to use a virtual environment.

The following instructions assume an installation of `conda` or `mamba` (a faster version of `conda`).

A virtual environment is an important tool for safely installing the dependencies of an application without inadvertently replacing existing versions that may be needed by other applications or programs.

## Users

1. create and activate the virtual environment

   <pre><code>mamba env create -f environment.yml</code></pre>
   <pre><code>mamba activate easyfermi</code></pre>

3. installing `easyfermi`

   <pre><code>pip install --no-deps easyfermi</code></pre>

- (ONLY FOR WINDOWS) Instal the `libgl1` package

   <pre><code>sudo apt-get install libgl1</code></pre>

- If you want, you can set easyfermi as an environmental variable with e.g. in Bash:

<pre><code>alias easyfermi="miniforge && mamba activate easyfermi && python -c 'import easyfermi'"</code></pre>

## Developers

1. clone this repository
2. create the virtual environment

   <pre><code>mamba env create -f environment.yml</code></pre>

3. install `easyfermi` in editable mode

   <pre><code>pip install -e .</code></pre>

# Upgrading

To upgrade your easyfermi installation to the latest version, run the following command in the easyfermi environment:

   <pre><code>pip install easyfermi --upgrade --no-deps</code></pre>

*Warning*: until a conda official package is available, you might need to update the conda environment first

   <pre><code>mamba env update --file environment.yml --prune</code></pre>

You can check your currently installed version of easyfermi with _pip show_:

   <pre><code>pip show easyfermi</code></pre>

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

<p align="center" width="100%">
 <img width="100%" height="100%" src="/code/resources/images/easyFermiWindow.png">
</p> 

# Tutorials

Check for tutorials on the easyfermi YouTube channel:

 <a href="https://www.youtube.com/channel/UCeLCfEoWasUKky6CPNN_opQ">

  ![youtube](https://img.shields.io/badge/YouTube-FF0000?style=for-the-badge&logo=youtube&logoColor=white)

 </a> 

- Learn how to analyze Fermi-LAT data in 10 min
  
 [![Watch the video](https://img.youtube.com/vi/Ny7aA9EBRUs/hqdefault.jpg)](https://youtu.be/Ny7aA9EBRUs)
- Goodness of fit and advanced configurations
  
 [![Watch the video](https://img.youtube.com/vi/OPMOsheCId8/hqdefault.jpg)](https://youtu.be/OPMOsheCId8)


# Acknowledgements

To acknowledge _easyfermi_ in a publication, please cite  [de Menezes, R (2022)](https://ui.adsabs.harvard.edu/abs/2022arXiv220611272D/abstract).

Since _easyfermi_ relies on _fermipy_, _gammapy_, _astropy_, and _emcee_, please also cite [Wood et al. 2017](https://ui.adsabs.harvard.edu/abs/2017ICRC...35..824W/abstract), [Donath et al. 2023](https://ui.adsabs.harvard.edu/abs/2023A%26A...678A.157D/abstract), [Astropy Collaboration et al. 2018](https://ui.adsabs.harvard.edu/abs/2018AJ....156..123A/abstract), and [Foreman-Mackey et al. 2013](https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract) 


I want to thank Clodomir Vianna for helping me design easyfermi, for making the easyfermi logo, and for the several hours of discussion about this project. Clodomir is the one responsible for making easyfermi user-friendly. Thanks to Fabio Cafardo, Lucas Costa Campos, and Raí Menezes for their help and strong support in this project. A big thanks to Alessandra Azzollini, Douglas Carlos, Kaori Nakashima, Lucas Siconato, Matt Pui, and Romana Grossova, the first users/testers of easyfermi.

 <br>
<p align="center" width="100%">
 <img height="200" src="https://github.com/clodoN1109/easyFermi/assets/104923248/a5fd6166-4dce-475b-92e6-78cbcbcd36af">
</p> 