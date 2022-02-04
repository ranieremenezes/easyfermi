# easyFermi
The easiest way to analyze Fermi-LAT data.

# Requirements
easyFermi relies on Python 3, Fermitools and Fermipy. 

We recommend the user to install Miniconda 3 or Anaconda 3 to proceed.

To install _Fermitools_ with conda, do:

<pre><code>$ conda create -n fermi -c conda-forge -c fermi fermitools python=3
</code></pre>

If you have problems with the installation of Fermitools, please take a look on the video tutorial here: https://fermi.gsfc.nasa.gov/ssc/data/analysis/video_tutorials/

Now that you installed the _Fermitools_, open the fermi environment in the terminal with:

<pre><code>$ conda activate fermi
</code></pre>

And then install _Fermipy_ with the following command:

<pre><code>$ pip install fermipy
</code></pre>


For more details, check the documentation of Fermipy here: https://fermipy.readthedocs.io/en/latest/install.html


# Installation 

Once you are in the fermi environment and have installed _Fermipy_, do:

<pre><code>$ pip install easyFermi
</code></pre>

# Usage

While in the fermi environment, do:

<pre><code>$ python
>>> import easyFermi
</code></pre>



# Tutorials

You can check _easyFermi_ tutorials on YouTube:

https://www.youtube.com/channel/UCeLCfEoWasUKky6CPNN_opQ
