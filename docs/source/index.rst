Welcome to easyfermi's documentation!
=====================================

The graphical interface `easyfermi <https://github.com/ranieremenezes/easyfermi>`_ provides a user-friendly way to analyze *Fermi*-LAT data, making use of `fermipy <https://github.com/fermiPy/fermipy>`_, `gammapy <https://gammapy.org/>`_, `astropy <https://www.astropy.org/>`_, and `emcee <https://emcee.readthedocs.io/en/stable/>`_.


**easyfermi** allows for several different analysis on Fermi-LAT data, such as:

* Modeling of the region of interest (RoI) using the `Binned likelihood analsyis <https://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/binned_likelihood_tutorial.html>`_.

* Generating TS and residual maps for the RoI.

* Extracting a spectral energy distribution (SED) for a target.

* Correcting an SED for EBL absorption.

* Generating a light curve for a target (adaptive or constant time bins).

* Localizing a source or fitting its spatial extension.

* Fitting several SED models with a MCMC approach.


Check out the :doc:`usage` section for further information, including
the GUI :ref:`installation`.



Acknowledging *easyfermi*
---------------------

To acknowledge *easyfermi* in a publication please cite `de Menezes, R. 2022 <https://ui.adsabs.harvard.edu/abs/2022A%26C....4000609D/abstract>`_.

Since easyfermi relies on fermipy, gammapy, astropy, and emcee, please also cite `Wood et al. 2017 <https://ui.adsabs.harvard.edu/#abs/2017ICRC...35..824W>`_, `Donath et al. 2023 <https://ui.adsabs.harvard.edu/abs/2023A%26A...678A.157D/abstract>`_, `Astropy Collaboration et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018AJ....156..123A/abstract>`_, and `Foreman-Mackey et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_.

The EBL models adopted in *easyfermi* are from:
 - `Franceschini et al. 2008 <http://adsabs.harvard.edu/abs/2008A%26A...487..837F>`_.
 - `Finke et al. 2010 <http://adsabs.harvard.edu/abs/2009arXiv0905.1115F>`_.
 - `Dominguez et al. 2011 <http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1007.1459>`_.
 - `Franceschini & Rodighiero 2017 <https://ui.adsabs.harvard.edu/abs/2017A%26A...603A..34F/abstract>`_.
 - `Saldana-Lopez et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.5144S/abstract>`_.

If you make use of one of these models in your publication, please cite the corresponding paper.


Contents
--------

.. toctree::

   usage
   Basic_analysis
   lightcurve
   SED
   api
