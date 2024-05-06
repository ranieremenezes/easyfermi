Custom analysis
===============

.. role::  raw-html(raw)
    :format: html

.. image:: ./easyfermi_window_custom.png
  :width: 700

.. _Custom:

Here we give you some examples on how to customize your analysis with ``easyfermi``. In summary, you simply have to modify a configuration file and upload it in the ``easyfermi`` window as indicated in the figure above.

If you don't know where to start with your own configuration file, you can add the coordinates (or name), desired energy and time ranges, and the paths to the data files under the "Standard" button in the ``easyfermi`` main window and then click in the button "Generate config file". This action will generate a 'config.yaml' file within the output directory, allowing you full flexibility to customize it according to your preferences.

.. note::

   If you check the boxes "Improve resolution" and/or "Improve sensitivity", the config.yaml file generated with the button "Generate config file" will be modified accordingly.


Below we give you some examples of customized configuration files.

Precise selection of time intervals 
-----------------------------------

`Tutorial available on YouTube <https://www.youtube.com/watch?v=BG3ldxJv7t4&ab_channel=easyFermi>`_.

Let's suppose you want to build the average SED for 3 specific time windows for BL Lac. In this case, the only feature you need to add to the "config.yaml" file is the line "filter", as indicated below.


.. code-block::

    data:
      evfile : /home/user/Documentos/GUI/Tutorials/BLLac/Output/list.txt
      scfile : /home/user/Documentos/spacecraft/L240206050150320729A098_SC00.fits

    binning:
      roiwidth   : 15
      binsz      : 0.1
      binsperdec : 8

    selection :
      emin : 100.0
      emax : 300000.0
      zmax    : 90
      evclass : 128
      evtype  : 3
      ra: 330.68038041666665
      dec: 42.277771944444446
      tmin: 638496005
      tmax: 655776005
      filter: '((START>6.40997722E8) && (STOP<6.41582362E8)) || ((START>6.4909916E8) && (STOP<6.4943324E8))  || ((START>6.4943324E8) && (STOP<6.4976732E8))'

    gtlike:
      edisp : True
      irfs : 'P8R3_SOURCE_V3'
      edisp_disable : ['isodiff']
      edisp_bins : -2

    model:
      src_roiwidth : 25
      galdiff  : '/home/user/Documentos/Background_Models/gll_iem_v07.fits'
      isodiff  : '/home/user/Documentos/Background_Models/iso_P8R3_SOURCE_V3_v1.txt'
      catalogs : ['4FGL-DR3']



Customized resolution and better SED on the Galactic plane
----------------------------------------------------------

At the cost of decreasing the sensitivity, you can cut the photons with worst positional reconstruction from your dataset by selecting only the photons lying within the best PSF quartiles (`details here <https://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Data/LAT_DP.html>`_.)

Below we give the example of a config.yaml file generated with the button "Generate config file" for Mrk 501, and then we discuss how you can modify it to analyze the data that better suits your goals.

.. note::

   Even if you are not interested in a better resolution, you can use this method to improve the quality of your low energy (i.e. < 500 MeV) SED data points. For instance, if your target is in the Galactic plane, where the contamination levels are very high at low energies, a standard analysis eventually gives you an SED where the lowest energy data points seem too high to be true (e.g. more than :math:`3\sigma` away from the fitted model). This happens because several badly reconstructed photons that do not belong to your target are being swallowed into your analysis. So if you are analyzing a strong source in the Galactic plane, it is typically a good idea to remove the low-energy photons with the worst reconstruction (i.e. PSF0) from your analysis.

.. code-block::

    data:
      evfile : /home/user/Documentos/GUI/easyFermi/code/LHAASO_counterparts/Output_Mrk501/list.txt
      scfile : /home/user/Documentos/GUI/easyFermi/code/LHAASO_counterparts/spacecraft/L240204110942320729A088_SC00.fits

    binning:
      roiwidth   : 15
      binsz      : 0.1
      binsperdec : 8

    selection :
      emin : 100.0
      emax : 800000.0
      zmax    : 90
      evclass : 128
      evtype  : 48
      ra: 253.46756916666664
      dec: 39.76016888888889
      tmin: 636249601
      tmax: 686275200

    gtlike:
      edisp : True
      irfs : 'P8R3_SOURCE_V3'
      edisp_disable : ['isodiff']
      edisp_bins : -2

    model:
      src_roiwidth : 25
      galdiff  : '/home/user/Documentos/Background_Models/gll_iem_v07.fits'
      isodiff  : '/home/user/Documentos/Background_Models/iso_P8R3_SOURCE_V3_v1.txt'
      catalogs : ['4FGL-DR3']

    components:
      - model:
          galdiff  : '/home/user/Documentos/Background_Models/gll_iem_v07.fits'
          isodiff  : '/home/user/Documentos/Background_Models/iso_P8R3_SOURCE_V3_v1.txt'
        selection:
          emin : 100.0
          emax : 500
          zmax : 90
          evtype : 48
      - model:
          galdiff  : '/home/user/Documentos/Background_Models/gll_iem_v07.fits'
          isodiff  : '/home/user/Documentos/Background_Models/iso_P8R3_SOURCE_V3_v1.txt'
        selection:
          emin : 500
          emax : 1000
          zmax : 100
          evtype : 56
      - model:
          galdiff  : '/home/user/Documentos/Background_Models/gll_iem_v07.fits'
          isodiff  : '/home/user/Documentos/Background_Models/iso_P8R3_SOURCE_V3_v1.txt'
        selection:
          emin : 1000
          emax : 300000.0
          zmax : 105
          evtype : 3

Let's suppose you prefer to include all photons with more than 500 MeV in your analysis (in the configuration file above, we use only the 3 best PSF quartiles, i.e. only 75% of the photons, equivalent to setting `evtype : 56`). The only thing you need to do is to modify the following part of the file:

.. code-block:: yaml
    :emphasize-lines: 7,9
    
    [...]
    - model:
          galdiff  : '/home/user/Documentos/Background_Models/gll_iem_v07.fits'
          isodiff  : '/home/user/Documentos/Background_Models/iso_P8R3_SOURCE_V3_v1.txt'
        selection:
          emin : 500
          emax : 1000
          zmax : 100
          evtype : 3
    [...]
 



