Usage
=====

.. _installation:

Installation
------------

To use ``easyfermi``, first install the fermitools environment using mamba (or conda):

.. code-block:: console

   $ mamba create --name easyfermi -c conda-forge -c fermi python=3.9 "fermitools>=2.2.0" "healpy=1.16.1" "gammapy=1.1" "scipy=1.10.1" "astropy=5.3.3" "pyqt=5.15.9" "astroquery=0.4.6" "psutil=5.9.8" "emcee=3.1.4" "corner=2.2.2"
   
Then activate the environment and install ``fermipy`` and ``easyfermi``:

.. code-block:: console

    $ mamba activate easyfermi
    $ pip install fermipy easyfermi

* (ONLY FOR WINDOWS) Install the libgl1 package:

.. code-block:: console

    $ sudo apt-get install libgl1
    
* If you want, you can set *easyfermi* as an environmental variable. For instance, if you use a Bash shell environment, you can open the .bashrc file in your home and set:

.. code-block:: console

    $ alias easyfermi="mamba activate easyfermi && python -c 'import easyfermi'"
    
substituting miniforge and mamba by e.g. anaconda and conda if needed. This line of command depends on which distribution of Python you installed and how you set up the mamba/conda environment.

Upgrading
----------------

You can check your currently installed version of easyfermi with pip show:

.. code-block:: console

    $ pip show easyfermi
    
If you have **easyfermi 2.0.X** installed, upgrade your installation to the latest version by running the following command in the easyfermi environment:

.. code-block:: console

    $ pip install easyfermi --upgrade --no-deps
    
If instead, you have **easyfermi 1.X.X** installed, please install easyfermi V2 following section :ref:`installation`.

Uninstalling
----------------

In the terminal, run:

.. code-block:: console

    $ mamba deactivate
    $ mamba env remove --name easyfermi


Running
----------------

If you defined the variable *easyfermi* in your shell environment (see Installing), simply type the following in the terminal:

.. code-block:: console

    $ easyfermi
    
Otherwise, type:

.. code-block:: console

    $ mamba activate easyfermi
    $ python -c "import easyfermi"
    
Substituting mamba by conda if this is the case for you.

