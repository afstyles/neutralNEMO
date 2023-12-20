Welcome to neutralNEMO's documentation!
===================================

``neutralNEMO`` is a Python library which helps users calculate neutral surfaces from NEMO model outputs.

This library makes it easier to use the highly-optimized algorithms of the `neutralocean <https://github.com/geoffstanley/neutralocean>`_ 
package (`Stanley et al., 2021 <https://doi.org/10.1029/2020MS002436>`_) with model outputs from NEMO.

The package contains routines to correctly load temperature and salinity data alongside the necessary horizontal and vertical grid information.

A streamlined routine is also included to calculate an approximately neutral density surface (an omega surface) which is fixed to a specific 
location in time and space. Additional options to calculate the Veronis density and potential density surfaces are also included.

Take a look at the :doc:`usage` section to find out how to get started with neutralNEMO.

.. note::

   This project is under active development. 
   
   neutralNEMO is at an early stage of development and I am keen to here how others get on with the 
   library when they use it on their model runs. If it works well for you, please spread the word and
   let me know as it will make my day! If you do encounter any problems, please contact me or raise
   an issue on GitHub, I will try to help in any way I can.

Contents
--------

.. toctree::

   usage
   api
