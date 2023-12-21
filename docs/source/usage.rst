Usage
=====

.. _installation:

Installation
------------

neutralNEMO can be installed using conda

.. code-block:: console

   (.venv) $ conda install -c conda-forge neutralNEMO

alternatively, neutralNEMO can be installed through pip

.. code-block:: console

   (.venv) $ pip install neutralNEMO


Example of use 
--------------

First, we will use the ``load_hgridata`` and ``load_zgriddata`` routines to load the necessary
data from the model grid file (typically a ``mesh_mask`` or ``domaincfg`` file)

.. code-block:: Python

   from neutralNEMO.grid import load_hgriddata, load_zgriddata

   hgd = load_hgriddata( "/path/to/file/mesh_mask.nc" )
   zgd = load_zgriddata( "/path/to/file/mesh_mask.nc" )

``hgd`` and ``zgd`` are dictionaries containing cell widths, cell thicknesses and T-point masks
contained within ``mesh_mask.nc``. The netcdf variable names for these terms can vary
between datasets and can be adjusted accordinly using keyword arguments. For example:

.. code-block:: Python

   hgd = load_hgriddata( "/path/to/file/mesh_mask.nc" , e1u_varname="e1u_0")

A complete list of the available keyword arguments can be found in the :doc:`api`.

We then need to create a ``neutralocean`` grid object. This can be simply done using the ``build_nemo_hgrid``
routine

.. code-block:: Python
   
   from neutralNEMO.grid import build_nemo_hgrid

   neutral_grid = build_nemo_hgrid( hgd, iperio=False, jperio = False)

The ``iperio`` and ``jperio`` keyword arguments specify the periodicity of the model domain. For example, if
``iperio = True`` then the model is periodic in the i-direction (this is the case for a global configuration).

The last step before calculating our neutral surface is the loading of the temperature and salinity data.

.. code-block:: Python
   
   from neutralNEMO.surf import load_tsdata

   tsd = load_tsdata( "/path/to/file/NEMO_grid-T.nc", to_varname="temperature" )

As shown above, the specific netcdf variable name can be specified as a keyword argument.

Now we have our temperature and grid information we can finally calculate our neutral surface. 

.. code-block:: Python

   from neutralNEMO.surf import find_omega_surfs

   zpins = [150., 300.]   # List of depths to pin each surface to
   ipins = [10, 10]       # List of i-indices to pin each surface to
   jpins = [9, 9]         # List of j-indices to pin each surface to
   tpins = [-1,-1]        # List of time indices to pin each surface to

   surf_dataset = find_omega_surfs( tsd, neutral_grid , zgd, [150., 300.], [10,10], [9,9], [-1,-1],
                                         eos="gsw", ITER_MAX=10, calc_veronis=True)

   #Save as netcdf (optional)
   surf_dataset.to_netcdf("my_surfs.nc")

In the above example, two neutral surfaces are calculated. The first surface is pinned to 150 m depth at (i=9,
j=9) in the final time step. The second surface is the same but pinned to 300 m depth. The surfaces depths,
temperatures, and salinities are outputted as an xarray DataSet and can be easilly saved to netcdf.

To calculate the initial potential density, the equation of state needs to be known. In this case, the ``gsw`` 
equation of state is adopted (see neutralocean documentation for specifics on the equation of state.)

``ITER_MAX=10`` sets the maximum number of iterations carried out by the neutralocean algorithm.

``calc_veronis=True``, enables the calculation of the Veronis density as a label for the density surfaces.

