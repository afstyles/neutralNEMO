ORCA example
=====

.. _orcaexample

Example files
--------------

The ORCA tripolar grid is similar to a latitude-longitude grid in places but splits the north pole 
singularity into two pieces which reside on land. Because of this, the ORCA grid has to be treated 
differently to other grids.

Example ORCA1 (1 degree equivalent) data can be downloaded from the following archive

`Link to Zenodo archive <https://doi.org/10.5281/zenodo.10639898>`_

You will need the following files

* ``domcfg_eORCA1v2.2x.nc`` - Grid information
* ``nemo_bu978o_1y_19771201-19781201_grid-T.nc`` - T point data containing temperature and salinity information


Loading grid data
--------------

First, we will use the ``load_hgridata`` and ``load_zgriddata`` routines to load the necessary
grid data from the model grid file (typically a ``mesh_mask`` or ``domaincfg`` file).

Loading the horizontal grid data is simple enough

.. code-block:: Python

   from neutralNEMO.grid import load_hgriddata

   hgd = load_hgriddata( "domcfg_eORCA1v2.2x.nc"  )

``hgd`` is a dictionary containing cell widths for the grid. The netcdf variable name for all of these terms
can vary between datasets and can be adjusted accordingly using keyword arguments. For example:

.. code-block:: Python

   hgd = load_hgriddata( "domcfg_eORCA1v2.2x.nc"  , e1u_varname="e1u_0")

When loading the vertical grid data, we load cell thicknesses alongside the two-dimensional (y,x) and three-dimensional 
(z,y,x) T point mask. Often the masks are not provided explicitly in the grid data but are applied to the data. This is the
case for the example dataset ``domcfg_eORCA1v2.2x.nc``, so additional keyword arguments are needed. For example:

.. code-block:: Python

   from neutralNEMO.grid import load_zgriddata

   zgd = load_zgriddata( "domcfg_eORCA1v2.2x.nc", infer_tmask2d=True, infer_tmask3d=True, 
                          infer_path="nemo_bu978o_1y_19771201-19781201_grid-T.nc", 
                          infer_varname="so", vert_dim="z" )

``zgd`` is a dictionary containing depth information, the two-dimensional T mask, and the three-dimensional T mask.

Because ``infer_mask2d = True``, the two-dimensional T mask is inferred from the bottom cell information in 
``domcfg_eORCA1v2.2x.nc``.

Because ``infer_mask3d = True``, the three-dimensional T mask is inferred from masked data in another file 
(nemo_bu978o_1y_19771201-19781201_grid-T.nc  defined by ``infer_path``). The three-dimensional variable 'so'
(defined by ``infer_varname``) is used to deduce the three-dimensional mask.

Variable names, coordinate names, and masking conventions can vary between NEMO datasets. To help with this, there
are plenty of options to correctly load the grid data in the :doc:`api`.

Create a grid object
--------------

We then need to create a ``neutralocean`` grid object. This can be simply done using the ``build_nemo_hgrid``
routine

.. code-block:: Python
   
   from neutralNEMO.grid import build_nemo_hgrid

   neutral_grid = build_nemo_hgrid(hgd, iperio=True, jperio=False, gridtype="orca")

The last step before calculating our neutral surface is the loading of the temperature and salinity data.

Load T-S data
--------------

.. code-block:: Python
   
   from neutralNEMO.surf import load_tsdata

   tsd = load_tsdata("nemo_bu978o_1y_19771201-19781201_grid-T.nc", zgd, to_varname="thetao")

As shown above, netcdf variable names can be specified as a keyword argument.

Define the equation of state
--------------

We also need to define the equation of state used in the NEMO simulation. The default equations 
of state used in NEMO are described in ``neutralNEMO.eos`` and can be called here. In this case, the
EOS80 equation of state is used.


.. code-block:: Python

   from neutralNEMO.eos import NEMO_eos

   eos = NEMO_eos( 'eos80' )

If we wished to use the TEOS10 equation of state (and it matched the NEMO simulation) we would simply 
do the following instead

.. code-block:: Python

   eos = NEMO_eos( 'teos10' )

We also have the option to use the configurable simplified equation of state (S-EOS). An example of this
can be found in the :doc:`rectexample`.

If you have used a non-standard equation of state then you will need to define the equation and its gradients 
with respect to T and S. Details on how to do this can be found in the `neutralocean documentation <https://neutralocean.readthedocs.io/en/latest/internals.html#equation-of-state>`_ .


Calculate surfaces
--------------

Now we have our temperature and grid information we can calculate our neutral surfaces. 

.. code-block:: Python

   from neutralNEMO.surf import find_omega_surfs

   zpins = [150., 300.]   # List of depths to pin each surface to
   ipins = [10, 10]       # List of i-indices to pin each surface to
   jpins = [9, 9]         # List of j-indices to pin each surface to
   tpins = [-1,-1]        # List of time indices to pin each surface to

   surf_dataset = find_omega_surfs( tsd, neutral_grid , zgd, zpins, ipins, jpins, tpins,
                                         eos=eos, ITER_MAX=10)

   #Save as netcdf (optional)
   surf_dataset.to_netcdf("my_surfs.nc")

In the above example, two neutral surfaces are calculated. The first surface is pinned to 150 m depth at (i=10,
j=9) in the final time step. The second surface is the same but pinned to 300 m depth. The surfaces depths,
temperatures, and salinities are outputted as an xarray DataSet and can be easilly saved to netcdf.

``ITER_MAX=10`` sets the maximum number of iterations carried out by the neutralocean algorithm.

If you are interested in calculating the associated Veronis density for the surfaces you calculate then ``calc_veronis=True``
enables the calculation of the Veronis density as a label for the density surfaces.

.. code-block:: Python

   ver_ipins = 15
   ver_jpins = 150

   surf_dataset = find_omega_surfs( tsd, neutral_grid , zgd, zpins, ipins, jpins, tpins,
                                         eos=eos, ITER_MAX=10, calc_veronis=True, 
                                         ver_ipins=ver_ipins, ver_jpins=ver_jpins)

If you want the calculated Veronis densities to be consistent between surfaces at a given time step, the Veronis density should be
calculate at the same point for all surfaces. In this case, the Veronis density is calculated at (i=15,j=150).
