Rectilinear example
=====

.. _rectexample

Example files
--------------

On a rectilinear grid, the grid cells are always rectangular and the grid axes are parallel to the
coordinate axes. Rectilinear grids are often used in NEMO test configurations such as GYRE_PISCES.

.. note::

   The ORCA grid is not rectilinear. If you are working with global model data, please refer to the 
   :doc:`orcaexample`.

Example rectilinear data can be downloaded from the following archive

`Link to Zenodo archive <https://doi.org/10.5281/zenodo.10639898>`_

You will need the following files

* ``mesh_mask.nc`` - Grid information
* ``GYRE_1m_04110101_04201230_grid_T.nc`` - T point data containing temperature and salinity information


Loading grid data
--------------

First, we will use the ``load_hgridata`` and ``load_zgriddata`` routines to load the necessary
grid data from the model grid file (typically a ``mesh_mask`` or ``domaincfg`` file).

Loading the horizontal grid data is simple enough

.. code-block:: Python

   from neutralNEMO.grid import load_hgriddata

   hgd = load_hgriddata( "mesh_mask.nc"  )

``hgd`` is a dictionary containing cell widths for the grid. The netcdf variable name for all of these terms
can vary between datasets and can be adjusted accordingly using keyword arguments. For example:

.. code-block:: Python

   hgd = load_hgriddata( "mesh_mask.nc"  , e1u_varname="e1u_0")

When loading the vertical grid data, we load cell thicknesses alongside the two-dimensional (y,x) and three-dimensional 
(z,y,x) T point mask. In this case, the necessary masks are contained in the ``mesh_mask.nc`` file.

.. code-block:: Python

   from neutralNEMO.grid import load_zgriddata

   zgd = load_zgriddata( "mesh_mask.nc", vert_dim="nav_lev" )

``zgd`` is a dictionary containing depth information, the two-dimensional T mask, and the three-dimensional T mask. Note that 
the vertical dimension name of the netcdf file had to be specified, because it was a non-default name.

In some cases, such as the :doc:`orcaexample`, the masks are not provided explicitly in the grid file but are applied to the data.
In this case, the ``infer_mask2d`` and/or ``infer_mask3d`` keyword arguments will be needed.

Variable names, coordinate names, and masking conventions can vary between NEMO datasets. To help with this, there
are plenty of options to correctly load the grid data in the :doc:`api`.

Create a grid object
--------------

We then need to create a ``neutralocean`` grid object. This can be simply done using the ``build_nemo_hgrid``
routine

.. code-block:: Python
   
   from neutralNEMO.grid import build_nemo_hgrid

   neutral_grid = build_nemo_hgrid(hgd, iperio=False, jperio=False, gridtype="rectilinear")

The last step before calculating our neutral surface is the loading of the temperature and salinity data.

Load T-S data
--------------

.. code-block:: Python
   
   from neutralNEMO.surf import load_tsdata

   tsd = load_tsdata("GYRE_1m_04110101_04201230_grid_T.nc", zgd, 
                     to_varname="votemper", so_varname="vosaline")

As shown above, netcdf variable names can be specified as a keyword argument. In this case, ``votemper`` is the variable name
for the potential temperature and ``vosaline`` is the practical salinity.

Define the equation of state
--------------

We also need to define the equation of state used in the NEMO simulation. The default equations 
of state used in NEMO are described in ``neutralNEMO.eos`` and can be called here. In this case, we
configure the simplified equation of state (S-EOS) to describe a linear equation of state using the 
same parameters found in the NEMO namelist.

.. code-block:: Python

   from neutralNEMO.eos import NEMO_eos

   eos = NEMO_eos( 'seos' , rn_a0 = 1.655e-1, rn_b0 = 7.655e-1,
                    rn_nu = 0., rn_lambda1 = 0., rn_lambda2 = 0.,
                    rn_mu1 = 0., rn_mu2 = 0. )

If we wished to use the TEOS10 equation of state (and it matched the NEMO simulation) we would simply 
do the following instead

.. code-block:: Python

   eos = NEMO_eos( 'teos10' )

We also have the option to use EOS80. An example of this
can be found in the :doc:`orcaexample`.

If you have used a non-standard equation of state then you will need to define the equation and its gradients 
with respect to T and S. Details on how to do this can be found in the `neutralocean documentation <https://neutralocean.readthedocs.io/en/latest/internals.html#equation-of-state>`_ .

Calculate surfaces
--------------

Now we have our temperature and grid information we can calculate our neutral surfaces. 

.. code-block:: Python

   from neutralNEMO.surf import find_omega_surfs

   zpins = [100., 200.]   # List of depths to pin each surface to
   ipins = [10, 10]       # List of i-indices to pin each surface to
   jpins = [5, 5]         # List of j-indices to pin each surface to
   tpins = [0, 0]        # List of time indices to pin each surface to

   surf_dataset = find_omega_surfs( tsd, neutral_grid , zgd, zpins, ipins, jpins, tpins,
                                         eos=eos, ITER_MAX=10)

   #Save as netcdf (optional)
   surf_dataset.to_netcdf("my_surfs.nc")

In the above example, two neutral surfaces are calculated. The first surface is pinned to 100 m depth at (i=10,
j=5) in the first time step. The second surface is the same but pinned to 200 m depth. The surfaces depths,
temperatures, and salinities are outputted as an xarray DataSet and can be easilly saved to netcdf.

``ITER_MAX=10`` sets the maximum number of iterations carried out by the neutralocean algorithm.

If you are interested in calculating the associated Veronis density for the surfaces you calculate then ``calc_veronis=True``
enables the calculation of the Veronis density as a label for the density surfaces.

.. code-block:: Python

   ver_ipins = 10
   ver_jpins = 10

   surf_dataset = find_omega_surfs( tsd, neutral_grid , zgd, zpins, ipins, jpins, tpins,
                                         eos=eos, ITER_MAX=10, calc_veronis=True, 
                                         ver_ipins=ver_ipins, ver_jpins=ver_jpins)

If you want the calculated Veronis densities to be consistent between surfaces at a given time step, the Veronis density should be
calculate at the same point for all surfaces. In this case, the Veronis density is calculated at (i=15,j=150).
