import xarray as xr
import copy
import numpy as np
from neutralocean.label import veronis_density
from neutralocean.surface.isopycnal import potential_surf
from neutralocean.surface import omega_surf
from neutralocean.traj import neutral_trajectory

def load_tsdata( path, zgriddata, to_varname="to", so_varname="so", vert_dim="deptht",
                 i_dim="x", j_dim="y", open_mf=False, **kwargs):
    
    """
    Create tsdata dictionary by loading from a standard NEMO output (typical a *grid-T.nc file)
    or a netcdf file containing the necessary temperature and salinity information.

    Parameters
    __________

    path: String or List of Strings
        Location(s) of the file(s) to read using xarray.open_dataset (default, open_mf=False) or
        xarray.open_mfdataset (open_mf=True)

        If open_mf=False: A string in the form `"/path/to/my/files/my_file.nc"`
        If open_mf=True: A string glob in the form `"/path/to/my/files/*.nc"`
                         or an explicit list of files to open

    zgriddata: Dictionary
        Must contain the following information about the vertical C-grid 
        discretization

        hgriddata["tmask3d"]: DataArray (z_t, y, x)
        3D mask for T-points using the NEMO masking convention (False = masked point)

        
    to_varname: String (optional, default = "to")
        Netcdf variable name for the potential or conservative temperature 
        Can be found using `ncdump -h your_file.nc`

    so_varname: String (optional, default = "so")
        Netcdf variable name for the practical or absolute salinity
        Can be found using `ncdump -h your_file.nc`
        
    vert_dim: String (optional, default = "deptht")
        Dimension name for the vertical (k) axis. 
        Can be found using `ncdump -h your_file.nc`
    
    i_dim: String (optional, default = "x")
        Dimension name for the i axis. 
        Can be found using `ncdump -h your_file.nc`

    j_dim: String (optional, default = "y")
        Dimension name for the j axis. 
        Can be found using `ncdump -h your_file.nc`

        
    Returns
    __________

    tsdata: Dictionary
        Contains the following information about the ocean temperature and salinity

        tsdata["to"]: DataArray (time, z_t, y, x)
        Potential or conservative temperature 

        tsdata["so"]: DataArray (time, z_t, y, x)
        Practical or absolute salinity 


    """

    if open_mf == False:
        dataset = xr.open_dataset(path, **kwargs)
    else:
        print("Opening multiple files as a single dataset (open_mf = {})".format(open_mf))
        dataset = xr.open_mfdataset(path, **kwargs)

    #Rename dimensions for consistency between ts, hgrid, and zgrid datasets
    if vert_dim != "z_t": dataset = dataset.rename({vert_dim:"z_t"}) 
    if i_dim != "x": dataset = dataset.rename({i_dim:"x"})
    if j_dim != "y": dataset = dataset.rename({j_dim:"y"})

    tsdata = {}

    try:
        tsdata["to"] = dataset[to_varname].where(zgriddata["tmask3d"],np.nan)
    except:
        print("Variable {} not found in: {}".format(to_varname, path))

    try:
        tsdata["so"] = dataset[so_varname].where(zgriddata["tmask3d"],np.nan)
    except:
        print("Variable {} not found in: {}".format(to_varname, path))

    return tsdata
    

def find_omega_surfs( tsdata, neutralgrid, zgriddata, 
                     zpins, ipins, jpins, tpins, eos="gsw" , eos_type="specvol",
                     ITER_MAX=10, 
                     calc_veronis=True, calc_potsurf=True,
                     ver_ipins=None, ver_jpins=None,
                     ver_ref=0., pot_ref=0., **kwargs):
    
    """
    Calculate omega surfaces with a fixed point at a specific locations in time and 
    space. There are options to also calculate the veronis density and/or the
    intersecting potential density surfaces for labelling or reference.

    Parameters
    __________

    tsdata: Dictionary
        Can be generated using `load_tsdata`,
        Must contain the following information about the ocean temperature and salinity.

            tsdata["to"]: DataArray (time, z_t, y, x)
            Potential or conservative temperature (depending on the equation of state, 
            see `eos` argument)

            tsdata["so"]: DataArray (time, z_t, y, x)
            Practical or absolute salinity (depending on the equation of state, 
            see `eos` argument)

    neutralgrid: Dictionary
        (Description below from NeutralOcean documentation)
        Can be generated using `build_nemo_hgrid`
        Must contain the following:

            edges : tuple of length 2
                Each element is an array of int of length `E`, where `E` is the number of
                edges in the grid's graph, i.e. the number of pairs of adjacent water
                columns (including land) in the grid.
                If `edges = (a, b)`, the nodes (water columns) whose linear indices are
                `a[i]` and `b[i]` are adjacent.

            dist : 1d array
                Horizontal distance between adjacent water columns (nodes).
                `dist[i]` is the distance between nodes whose linear indices are
                `edges[0][i]` and `edges[1][i]`.

            distperp : 1d array
                Horizontal distance of the face between adjacent water columns (nodes).
                `distperp[i]` is the distance of the interface between nodes whose
                linear indices are `edges[0][i]` and `edges[1][i]`.


    zgriddata: Dictionary
        Can be generated using `load_zgriddata`
        Must contain the following information about the vertical C-grid 
        discretization

            hgriddata["deptht"]: DataArray (z_t, y, x)
            Depth of the T-points 


    zpins: List of floats
        Depths of the fixed points

    ipins: List of integers
        i-indices of the fixed point

    jpins: List of integers
        j-indices of the fixed point

    tpins: List of integers
        time indices of the fixed point

    eos: String or Function or Tuple of Functions (optional, default="gsw")
        Specification for the equation of state.

        The definitions of temperature (potential vs conservative) and
        salinity (practical vs absolute) must match the equation of state.
        By default the "gsw" equation of state is used which depends on#
        absolute salinity and conservative temperature.

        Two premade options from NeutralOcean are described below:
            "gsw" - TEOS-10 Gibbs Sea Water
            "jmd95" - Jackett and McDougall 1995 equation of state


        (Description below from NeutralOcean documentation)
        If a str, can be any of the strings accepted by
        `neutralocean.eos.tools.make_eos`,
        e.g. `'jmd95'`, `'jmdfwg06'`, `'gsw'`.

        If a function, must take three inputs corresponding to `S`, `T`, and
        `P`, and output the density (or specific volume).  This form is not
        allowed when `diags` is `True`.  This can be made as, e.g.,
        `eos = neutralocean.eos.make_eos('gsw')`
        for a non-Boussinesq ocean, or as
        `eos = neutralocean.eos.make_eos('gsw', grav, rho_c)`
        for a Boussinesq ocean with `grav` and `rho_c` (see inputs below).

        If a tuple of functions, the first element must be a function for the
        equation of state as above, and the second element must be a function
        taking the same three inputs as above and returning two outputs, namely
        the partial derivatives of the equation of state with respect to `S`
        and `T`.  The second element can be made as, e.g.,
        `eos_s_t = neutralocean.eos.make_eos_s_t('gsw', grav, rho_c)`

        The function (or the first element of the tuple of functions) should be
        `@numba.njit` decorated and need not be vectorized -- it will be called
        many times with scalar inputs.

    EOS_TYPE: String
        String describing the quantity returned by the equation of state. The
        options are

        'specvol' : Specific volume used by gsw
        'insitu' : In-situ density 

    ITER_MAX: Integer (option, default=10)
        Maximum number of iterations for omega surface calculation.

    ITER_MIN : Integer (option, default=1)
        Minimum number of iterations for omega surface calculation.
    
    calc_veronis: Logical (default=True)
        Calculate the Veronis density of each surface at their respective
        fixed location in time and space

    calc_potsurf: Logical (default=True)
        Calculate the potential density surface which intersects the fixed point
        location in space and time. If you set ITER_MAX=0 while calc_potsurf=True,
        only potential density surfaces will be calculated.        

    ver_ipins: None or Integer or list of Integers (default=None)
        i-index (or indices) to calculate the Veronis density (if 
        calc_veronis = True). If a single integer, then the Veronis density 
        is calculate at the same position (ver_ipin, ver_jpin) for all surfaces
        (the preferred method for consisten Veronis density)
        
        If None, then ver_ipins = ipins. 
        
        If a list of integers, then the length of the list must match ipins. 

    ver_jpins: None or Integer or list of Integers (default=None)
        j-index (or indices) to calculate the Veronis density (if 
        calc_veronis = True). If a single integer, then the Veronis density 
        is calculate at the same position (ver_ipin, ver_jpin) for all surfaces
        (the preferred method for consisten Veronis density)
        
        If None, then ver_jpins = jpins. 
        
        If a list of integers, then the length of the list must match jpins. 

    ver_ref:
        Reference depth for the calculation of the Veronis density (if 
        calc_veronis = True)

    pot_ref:
        Reference depth for the calculation of the potential density (if 
        calc_potsurf = True)

    kwargs:
        Additional keyword arguments for `neutralocean.surface.omega_surf`

            
    Returns
    __________

    surf_dataset: Dataset containing all of the surfaces requested.

    """
    
    final_set = []

    for zpin, ipin, jpin, tpin, n in zip(zpins, ipins, jpins, tpins, range(len(zpins))):

        shared_attrs = {"ipin":ipins, "jpin":jpins, "pot_ref": pot_ref, }#"eos":eos}

        if calc_veronis==True:

            if ver_ipins is None: ver_ipin = ipin   #If None, use ipin
            elif type(ver_ipins) is list:  ver_ipin = ver_ipins[n] #If list, use nth element
            else: ver_ipin = ver_ipins #If integer, use same point every time

            if ver_jpins is None: ver_jpin = jpin                  #If None, use jpin
            elif type(ver_jpins) is list:  ver_jpin = ver_jpins[n] #If list, use nth element
            else: ver_jpin = ver_jpins                             #If integer, use same point every time

            if ver_ipins is None: shared_attrs["ver_ipins"] = ipins
            else: shared_attrs["ver_ipins"] = ver_ipins

            if ver_jpins is None: shared_attrs["ver_jpins"] = jpins
            else: shared_attrs["ver_jpins"] = ver_jpins



            
            rho_ver = veronis_density(tsdata["so"].to_numpy()[tpin,:,ver_jpin,ver_ipin],
                                        tsdata["to"].to_numpy()[tpin,:,ver_jpin,ver_ipin],
                                        zgriddata["deptht"].where(zgriddata["tmask3d"]).to_numpy()[:,ver_jpin,ver_ipin],
                                        float(zpin), eos=eos, p_ref=ver_ref)

            if eos_type == 'specvol': rho_ver = 1/rho_ver


        if calc_potsurf==True:
            s,t,z,d = potential_surf( tsdata["so"][tpin,:,:,:],
                                    tsdata["to"][tpin,:,:],
                                    zgriddata["deptht"][:,:,:],
                                    grid=neutralgrid, vert_dim="z_t",
                                    ref=pot_ref, pin_cast=(jpin,ipin), 
                                    pin_p=float(zpin), eos=eos )
            
            rho_pot = d["isoval"]
            if eos_type == 'specvol' : rho_pot = 1/rho_pot

            #Prepare depth data in DataArray
            surf_z = xr.DataArray(data=z, dims=["y","x"], name="potsurf_z",
                                     attrs={**shared_attrs, "units":"m", 
                                    "long_name":"Depth of initial potential density surface"})

            #Prepare salinity data in DataArray
            surf_s = xr.DataArray(data=s, dims=["y","x"], name="potsurf_s",
                                     attrs={**shared_attrs, "units":tsdata["so"].units, 
                                    "long_name":"Salinity of initial potential density surface"})

            #Prepare temperature data in DataArray
            surf_t = xr.DataArray(data=t, dims=["y","x"], name="potsurf_t",
                                     attrs={**shared_attrs, "units":tsdata["to"].units, 
                                    "long_name":"Temperature of initial potential density surface"})
            

            potsurf_set = xr.merge([surf_z, surf_s, surf_t])
            potsurf_set = potsurf_set.expand_dims(dim={"rho_pot":np.array([rho_pot])})
            if calc_veronis == True: 
                potsurf_set = potsurf_set.assign_coords(rho_ver=("rho_pot", np.array([rho_ver])))

            final_set = final_set + [potsurf_set]

        if ITER_MAX >= 1:
            s,t,z,d = omega_surf(tsdata["so"][tpin,:,:,:],
                                 tsdata["to"][tpin,:,:,:],
                                 zgriddata["deptht"][:,:,:],
                                 neutralgrid, vert_dim="z_t",
                                 ref=pot_ref, pin_cast=(jpin,ipin),
                                 pin_p=float(zpin), eos=eos, ITER_MAX=ITER_MAX,
                                 **kwargs)
            
            #Prepare depth data in DataArray
            surf_z = xr.DataArray(data=z, dims=["y","x"], name="omsurf_z",
                                     attrs={**shared_attrs, "units":"m", 
                                    "long_name":"Depth of omega density surface"})

            #Prepare salinity data in DataArray
            surf_s = xr.DataArray(data=s, dims=["y","x"], name="omsurf_s",
                                     attrs={**shared_attrs, "units":tsdata["so"].units, 
                                    "long_name":"Salinity of omega density surface"})

            #Prepare temperature data in DataArray
            surf_t = xr.DataArray(data=t, dims=["y","x"], name="omsurf_t",
                                     attrs={**shared_attrs, "units":tsdata["to"].units, 
                                    "long_name":"Temperature of omega density surface"})
            
            omsurf_set = xr.merge([surf_z, surf_s, surf_t])

            if calc_potsurf == True:
                omsurf_set = omsurf_set.expand_dims(dim={"rho_pot":np.array([rho_pot])})
                omsurf_set = omsurf_set.assign_coords(surfno=("rho_pot", np.array([n])))
                if calc_veronis == True:
                    omsurf_set = omsurf_set.assign_coords(rho_ver=("rho_pot", np.array([rho_ver])))

            elif calc_veronis == True:
                omsurf_set = omsurf_set.expand_dims(dim={"rho_ver":np.array([rho_ver])})
                omsurf_set = omsurf_set.assign_coords(surfno=("rho_ver", np.array([n])))
            else:
                omsurf_set = omsurf_set.expand_dims(dim={"surfno":np.array([n])})

            final_set = final_set + [omsurf_set]
            
    

    surf_dataset = xr.merge(final_set)


    return surf_dataset

def find_evolving_omega_surfs( tsdata, neutralgrid, zgriddata, 
                               zpins, ipins, jpins, tpins, calc_veronis=True, calc_potsurf=True, itmin = None, itmax = None,
                               **kwargs ):

    nt = len(tsdata['to']['time_counter'])
    if itmin is None: itmin = 0
    if itmax is None: itmax = nt-1


    # Calculate the pinned surfaces first (t=tpin)
    surf_pin = find_omega_surfs( tsdata, neutralgrid , zgriddata, zpins, ipins, jpins, tpins, 
                                        calc_veronis=calc_veronis, calc_potsurf=calc_potsurf, **kwargs)

    if (calc_veronis == True) or (calc_potsurf == True):
        rho_dimname = [x for x in surf_pin.dims if 'rho_' in x][0]
        surf_pin = surf_pin.swap_dims({rho_dimname:'surfno'})


    n_S = []
    n_T = []
    n_Z = []
    for zpin, ipin, jpin, tpin, n in zip(zpins, ipins, jpins, tpins, range(len(zpins))):

        # Define the Temp, Sal, and Depth (Pressure) cast
        castT = tsdata['to'].isel(x=ipin, y=jpin)
        castS = tsdata['so'].isel(x=ipin, y=jpin)
        castZ = zgriddata['deptht'].where(zgriddata['tmask3d']).isel(x=ipin, y=jpin)

        # Find the neutral trajectory AFTER (and including) the time pin
        castT_after = castT.isel(time_counter=np.arange(tpin,nt,1,dtype=int))
        castS_after = castS.isel(time_counter=np.arange(tpin,nt,1,dtype=int))
        castZ_after = castZ.broadcast_like(castT_after)

        n_after = neutral_trajectory(S=castS_after, T=castT_after, P=castZ_after, p0=float(zpin), vert_dim='z_t')

        # Find the neutral trajectory BEFORE (and including) the time pin
        castT_before = castT.isel(time_counter=np.arange(tpin,-1,-1,dtype=int))
        castS_before = castS.isel(time_counter=np.arange(tpin,-1,-1,dtype=int))
        castZ_before = castZ.broadcast_like(castT_before)

        n_before = neutral_trajectory(S=castS_before, T=castT_before, P=castZ_before, p0=float(zpin), vert_dim='z_t')

        # Combine the before and after trajectories to provide a full time series of T, S, and Z
        # n_S = n_S + [np.append(n_before[0][::-1], n_after[0][1:])]
        # n_T = n_T + [np.append(n_before[1][::-1], n_after[1][1:])]
        n_Z = n_Z + [np.append(n_before[2][::-1], n_after[2][1:])]

    n_Z = np.array(n_Z)

    for it in range(itmin, itmax+1):

        time_value = list(tsdata["to"]["time_counter"].isel(time_counter=[it]).values)
        zpins_it = list(n_Z[:,it])
        tpins_it = [it]*len(zpins_it)

        if it != tpin:
            surf_tmp = find_omega_surfs( tsdata, neutralgrid , zgriddata, zpins_it, ipins, jpins, tpins_it , 
                                            calc_veronis=False, calc_potsurf=calc_potsurf, **kwargs)
            
            if calc_potsurf == True: 
                surf_tmp = surf_tmp.drop_vars(['rho_pot'])
                surf_tmp = surf_tmp.swap_dims({'rho_pot':'surfno'})

        else:
            surf_tmp = copy.deepcopy( surf_pin )

        surf_tmp = surf_tmp.expand_dims( {'time_counter':time_value})

        if it == itmin:
            surf_out = copy.deepcopy(surf_tmp)
        else:
            surf_out  = xr.concat((surf_out, surf_tmp), 'time_counter', compat='override', coords='minimal')


    return surf_out, surf_pin