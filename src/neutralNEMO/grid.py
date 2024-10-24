import xarray as xr

def build_nemo_hgrid( hgriddata, iperio = True, jperio = False, gridtype="rectilinear" ):
    """
    Create NeutralOcean grid object from a dictionary of horizontal grid variables.

    Parameters
    __________

    hgriddata: Dictionary of DataArrays (xarray) 
        Describes the horizontal discretization.
        The dictionary can be automatically generated from standard NEMO outputs using
        grid.load_hgriddata or created by hand. Essential entries in the dictionary are
        the following:

               hgriddata["e1u"] -> U-cell widths in i direction
               hgriddata["e1v"] -> U-cell widths in j direction
               hgriddata["e2u"] -> V-cell widths in i direction
               hgriddata["e2v"] -> V-cell widths in j direction

    iperio: Logical (optional, default=True)
        Specifies periodicity in i-direction
        iperio = True -> i-coordinate is periodic (true on global ORCA grid)

    jperio: Logical (optional, default=False)
        Specifies periodicity in j-direction
        jperio = True -> j-coordinate is periodic (false on global ORCA grid)

        
    Returns
    __________

    neutralgrid: Dictionary
        (Description below from NeutralOcean documentation)
        Containing the following:

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

    """
        


    e1u = hgriddata["e1u"]
    e1v = hgriddata["e1v"]
    e2u = hgriddata["e2u"]
    e2v = hgriddata["e2v"]

    shape = ( e1u.shape[-2], e1u.shape[-1] )

    if gridtype.lower() == 'rectilinear':
        from neutralocean.grid.rectilinear import build_grid
        neutralgrid = build_grid( shape, (jperio,iperio), dyC=e1u.to_numpy(), dxC=e2v.to_numpy(), 
                                  dyG=e1v.to_numpy(), dxG=e2u.to_numpy()  )

    elif gridtype.lower() == 'orca':
        from neutralocean.grid.tripolar import build_grid
        neutralgrid = build_grid( shape, e1u=e1u.to_numpy(), e2u=e2u.to_numpy(),
                                  e1v=e1v.to_numpy(), e2v=e2v.to_numpy())

    else:
        raise ValueError("{} is not a supported grid type".format(gridtype))
        
    return neutralgrid

def load_hgriddata( path, e1u_varname = "e1u", e2u_varname = "e2u", e1v_varname = "e1v",
                    e2v_varname = "e2v", i_dim="x", j_dim="y",
                    open_mf=False, **kwargs):
    """
    Create hgriddata dictionary by loading from a standard NEMO output (either a 
    mesh_mask or domaincfg file) or a netcdf file containing the necessary C-grid
    information

    Parameters
    __________

    path: String or List of Strings
        Location(s) of the file(s) to read using xarray.open_dataset (default, open_mf=False) or
        xarray.open_mfdataset (open_mf=True)

        If open_mf=False: A string in the form `"/path/to/my/files/my_file.nc"`
        If open_mf=True: A string glob in the form `"/path/to/my/files/*.nc"`
                         or an explicit list of files to open

    e1u_varname: String (optional, default = "e1u")
        Netcdf variable name for the U-cell width in the i-direction. 
        Can be found using `ncdump -h your_file.nc`

    e2u_varname: String (optional, default = "e2u")
        Netcdf variable name for the U-cell width in the j-direction. 
        Can be found using `ncdump -h your_file.nc`

    e1v_varname: String (optional, default = "e1v")
        Netcdf variable name for the V-cell width in the i-direction. 
        Can be found using `ncdump -h your_file.nc`

    e2v_varname: String (optional, default = "e2v")
        Netcdf variable name for the V-cell width in the j-direction. 
        Can be found using `ncdump -h your_file.nc`

    i_dim: String (optional, default = "x")
        Dimension name for the i axis. 
        Can be found using `ncdump -h your_file.nc`

    j_dim: String (optional, default = "y")
        Dimension name for the j axis. 
        Can be found using `ncdump -h your_file.nc`

    open_mf: Logical (optional, default = False)
        open_mf = False -> Load a dataset from a single netcdf file 
                           (see xarray.open_dataset)
        open_mf = True -> Load a dataset from multiple netcdf files 
                          (see xarray.open_mfdataset)


    kwargs: Keyword arguments
        Keyword arguments for xarray.open_dataset (default, open_mf=False) 
        or xarray.open_mfdataset (default, open_mf=True)

        
    Returns
    __________

    hgriddata: Dictionary
        Contains the following information about the horizontal C-grid 
        discretization

        hgriddata["e1u"]: DataArray
        U-cell width in the i-direction

        hgriddata["e1v"]: DataArray
        V-cell width in the i-direction

        hgriddata["e2u"]: DataArray
        U-cell width in the j-direction

        hgriddata["e2v"]: DataArray
        V-cell width in the j-direction

        If any of the above variables cannot be found in the netcdf file(s) then a warning
        message will appear and the variable will be missing in hgriddata.
    
    """


    if open_mf == False:
        dataset = xr.open_dataset(path, **kwargs)
    else:
        print("Opening multiple files as a single dataset (open_mf = {})".format(open_mf))
        dataset = xr.open_mfdataset(path, **kwargs)

    if i_dim != "x": dataset = dataset.rename({i_dim:"x"})
    if j_dim != "y": dataset = dataset.rename({j_dim:"y"})

    hgriddata = {}

    try:
        hgriddata["e1u"] = dataset[e1u_varname].squeeze()
    except:
        print("Variable {} not found in: {}".format(e1u_varname, path))
    
    try:
        hgriddata["e2u"] = dataset[e2u_varname].squeeze()
    except:
        print("Variable {} not found in: {}".format(e2u_varname, path))

    try:
        hgriddata["e1v"] = dataset[e1v_varname].squeeze()
    except:
        print("Variable {} not found in: {}".format(e1v_varname, path))
    
    try:
        hgriddata["e2v"] = dataset[e2v_varname].squeeze()
    except:
        print("Variable {} not found in: {}".format(e2v_varname, path))
    
    return hgriddata



def load_zgriddata( path, deptht_varname="gdept_0", tmask3d_varname="tmask", tmask2d_varname="tmaskutil", 
                    vert_dim="nav_lev", i_dim="x", j_dim="y", open_mf=False, 
                    infer_tmask2d=False, blev_varname="bottom_level", infer_tmask3d=False, 
                    infer_path='nemo-T.nc', infer_varname='so', 
                    infer_tdim="time_counter", infer_idim="x", infer_jdim="y", infer_zdim="deptht",
                    infer_val=None, calc_deptht=False, e3w_varname='e3w_0', **kwargs ):
    """
    Create zgriddata dictionary by loading from a standard NEMO output (either a 
    mesh_mask or domaincfg file) or a netcdf file containing the necessary C-grid
    information

    Parameters
    __________

    path: String or List of Strings
        Location(s) of the file(s) to read using xarray.open_dataset (default, open_mf=False) or
        xarray.open_mfdataset (open_mf=True)

        If open_mf=False: A string in the form `"/path/to/my/files/my_file.nc"`
        If open_mf=True: A string glob in the form `"/path/to/my/files/*.nc"`
                         or an explicit list of files to open

    deptht_varname: String (optional, default = "gdept_0")
        Netcdf variable name for the T-cell depth. 
        Can be found using `ncdump -h your_file.nc`

    tmask3d_varname: String (optional, default = "tmask")
        Netcdf variable name for the 3-dimensional (z,y,x) tmask. The mask must use 
        the NEMO masking convention, 0 (or False) = Masked.
        Can be found using `ncdump -h your_file.nc`

    tmask2d_varname: String (optional, default = "tmaskutil")
        Netcdf variable name for the 2-dimensional (y,x) tmask. The mask must use 
        the NEMO masking convention, 0 (or False) = Masked.
        Can be found using `ncdump -h your_file.nc`

    vert_dim: String (optional, default = "nav_lev")
        Dimension name for the vertical (k) axis. 
        Can be found using `ncdump -h your_file.nc`
    
    i_dim: String (optional, default = "x")
        Dimension name for the i axis. 
        Can be found using `ncdump -h your_file.nc`

    j_dim: String (optional, default = "y")
        Dimension name for the j axis. 
        Can be found using `ncdump -h your_file.nc`

    calc_deptht: Logical (optional, default = False)
        calc_deptht = False -> Load deptht from the file using the variable name
                               defined by deptht_varname
        calc_deptht = True -> Calculate deptht using cell thicknesses in the file

    e3w_varname: String (optional, default = 'e3w_0')
        Only used if calc_deptht == True
        Variable name for W cell thickness in the file.

    infer_tmask2d: Logical (optional, default = False)
        infer_tmask2d = False -> Load the two-dimensional T mask from the file
                                 using the variable name defined by tmask2d_varname
        infer_tmask2d = True  -> Deduce the two-dimensional T mask from the bottom 
                                 level variable in the file. The bottom level variable
                                 name is defined by blev_varname

    blev_varname: String (optional, default = 'bottom_level')
        Only used if infer_tmask2d == True
        Variable name for bottom level in the file

    infer_tmask3d: Logical (optional, default = False)
        infer_tmask3d = False -> Load the three-dimensional T mask from the file
                                 using the variable name defined by tmask3d_varname
        infer_tmask3d = True  -> Deduced the three-dimensional T mask from a file
                                 defined by infer_path (infer_path can equal path), using
                                 the variable name defined by infer_varname. The dimension
                                 names for the variable are defined by infer_tdim, infer_idim,
                                 infer_jdim, infer_zdim (details below). The method for deducing
                                 the mask depends on the argument infer_val (see below).

    infer_val: String or Float (optional, default = None)
        Only used if infer_tmask3d == True
        infer_val = None  -> The three-dimensional T mask is deduced from the infer variable by
                             identifying null values (e.g. Nan).
        infer_val = Float -> The three-dimensional T mask is deduced from the infer variable by
                             finding occurences of infer_val. Where the infer variable == infer_val,
                             tmask3d = False.

    infer_path: String (optional, default = "nemo-T.nc")
        Only used if infer_tmask3d == True
        Location of the file to read when deducing the three-dimensional T mask.

    infer_varname: String (optional, default = "so")
        Only used if infer_tmask3d == True
        Name of the variable to use when deducing the three-dimensional T mask.

    infer_tdim: String (optional, default = "time_counter")
        Only used if infer_tmask3d == True
        Dimension name for the time axis in the infer variable file (defined by infer_path)
        Can be found using `ncdump -h your_file.nc`

    infer_idim: String (optional, default = "x")
        Only used if infer_tmask3d == True
        Dimension name for the i axis in the infer variable file (defined by infer_path)
        Can be found using `ncdump -h your_file.nc`

    infer_jdim: String (optional, default = "y")
        Only used if infer_tmask3d == True
        Dimension name for the i axis in the infer variable file (defined by infer_path)
        Can be found using `ncdump -h your_file.nc`

    infer_zdim: String (optional, default = "deptht")
        Only used if infer_tmask3d == True
        Dimension name for the z axis in the infer variable file (defined by infer_path)
        Can be found using `ncdump -h your_file.nc`                          

    open_mf: Logical (optional, default = False)
        open_mf = False -> Load a dataset from a single netcdf file 
                           (see xarray.open_dataset)
        open_mf = True -> Load a dataset from multiple netcdf files 
                          (see xarray.open_mfdataset)        

    kwargs: Keyword arguments
        Keyword arguments for xarray.open_dataset (default, open_mf=False) 
        or xarray.open_mfdataset (default, open_mf=True)


    Returns
    __________

    zgriddata: Dictionary
        Contains the following information about the vertical C-grid 
        discretization

        hgriddata["deptht"]: DataArray (z_t, y, x)
        Depth of the T-points 

        hgriddata["tmask3d"]: DataArray (z_t, y, x)
        3D mask for T-points using the NEMO masking convention (False = masked point)

        hgriddata["tmask2d"]: DataArray (y, x)
        2D mask for T-points using the NEMO masking convention (False = masked point)

        If any of the above variables cannot be found in the netcdf file(s) then a warning
        message will appear and the variable will be missing in hgriddata.
    

    """

    if open_mf == False:
        dataset = xr.open_dataset(path, **kwargs)
    else:
        print("Opening multiple files as a single dataset (open_mf = {})".format(open_mf))
        dataset = xr.open_mfdataset(path, **kwargs)

    #Rename vertical coordinate for consistency between datasets
    if vert_dim != "z_t": dataset = dataset.rename({vert_dim:"z_t"}) 
    if i_dim != "x": dataset = dataset.rename({i_dim:"x"})
    if j_dim != "y": dataset = dataset.rename({j_dim:"y"})

    if infer_tmask3d == True:
        infer_dataset = xr.open_dataset(infer_path,**kwargs)
        if infer_zdim != "z_t": infer_dataset = infer_dataset.rename({infer_zdim:"z_t"}) 
        if infer_idim != "x": infer_dataset = infer_dataset.rename({infer_idim:"x"})
        if infer_jdim != "y": infer_dataset = infer_dataset.rename({infer_jdim:"y"})
    
    zgriddata = {}

    if calc_deptht == False:
        try:
            zgriddata["deptht"] = dataset[deptht_varname].squeeze()
        except:
            print("Variable {} not found in: {}".format(deptht_varname, path))
            print("Try running with calc_deptht = True")
            print('Or add the depth manually using the key "deptht"')
    else:
        e3w = dataset[e3w_varname]
        zgriddata["deptht"] = e3w.cumsum(dim='z_t').squeeze()

    if infer_tmask3d == False:
        print("Loading masks directly from files")
        try:
            zgriddata["tmask3d"] = dataset[tmask3d_varname].squeeze().astype(bool)
        except:
            print("Variable {} not found in: {}".format(tmask3d_varname, path))
            print("Try running with infer_tmask3d = True")
            print('Or add the mask manually using the key "tmask3d"')

    else:
        print("Inferring tmask3d from other data field [Variable {} in File {}]".format(infer_varname, infer_path))
        infer_var = infer_dataset[infer_varname]
        if infer_tdim in infer_var.dims: infer_var = infer_var.isel({infer_tdim:-1})
        
        if infer_val is None:
            zgriddata["tmask3d"] = infer_var.notnull()
        else:
            zgriddata["tmask3d"] = xr.where(infer_var == infer_val, False, True).astype(bool)

    if infer_tmask2d == False:
        try:
            zgriddata["tmask2d"] = dataset[tmask2d_varname].squeeze().astype(bool)
        except:
            print("Variable {} not found in: {}".format(tmask2d_varname, path))
            print("Try running with infer_tmask2d = True")
            print('Or add the tmask2d manually using the key "tmask2d"')

    else:
        print("Inferring tmask2d from bottom level data [Variable {} in file {}]".format(blev_varname, path))
        blev = dataset[blev_varname].squeeze().astype(int)
        zgriddata["tmask2d"] = xr.where(blev == 0, False, True).astype(bool)

    

    return zgriddata

