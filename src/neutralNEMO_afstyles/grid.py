from neutralocean.grid.rectilinear import build_grid
import xarray as xr

def build_nemo_grid( hgriddata, iperio = True, jperio = False):
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

    grid: Dictionary
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

    grid = build_grid( shape, (jperio,iperio), dyC=e1u.to_numpy(), dxC=e2v.to_numpy(), 
                      dyG=e1v.to_numpy(), dxG=e2u.to_numpy()  )

    return grid

def load_hgriddata( path, e1u_varname = "e1u", e2u_varname = "e2u", e1v_varname = "e1v",
                    e2v_varname = "e2v", open_mf=False, **kwargs):
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

