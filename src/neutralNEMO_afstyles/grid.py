from neutralocean.grid.rectilinear import build_grid
import xarray as xr

def build_nemo_grid( hgriddata, iperio = True, jperio = False, verbose=False):

    e1u = hgriddata["e1u"]
    e1v = hgriddata["e1v"]
    e2u = hgriddata["e2u"]
    e2v = hgriddata["e2v"]

    shape = ( e1u.shape[-2], e1u.shape[-1] )

    grid = build_grid( shape, (jperio,iperio), dyC=e1u.to_numpy(), dxC=e2v.to_numpy(), dyG=e1v.to_numpy(), dxG=e2u.to_numpy()  )

    return grid

def load_hgriddata( path, e1u_varname = "e1u", e2u_varname = "e2u", e1v_varname = "e1v", e2v_varname = "e2v" ):

    dataset = xr.open_dataset(path)

    hgriddata = {}

    try:
        hgriddata["e1u"] = dataset[e1u_varname]
    except:
        print("Variable ", e1u_varname, " not found in: ", path)
    
    try:
        hgriddata["e2u"] = dataset[e2u_varname]
    except:
        print("Variable ", e2u_varname, " not found in: ", path)

    try:
        hgriddata["e1v"] = dataset[e1v_varname]
    except:
        print("Variable ", e1v_varname, " not found in: ", path)
    
    try:
        hgriddata["e2v"] = dataset[e2v_varname]
    except:
        print("Variable ", e2v_varname, " not found in: ", path)
    
    return hgriddata

