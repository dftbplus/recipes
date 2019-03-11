cubefile = "wp-1-1-743-real.cube"
isovalue = 0.02
load @cubefile
isosurface plus @isovalue @cubefile
isosurface fill translucent
isosurface minus @{isovalue * -1} @cubefile
isosurface color red
isosurface fill translucent 

