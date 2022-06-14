function project_to_ice_edge, lat, lon, ie_lat, ie_lon, stepdist,  mindist = mindist ; stepdist in m

  R = 6356.75*1000.d ; polar earth radius (in m)
  ;R = 6371.; mean earth radius (in m) for checking
  
  ; find points to take bearing between (more than 100km away if possible, if not return the farthest dist) 
;  distgt100km = where(distance gt 100000., /null)
;  if n_elements(distgt100km) ge 1. then dist2ind = distgt100km[0] else 
;  
;  
    distance = ll2dist(lat[0], lon[0], lat, lon)*1000
    dist2test = where(distance ge 100000, /null) ;  
    if n_elements(dist2test) ge 1 then dist2ind = dist2test[0] else dist2ind = float(n_elements(distance))-1.
    total_dist = distance[dist2ind] ; in m 
  

;dist2ind = float(n_elements(distance))-1.
  
  lat1 = lat[dist2ind] 
  lon1 = lon[dist2ind]
  lat2 = lat[0] 
  lon2 = lon[0]
  
;    lat1 = -64.9 ; some parameters used for testing
;    lon1 = -132.1
;    lat2 = -65.2
;    lon2 = -132.2
    
    lat1rad = lat1*!DTOR
    lon1rad = lon1*!DTOR
    lat2rad = lat2*!DTOR
    lon2rad = lon2*!DTOR
  
  ;step_dist_start = 6250.d
  step_dist_start = double(stepdist)
  step_dist = step_dist_start
  
  delta = total_dist/R
  max_lat = max(ie_lat)
  
  proj_coords_all = []
  mindist = []
  
  counter = 0 

  repeat begin
    f = 1.+step_dist/total_dist
    a = sin((1-f)*delta)/sin(delta)
    b = sin(f*delta)/sin(delta)
    x = a*cos(lat1rad)*cos(lon1rad)+b*cos(lat2rad)*cos(lon2rad)
    y = a*cos(lat1rad)*sin(lon1rad)+b*cos(lat2rad)*sin(lon2rad)
    z = a*sin(lat1rad)+b*sin(lat2rad)
    outlatrad = atan(z, sqrt(x^2+y^2))
    outlonrad = atan(y,x)
    outlat = outlatrad*!RADEG
    outlon = outlonrad*!RADEG
    proj_coords_all = [[proj_coords_all],[outlat,outlon]]
    step_dist = step_dist + step_dist_start 
    mindist = [mindist, min(ll2dist(outlat,outlon,ie_lat,ie_lon)*1000.)]
    counter = counter + 1
    
  endrep until outlat gt max_lat

return, proj_coords_all

end
