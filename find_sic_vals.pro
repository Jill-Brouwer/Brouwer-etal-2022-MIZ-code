function find_sic_vals, lat, lon, SICeasting, SICnorthing, SICmap, eastingmatch = eastingmatch, northingmatch = northingmatch

  all_SIC_vals = fltarr(n_elements(lat), /nozero)
  eastingmatch = fltarr(n_elements(lat), /nozero)
  northingmatch = fltarr(n_elements(lat), /nozero)

  for i = 0, n_elements(lat)-1 do begin
    xy = ps2xy_chadstyle2(lat[i],lon[i],-70.0,0.0)
    eastdiff = abs(SICeasting-xy[0]) ;find nearest pixel on SIC map using SIC easting/northing tables
    northdiff = abs(SICnorthing-xy[1])
    indeast = where(eastdiff eq min(eastdiff))
    indnorth = where(northdiff eq min(northdiff))
    eastingmatch[i] = indeast
    northingmatch[i] = indnorth
    ;retrieve SIC conc for each midpoint
    thisSIC = SICmap[indeast, indnorth]
    all_SIC_vals[i] = thisSIC
  endfor

  return, all_SIC_vals

end