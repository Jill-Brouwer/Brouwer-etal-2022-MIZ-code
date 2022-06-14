
function prelim_filter_cdaso, nanthresh, filename, sicmap, mizedge_sic, sicfile
  
  ; read in preprocessed sic data eastings and northings
  restore, sicfile
  
  ; set good_start_flag and good_end_flag to be zero initially
  good_start_flag = 0
  good_end_flag = 0

      ; calc SIC ice edge info
      ice_edge = where(mizedge_sic eq 1)
      ie_lat = SIClats[ice_edge]
      ie_lon = SIClons[ice_edge]

        print, filename
        ; read in data H5F format
        file_id=H5F_OPEN(filename)

        ; need to first test whether it is l or right that is the strong beam
        ; 0 backward = strong on left
        ; 1 forward = strong on right
        other_path='orbit_info/sc_orient'
        other_id=H5D_OPEN(file_id, other_path)
        sc_orient=H5D_READ(other_id)
        H5D_CLOSE, other_id
        
        
        ; make sure we catch the condition where it is transitioning and discard/don't analyse these ones 
        ; data dictionary: 'transition when changing between the two orientations, scientific quality may be degraded'
        leftright = 'l'
        if n_elements(sc_orient) eq 1 then begin
          if (sc_orient eq 0) or (sc_orient eq 1) then begin
              if (sc_orient EQ 1) then leftright='r'
    
              s = 'gt2'+leftright ; only use beam 2 data
              filepaths=[s+'/sea_ice_segments/heights/height_segment_height', + $
                s+'/sea_ice_segments/heights/height_segment_length_seg', + $
                s+'/sea_ice_segments/latitude', + $
                s+'/sea_ice_segments/longitude'] ;, + $
              vars=['height', 'seglen', 'la', 'lo']
          
              for i = 0,n_elements(filepaths)-1 do begin
                path = filepaths[i]
                path_id = H5D_OPEN(file_id,path)
                result = execute(vars[i] + '=H5D_READ(path_id)')
                H5D_CLOSE, path_id
              endfor

	h5f_close, file_id  ; !!! CLOSEEEE
         
              ; crop out height outliers (greater than 100m)
              height_masked = height
              height_masked[where(height gt 100.)]=!Values.F_nan
              seglen[where(height gt 100.)]=!Values.F_nan ; take these out of seglen too
      
              ;cut out points where points are more than 200m apart (to prevent interpolat_partion errors)
              distance = (ll2dist(la[0], lo[0], la, lo))*1000
              distdiff = ts_diff(distance, 1)
              height_masked[where(distdiff gt 200.)]=!Values.F_nan
         
              for p = 0, 1 do begin   
                      
                      minlat = where(la eq min(la))
            
                      ind = minlat[0] ; this makes sure there is only one - all the same lat_parts/lon_parts are nans
                      if p eq 0 then inds = [0:ind]
                      if p eq 1 then inds = reverse([ind+1:n_elements(la)-1])
      
                      height_masked_part = height_masked[inds]
                      lat_part = la[inds]
                      lon_part = lo[inds]
                      seglen_part = seglen[inds]
                      
                      ; calculate closest mizedge point for track start and end
                      diststart = ll2dist(lat_part[0], lon_part[0], ie_lat, ie_lon)*1000
                      minstart = min(diststart, location)
                      minstartind = array_indices(diststart, location)
                      minstartcoord = [ie_lat[minstartind], ie_lon[minstartind]]
                      
                      ; if lat is higher in IS-2 than ice edge, calculate number of points from the start of track
                      if lat_part[0] gt minstartcoord[0] then begin
                        print, 'beam 2 IS-2 lat start is greater than ice edge'
                        minstartcoord = [lat_part[0],lon_part[0]]
                      endif
      
                      ; calculat_parte distance from those points for each line (in m!!)
                      distfromstart_ie = ll2dist(minstartcoord[0], minstartcoord[1], lat_part, lon_part)*1000.
                      
                      ; select those points within 500km of ie and corresponding heights
                      start500mind = where(distfromstart_ie le 500000.)       
                      height_masked_500m_start = height_masked[start500mind]
                        ; calc how many nans in these sections
                      startnans = float(n_elements(where(~finite(height_masked_500m_start), /null)))
                      ; assume mean seglength of 17m, there should be 400000/17 points 
                      nonanseglen = seglen_part[where(seglen_part ge 0)]
                      nexp = 500000./mean(seglen_part, /nan)
                      startpoints = float(n_elements(start500mind))-startnans
                      percentstartdata = startpoints/nexp*100.
                      
                      ; add in another one for within the first 100m 
                      start200mind = where(distfromstart_ie le 100000.)       
                      height_masked_200m_start = height_masked[start200mind]
                        ; calc how many nans in these sections
                      startnans2 = float(n_elements(where(~finite(height_masked_200m_start), /null)))
                      ; assume mean seglength of 17m, there should be 400000/17 points 
                      nexp2 = 100000./mean(seglen_part, /nan)
                      startpoints2 = float(n_elements(start200mind))-startnans2
                      percentstartdata2 = startpoints2/nexp2*100.
                      
                      ; set 50% nan threshold and save as good flag             
                         
                      if (p eq 0) and ((percentstartdata ge nanthresh) or (percentstartdata2 ge nanthresh)) then good_start_flag = 1
                      if (p eq 1) and ((percentstartdata ge nanthresh) or (percentstartdata2 ge nanthresh)) then good_end_flag = 1 
                      
                      
;                                      mPolStereo = MAP('Polar Stereographic',  LIMIT=[-50,180,-90,-180])
;                                      ;mPolStereo = MAP('Polar Stereographic',  LIMIT=[-59,110,-61,115], /current)
;                                      mc = MAPCONTINENTS(/FILL_BACKGROUND, FILL_COLOR='gray')
;                      
;                                      ;   c = scatterplot(projected_coords[1,*], projected_coords[0,*], /overplot, symbol = 'dot', sym_size =1, sym_filled = 1, sym_fill_color = 'black')
;                                      map2 = SCATTERPLOT(ie_lon, ie_lat, /OVERPLOT, SYMBOL='dot', SYM_SIZE=1, SYM_FILLED=1, SYM_FILL_COLOR='blue')
;                      
;                                      map3 = SCATTERPLOT(lon_part, lat_part, /OVERPLOT, SYMBOL='dot', SYM_SIZE=1, SYM_FILLED=1, SYM_FILL_COLOR='light grey')
;                                       map6 = scatterplot(minstartcoord[1], minstartcoord[0], /overplot, symbol = 'star', sym_size = 2, sym_filled = 1, sym_fill_color = 'green')
;                       
                          
              endfor
          endif
     endif
         

return, [good_start_flag, good_end_flag]

end
