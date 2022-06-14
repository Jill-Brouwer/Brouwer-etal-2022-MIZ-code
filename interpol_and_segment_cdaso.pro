pro interpol_and_segment_cdaso, filename, good_start_flag, good_end_flag, sicmap, mizedge_sic, sicpath, outpath,plotpath, sec_len, sample_int, nan_thresh, do_plots = do_plots

    restore, sicpath ; ancillary sic grid coord information

    sfi = good_start_flag
    efi = good_end_flag
    print, 'start interpol'+string(filename) + string(sfi) + string(efi)
    split = strsplit(filename, '/')
    split2 = filename.Substring(split[n_elements(split)-1])
    atldate = split2.substring(9,16)
    ID = split2.substring(9,31)

    ; get SIC ice edge data
    ice_edge = where(mizedge_sic eq 1)
    ie_lat = SIClats[ice_edge]
    ie_lon = SIClons[ice_edge]
    
    ; get julian day from atldate for labelling
    atlyear = split2.substring(9,12)
    atlmonth = split2.substring(13,14)
    atlday = split2.substring(15,16)
    atljuliandate = julday(atlmonth,atlday, atlyear)

    ; folder to save plots
    plotsavefolder = plotpath

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
      if sc_orient eq 0 or sc_orient eq 1 then begin
        if (sc_orient EQ 1) then leftright='r'

        s = 'gt'+strtrim(indgen(3)+1,2)+leftright
        m = 'gt'+strtrim(indgen(3)+1,2)
        filepaths=[s+'/sea_ice_segments/heights/height_segment_height', + $
                   s+'/sea_ice_segments/heights/height_segment_length_seg', + $
                   s+'/sea_ice_segments/latitude', + $
                   s+'/sea_ice_segments/longitude', + $
                   s+'/sea_ice_segments/delta_time', + $
                   s+'/sea_ice_segments/heights/height_segment_surface_error_est', + $
                   s+'/sea_ice_segments/heights/height_segment_type'] 
        vars=[m+'height', m+'seglen', m+'lat', m+'lon', + m+'time', +$
          m+'error', m+'segtype']

        for i = 0,n_elements(filepaths)-1 do begin
          path = filepaths[i]
          path_id = H5D_OPEN(file_id,path)
          result = execute(vars[i] + '=H5D_READ(path_id)')
          H5D_CLOSE, path_id
        endfor

        h5f_close, file_id ; need to close the file id too !!
        
    
        for part = 0,1 do begin ; do everything for start and end of track
          ;part = 1
        ;check the good start and end flag and only choose those that meet the criteria (value of 1)
          if part eq 0 and sfi eq 0 then continue
          if part eq 1 and efi eq 0 then continue
          
          if part eq 0 then plab = 'start' else plab = 'end' ; label for plots
    
          for beam = 0,2 do begin ; get variables and sort data individually for each beam
            if beam eq 0 then begin
              height = gt1height & seglen = gt1seglen & lat = gt1lat & lon = gt1lon & error = gt1error & time = gt1time
            endif
            if beam eq 1 then begin
              height = gt2height & seglen = gt2seglen & lat = gt2lat & lon = gt2lon & error = gt2error & time = gt2time
            endif
            if beam eq 2 then begin
              height = gt3height & seglen = gt3seglen & lat = gt3lat & lon = gt3lon & error = gt3error & time = gt3time
            endif
    
            ; crop out height outliers (greater than 100m), and seglengths that match with these (they have really large seglenths and we dont want to include these in the seglen stats)
            height_masked = height & seglen_masked = seglen & error_masked = error
            height_masked[where(height gt 100.)] = !Values.F_nan
            seglen_masked[where(height gt 100.)] = !Values.F_nan
            error_masked[where(height gt 100.)] = !Values.F_nan
    
            ; split the track into start and end based on min lat
            ; add in a clause here based on prelim selection flags in get_seglen
    
            minlat = where(lat eq min(lat))
    
            ind = minlat[0] ; this makes sure there is only one - all the same lats/lons are nans
            if part eq 0 then inds = [0:ind-1] ; did minus one in case lowest lat is on other side of Antarctica
            if part eq 1 then inds = reverse([ind+1:n_elements(lat)-1])
    
            height_masked_part = height_masked[inds]
            lat_part = lat[inds]
            lon_part = lon[inds]
            seglen_part = seglen_masked[inds]
            error_part = error_masked[inds]
            time_part = time[inds] ; save time vals for getting midpoint times
            
            if n_elements(lat_part) le 300. then data_flag = 1 else data_flag = 0 ; 300 points is approximately 1 section length assuming mean seglength of 20m
    
            ; calculate along-track distance
            ; correct for interpol errors: set points where distance to next point gt 200m as nan
            distance = (ll2dist(lat_part[0], lon_part[0], lat_part, lon_part))*1000. ; distance in meters
            distdiff = abs(ts_diff(distance,1))
            height_masked_part[where(distdiff gt 200.)] = !Values.f_nan
    
            ; interpolate (currently using least squares quadratic with 4 point neighbourhood)
            n_interp=(max(distance)-min(distance))/sample_int
            new_x = (dindgen(n_interp)*sample_int)+min(distance)
            height_reg = interpol(height_masked_part, distance, new_x, /spline) ; changed to spline 8 sept
    
    
            ; save variables
            if beam eq 0 then begin
              b1_dist = distance & b1_newx = new_x & b1_height = height_reg
              b1_lat = lat_part & b1_lon = lon_part
              b1_seglen = seglen_part & b1_error = error_part & b1_time = time_part
            endif
            if beam eq 1 then begin
              b2_dist = distance & b2_newx = new_x & b2_height = height_reg
              b2_lat = lat_part & b2_lon = lon_part
              b2_seglen = seglen_part & b2_error = error_part & b2_time = time_part
            endif
            if beam eq 2 then begin
              b3_dist = distance & b3_newx = new_x & b3_height = height_reg
              b3_lat = lat_part & b3_lon = lon_part
              b3_seglen = seglen_part & b3_error = error_part & b3_time = time_part
            endif
    
          endfor
          
          if data_flag eq 1 then continue
          
          print,'finish interpol' 
    
          ; get spatially matching sections based on beam 2
    
           ; find lowest lat start point and closest points for other beams
          startlats = [b1_lat[0],b2_lat[0],b3_lat[0]]
          minlat = min(startlats)
          startbeam = where(startlats eq minlat)+1
          otherbeamind = where(startlats ne minlat)
          endlats = [b1_lat[n_elements(b1_lat)-1], b2_lat[n_elements(b2_lat)-1], b3_lat[n_elements(b3_lat)-1]]
          
          ; add criteria if there are no overlapping 
          if (endlats[otherbeamind[0]] gt minlat) or (endlats[otherbeamind[1]] gt minlat) then continue
    
          ; interpol lat lons to find approximate start point (needs to be within 3.5km )
          ; first need to get coords of interpolated values for the start of the track to compare amongst the beams
          
          b1interpcoord = dblarr(2,100000, /nozero)
          b2interpcoord = dblarr(2,100000, /nozero)
          b3interpcoord = dblarr(2,100000, /nozero)
          
          for beam = 0, 2 do begin
            if beam eq 0 then begin
              lat = b1_lat & lon = b1_lon & dist = b1_dist & new_x = b1_newx
            endif
            if beam eq 1 then begin
              lat = b2_lat & lon = b2_lon & dist = b2_dist & new_x = b2_newx
            endif
            if beam eq 2 then begin
              lat = b3_lat & lon = b3_lon & dist = b3_dist & new_x = b3_newx
            endif
    
            ltind = long(where(lat ge minlat))
            gtind = long(where(lat lt minlat)) ; choose one lower than the minlat
            interpol_inds = [ltind,gtind[0]]
            
            lats_to_interpol = lat[interpol_inds] 
            lons_to_interpol = lon[interpol_inds]
            
            
            xcount = 0
    
            for l = 0, n_elements(interpol_inds)-2 do begin
              lat1 = lats_to_interpol[l] & lat2 = lats_to_interpol[l+1]
              lon1 = lons_to_interpol[l] & lon2 = lons_to_interpol[l+1]
              if lat1 ne lat2 and lon1 ne lon2 then begin ; if they are same coords don't calculate
                startdist = dist[l]
                enddist = dist[l+1]
                totdist = enddist-startdist
                x_ind = where(new_x lt enddist and new_x ge startdist, /null)
                interp_n = n_elements(x_ind)
                if interp_n gt 0 then begin ; if there are no points in between then don't calculate
                  for j = 0, interp_n-1 do begin
                    f = (new_x[x_ind[j]]-startdist)/totdist
                    coord = getlatlon2(lat1, lon1, lat2, lon2, f)              
                    if beam eq 0 then b1interpcoord[*,xcount] = coord
                    if beam eq 1 then b2interpcoord[*,xcount] = coord
                    if beam eq 2 then b3interpcoord[*,xcount] = coord            
                    xcount = xcount + 1
                  endfor
                endif
              endif
            endfor
            
            if beam eq 0 then b1interpcoord = b1interpcoord[*,0:xcount-1]
            if beam eq 1 then b2interpcoord = b2interpcoord[*,0:xcount-1]
            if beam eq 2 then b3interpcoord = b3interpcoord[*,0:xcount-1]
            
          endfor ; end beam loop
    
          ; find closest starting coord to lowest beam start for the other beams using the interpol coords
    
          if startbeam eq 1 then begin
            b1startind = 0
            b1mindist = 0
    
            disttob2 = ll2dist(b1_lat[0], b1_lon[0], b2interpcoord[0,*], b2interpcoord[1,*])
            b2mindist = min(disttob2)
            b2startind = long(where(disttob2 eq b2mindist))
    
            disttob3 = ll2dist(b1_lat[0], b1_lon[0], b3interpcoord[0,*], b3interpcoord[1,*])
            b3mindist = min(disttob3)
            b3startind = long(where(disttob3 eq b3mindist))
          endif
    
          if startbeam eq 2 then begin
            b2startind = 0
            b2mindist = 0
    
            disttob1 = ll2dist(b2_lat[0], b2_lon[0], b1interpcoord[0,*], b1interpcoord[1,*])
            b1mindist = min(disttob1)
            b1startind = long(where(disttob1 eq b1mindist))
    
            disttob3 = ll2dist(b2_lat[0], b2_lon[0], b3interpcoord[0,*], b3interpcoord[1,*])
            b3mindist = min(disttob3)
            b3startind = long(where(disttob3 eq b3mindist))
          endif
    
          if startbeam eq 3 then begin
            b3startind = 0
            b3mindist = 0
    
            disttob1 = ll2dist(b3_lat[0], b3_lon[0], b1interpcoord[0,*], b1interpcoord[1,*])
            b1mindist = min(disttob1)
            b1startind = long(where(disttob1 eq b1mindist))
    
            disttob2 = ll2dist(b3_lat[0], b3_lon[0], b2interpcoord[0,*], b2interpcoord[1,*])
            b2mindist = min(disttob2)
            b2startind = long(where(disttob2 eq b2mindist))
          endif
    
          ; check all are not greater than 3.5km away
          b2starttob1start = ll2dist(b2interpcoord[0,b2startind], b2interpcoord[1,b2startind], b1interpcoord[0,b1startind],b1interpcoord[1,b1startind])
          b2starttob3start = ll2dist(b2interpcoord[0,b2startind], b2interpcoord[1,b2startind], b3interpcoord[0,b3startind],b3interpcoord[1,b3startind])
          bad_start_flag = 0
          if b2starttob1start gt 3.5 or b2starttob3start gt 3.5 then begin
            print, 'min dist of start point greater than 3.5km'
            bad_start_flag = 1
          endif
          
          if bad_start_flag eq 1 then break
          
          n1= float(n_elements(b1_newx))-float(b1startind) ; account for different lengths if start has been cut off
          n2= float(n_elements(b2_newx))-float(b2startind)
          n3=float(n_elements(b3_newx))-float(b3startind)
    
          n_in_data=min([n1,n2,n3])
          kernel_centre=sec_len/2 ; works best for odd-span kernels
          dist_to_first=-kernel_centre ; a negative number
          dist_to_last=sec_len-1-kernel_centre ; a positive number
    
          ; select sections based on threshold - 6.25km sliding window with 1km step
          maxsec = min([max(b1_newx), max(b2_newx), max(b3_newx)])/1000.
          secindall_large = lonarr(sec_len, 3, maxsec)
          secind_large = secindall_large
          secindall_large[*] = !Values.F_nan
          secind_large[*] = !Values.F_nan
    
          ; save perc nans here (mdps in FFT section)
          sec_num = long(0)
          sec_counter = 0l
          goodsec = []
          for i=long(kernel_centre),(n_in_data-1)-dist_to_last do begin
            starti = long(i+dist_to_first) & endi = long(i+dist_to_last) ; if they were all starting at 0
            b1start = starti+b1startind & b1end = endi+b1startind ; need to account for the different start inds
            b2start = starti+b2startind & b2end = endi+b2startind
            b3start = starti+b3startind & b3end = endi+b3startind
            ; get height data
            b1_hdata = b1_height[b1start:b1end]
            b2_hdata = b2_height[b2start:b2end]
            b3_hdata = b3_height[b3start:b3end]
            ; test whether they meet 50% data criteria (either of the two different ones for now)
            nans1 = float(n_elements(where(~finite(b1_hdata), /null)))
            nans2 = float(n_elements(where(~finite(b2_hdata), /null)))
            nans3 = float(n_elements(where(~finite(b3_hdata), /null)))
            tot_perc_nans= (nans1+nans2+nans3)/(sec_len*3)*100.
            b1_perc_nans = nans1/sec_len*100.
            b2_perc_nans = nans2/sec_len*100.
            b3_perc_nans = nans3/sec_len*100.
            ; if acceptable save the heightsecs and other variables
            ;if ((tot_perc_nans ) lt nan_thresh) then begin ; or
            if (b1_perc_nans lt nan_thresh) and (b2_perc_nans lt nan_thresh) and (b3_perc_nans lt nan_thresh) then begin
              sec_num = sec_num+1 
              goodsec = [goodsec, sec_counter]       
            endif   
            secindall_large[*,0,sec_counter] = [b1start:b1end]
            secindall_large[*,1,sec_counter] = [b2start:b2end]
            secindall_large[*,2,sec_counter] = [b3start:b3end]      
            i = i + 124 ; 1km step (294)
            sec_counter = sec_counter+1.
          endfor
          
          secind = secindall_large[*,*,goodsec]
          secindall = secindall_large[*,*,0:sec_counter-1]
          
          secindall_large = 0
          
          if sec_num gt 0 then begin ; if sec num is 0 then it will break (undefined vars)
    
            
            ;for beam = 0, 2 do begin
              beam = 1 ; only do the distance calcs for middle beam (faster)
                 
                if beam eq 0 then begin
                  distance = b1_dist & lat = b1_lat & lon = b1_lon & newx = b1_newx & startind = b1startind 
                endif 
                if beam eq 1 then begin
                  distance = b2_dist & lat = b2_lat & lon = b2_lon & newx = b2_newx & startind = b2startind 
                endif
                if beam eq 2 then begin
                  distance = b3_dist & lat = b3_lat & lon = b3_lon & newx=b3_newx & startind = b3startind 
                endif
                
                ; corrected distance calculation steps are as follows, 
                ;1. interpolate alongtrack coords every 1km
                ;2. extrapolate to ice edge (if required) every 1km
                ;3. combine grids and get SIC for all coords
                ;4. calculate corrected distance for sections by integrating the combined function of SIC and distance
    
                
                interp_num = long(n_elements(secindall[0,0,*]))
                interpcoords = dblarr(2,interp_num, /nozero)
                interpinds = []
                interp_dist = newx[secindall[sec_len/2., beam, *]]
                og_interp_dist = interp_dist
                interp_dist = reform(interp_dist - min(interp_dist)) ; change this distance so everything is relative to first midpoint position
               
                ; get mplat for interpolated sections
                for i = 0, interp_num-1 do begin
                  thisnewx = newx[secindall[sec_len/2.,beam,i]] ; find new x at midpoint
                  distdiff = abs(distance-thisnewx)
                  closest_dist_ind = long(where(distdiff eq min(distdiff)))
                  if n_elements(closest_dist_ind) gt 1 then closest_dist_ind = closest_dist_ind[0]
                  interpinds = [interpinds, closest_dist_ind] ; for querying time later
                  dist1 = distance[closest_dist_ind]
                  n = 1
                  repeat begin
                      if dist1 lt thisnewx then closest_dist_ind2 = closest_dist_ind+n else closest_dist_ind2 = closest_dist_ind-n
                      dist2 = distance[closest_dist_ind2]
                      totdist = abs(dist1-dist2)
                      f = abs(thisnewx - dist1)/totdist
                      lat1 = lat[closest_dist_ind] & lon1 = lon[closest_dist_ind]
                      lat2 = lat[closest_dist_ind2] & lon2 = lon[closest_dist_ind2]
                      n = n + 1
                  endrep until (lat1 ne lat2) and (lon1 ne lon2)
                   
                  thiscoord = getlatlon2(lat1, lon1, lat2, lon2, f)
                  interpcoords[*,i] = thiscoord
                endfor
              
                  ; if further than 6.25km away project coordinates to ice edge (at ~1km interval) and find the closest point 
                  ; projection is from the first interpolated midpoint to ensure a continuous ~1km grid
                  projected_coords = project_to_ice_edge(interpcoords[0,*],interpcoords[1,*],ie_lat,ie_lon, 1000., mindist = thismindist)
                  projlats = reform(projected_coords[0,*])
                  projlons = reform(projected_coords[1,*])
                  
                  projected_coords = 0
                  
                   ; 07/10/2020 no ie needed (no points within 6.25km of an ice edge point)
                    ; need to combine the whole line and proj coords to find our own ice edge
                    ; change coords so going from ice edge in (north to south) 
                    projlatr = reverse(projlats)
                    projlonr = reverse(projlons)        
                    extrapnum = n_elements(projlatr) & interpnum = n_elements(interpcoords[0,*])
                    ;make combined array      
                    all_coords = dblarr(2, extrapnum+interpnum)
                    all_coords[0,0:extrapnum-1] = projlatr
                    all_coords[1,0:extrapnum-1] = projlonr
                    all_coords[*,extrapnum:*] = interpcoords
        
                    ; get sic coords for whole trackline
                    ; find sic vals on the median smoothed map
                    sicmapnonan  = sicmap 
                    sicmapnonan[where(~finite(sicmapnonan))] = 0
                    sicmapsmooth = estimator_filter(SICmapnonan, 9, /median)
                    testsic = find_sic_vals(all_coords[0,*], all_coords[1,*], SICeasting, SICnorthing, SICmapsmooth)                 
                    findice = where(testsic ge 15, /null)
                    sicmapnonan = 0
                    sicmapsmooth = 0
                    
                    if n_elements(findice) ge 1 then begin
                                       
                    icestartind = findice[0]
                                     
                    all_sic_coords = all_coords[*,icestartind:*]
                    ie_start_coord = all_coords[*,icestartind]
                    
                    
                    all_SIC_vals = find_SIC_vals(all_SIC_coords[0,*],all_SIC_coords[1,*], SICeasting, SICnorthing, SICmap, eastingmatch=all_sic_easting, northingmatch=all_SIC_northing)
                    n_proj = n_elements(all_sic_vals) - n_elements(interpcoords[0,*])
                    
                    ; check if snipped to only interpol or if projected are still included to get the indexes
                    if n_proj gt 0 then begin
                      extrapind = n_proj-1
                      sic_extrapol_dist = reform(ll2dist(all_SIC_coords[0,0], all_SIC_coords[1,0], all_SIC_coords[0,0:extrapind], all_SIC_coords[1,0:extrapind])*1000.)
                      min_ie_dist = ll2dist(all_SIC_coords[0,0], all_SIC_coords[1,0], interpcoords[0,0], interpcoords[1,0])*1000.
                      dist_from_ie = interp_dist + min_ie_dist[0]
                      all_SIC_dist = [SIC_extrapol_dist, dist_from_ie] 
                      good_mpinds = goodsec+n_proj
                      extrapol_inds = [0:extrapind]                    
                    endif else begin
                      interpind = (n_elements(interpcoords[0,*]) - n_elements(all_SIC_vals)) ; if you take one away itll be one off not sure why
                      all_SIC_dist = interp_dist[interpind:*] - interp_dist[interpind]
                      ie_start_coord = interpcoords[*,interpind]
                      newsecindall = secindall[*,*,interpind:*]
                      this_good_mpinds = goodsec - interpind
                      newgoodsec = where(this_good_mpinds ge 0)           
                      newsecind = secindall[*,*,goodsec[newgoodsec]]
                      good_mpinds = this_good_mpinds[newgoodsec]     
                      ; overwrite arrays and clear old variables              
                      secind = newsecind
                      secindall = newsecindall
                      sec_num = n_elements(newgoodsec)
                      extrapol_inds = !Values.F_nan
                      newsecind = 0 & newsecindall = 0 & newgoodsec = 0
                      
                    endelse   
                    interpcoords = 0 & all_coords = 0
                                                   
                ;calculate corrected distance from ice edge (Wadhams 1975)
                ;17/10/20 changed to this function instead - the integrate wasnt working well for really low concentrations
                
                
                stop
                
                ; CHANGE THIS TO EXCLUDE LINES WITH NANS OR ONLY INCLUDE UP TO THE NAN _ LAND MASK
                
                all_SIC_vals[where(~finite(all_SIC_vals))] = 100. ; make nans 100% - this happens when there is land, nans do not affect dist calcs
                sicfraction = all_SIC_vals/100.
                all_SIC_corrdist = correct_dist(all_SIC_dist, sicfraction, 1000.)
                
               
;                    all_SIC_corrdist = dblarr(n_elements(all_SIC_dist), /nozero)
;                    for i = 0, n_elements(all_SIC_dist)-1 do begin
;                      thisdist = all_SIC_dist[0:i]
;                      thisSIC = sicfraction[0:i]
;                      if i eq 0 then integrate = 0 else integrate = int_tabulated(thisdist, thisSIC)    ; initial distance should be 0 in all cases
;                      all_SIC_corrdist[i] = integrate
;                    endfor
;                    
     
                ; example of how to get these vals out of big arrays in the sav file:
    ;            mp_dist_from_ie_SICcorr = all_SIC_corrdist[good_mpinds]
    ;            mpSIC = all_SIC_vals[good_mpinds]
                mpinterpeasting = all_SIC_easting[good_mpinds]
                mpinterpnorthing = all_SIC_northing[good_mpinds]
                
    
                lineid = ID+'_'+plab
;    
do_plots=1

                IF n_elements(DO_plots) eq 1 then begin
                  
                    thiswin = window()
                    mpolstereo = Map('Polar Stereographic', LIMIT=[max(b2_lat)+1, max(b2_lon)+1, min(b2_lat)-1, min(b2_lon)-1], /current)
                    ;mPolStereo = MAP('Polar Stereographic',  LIMIT=[-50,180,-90,-180], /current)
                    mc = MAPCONTINENTS(/FILL_BACKGROUND, FILL_COLOR='gray')
                    sp = SCATTERPLOT(ie_lon, ie_lat, /OVERPLOT, SYMBOL='dot', SYM_SIZE=1, SYM_FILLED=1, SYM_FILL_COLOR='blue')
                    sp = SCATTERPLOT(b2_lon, b2_lat, /OVERPLOT, SYMBOL='dot', SYM_SIZE=1, SYM_FILLED=1, SYM_FILL_COLOR='light grey')
                    sp = scatterplot(ie_start_coord[1], ie_start_coord[0], /overplot, symbol = 'star', sym_size = 2, sym_filled = 1, sym_fill_color = 'green')
                    thiswin.save, plotsavefolder +lineid+ '-disttest.png'           
                    thiswin.close
                    undefine, thiswin, mPolStereo, mc, sp ; undefine them all (might help solve automation issues !)
                    
                    thiswin = window()
                    im2 = image(bytscl(sicmap+120.*mizedge_sic, min=0, max=120), image_dimensions=[1264,1328], /current)
                    label2 = symbol(mpinterpeasting, mpinterpnorthing , /DATA, 'Dot', SYM_SIZE = 1, SYM_COLOR = 'red')
                    tex = text(0,0,'beam2' +plab+ strtrim(filename,2))           
                    thiswin.save, plotsavefolder +lineid+ '-trackmap.png'
                    thiswin.close
                    undefine, thiswin, im2, label2, tex
                    
                    stop
                    
    ;                croppedim = sicmap[632:*,200:664]
    ;                croppedie = mizedge_sic[632:*,200:664]
    ;                im = image(bytscl(croppedim+120.*croppedie, min=0, max=120))
                    
                    win = window(dimension = [1000,400]) 
                    pl1 = plot(all_SIC_dist/1000., all_SIC_vals, /current, yrange = [min(all_SIC_vals),105],xtitle = 'distance from ice edge (km)', ytitle = 'SIC (%)', name = 'uncorrected distance')
                    pl2 = plot(all_SIC_corrdist/1000., all_SIC_vals, /overplot, color = 'deep sky blue', name = 'corrected distance')
                    pl3 = scatterplot(all_SIC_dist[good_mpinds]/1000., all_SIC_vals[good_mpinds], symbol = 'x', name = 'selected sections', /overplot)
                    pl4 = scatterplot(all_SIC_corrdist[good_mpinds]/1000., all_SIC_vals[good_mpinds], symbol = 'x', sym_color = 'deep sky blue', /overplot)
                    leg = legend(target = [pl1, pl2, pl3], position = [0.9,0.4])
                    
                    win.save, plotsavefolder+lineid+'-sicdist.png'
                    win.close
                    undefine, win, pl1, pl2, pl3, pl4, leg
                endif
        
                mean_error_sec = fltarr(sec_num,3, /nozero)
                std_error_sec = fltarr(sec_num,3, /nozero)
                mean_seglen_sec = fltarr(sec_num,3, /nozero)
                minmax_seglen_sec = fltarr(2,sec_num,3, /nozero)
                mptime = dblarr(sec_num,3, /nozero)
        
                ; calculate mean height error stats, seglen and time for each section
                for beam = 0, 2 do begin
                  if beam eq 0 then begin
                    distance = b1_dist & error = b1_error & seglen = b1_seglen & newx = b1_newx & time = b1_time
                  endif
                  if beam eq 1 then begin
                    distance = b2_dist & error = b2_error & seglen = b2_seglen & newx = b2_newx & time = b2_time
                  endif
                  if beam eq 2 then begin
                    distance = b3_dist & error = b3_error & seglen = b3_seglen & newx = b3_newx & time = b3_time
                  endif
                  mpind = dblarr(sec_num, /nozero)
                  for j = 0, sec_num-1 do begin
                    secstart =  newx[secind[0,beam,j]] 
                    secend = newx[secind[sec_len-1,beam,j]]
                    startdistdiff = abs(secstart-distance)
                    enddistdiff = abs(secend-distance)
                    thisstartind = where(startdistdiff eq min(startdistdiff)) & thisstartind = thisstartind[0] ; catch situations where there are multiple of same dist- usually height is nan here
                    thisendind = where(enddistdiff eq min(enddistdiff)) & thisendind = thisendind[0]
                    thiserror = error[thisstartind:thisendind]
                    thisseglen = seglen[thisstartind:thisendind]
                    mean_error_sec[j,beam] = mean(thiserror[where(thiserror lt 1000.)], /nan) ; some error values x10^38 (missing?) don't include these in calcs
                    std_error_sec[j,beam] = stddev(thiserror[where(thiserror lt 1000.)], /nan)
                    mean_seglen_sec[j,beam] = mean(thisseglen, /nan)
                    minmax_seglen_sec[*,j,beam] = [min(thisseglen), max(thisseglen)]
                    
                    secmp = newx[secind[sec_len/2.,beam,j]]
                    mpdistdiff = abs(secstart-distance)
                    thismpind = where(mpdistdiff eq min(mpdistdiff)) & thismpind = thismpind[0]
                    mpind[j] = thismpind          
                  endfor
                  mptime[*,beam] = time[mpind]
               endfor
    
               save, atljuliandate, lineid, plab, $ 
                  b1_height,b2_height,b3_height,b1_newx,b2_newx,b3_newx, b1_dist,b2_dist,b3_dist, b1_lat, b2_lat, b3_lat, b1_lon, b2_lon, b3_lon,$
                  b1_seglen,b2_seglen,b3_seglen,b1_error,b2_error,b3_error,$
                  sec_num, secind, secindall, mptime, mean_error_sec, std_error_sec, mean_seglen_sec, minmax_seglen_sec,$
                  ;mp_dist_from_ie,mp_dist_from_ie_SICcorr, mpSIC$          
                  all_SIC_vals, all_SIC_dist,all_SIC_coords,all_SIC_corrdist, good_mpinds, extrapol_inds, ie_start_coord, $
                  filename = outpath+lineid+'_interp.sav', /compress
                  
                ; to access mp lat and lon and SIC from big arrays (all_SIC_*) use good_mpinds on all_SIC_*
                ; extrapol inds are the index for projected coords (if ice edge to the north of the track) and info about these points
             
             endif ; end if findice undefined   
         endif ; end if sec_num > 0 
      endfor ; end part loop
    endif ; end sc orient if condition 2
endif ; end sc orient if condition 1

end

