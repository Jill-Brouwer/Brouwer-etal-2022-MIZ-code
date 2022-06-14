pro get_into_R


  ; read in Hs data and matching interpol
  monthstring = ['02','05','09','12'] ;'01',
  monthstring = ['09']
  monthstring = ['05']
  yearstring = '2019'
  processing_date = '18_10_2020'

;monthstring = '09'

for month = 0,n_elements(monthstring)-1 do begin
  
  thismonthstring = monthstring[month]
  print, monthstring[month]
  
  file_path = '/jill_hdd/Batch_process/'+yearstring+'_'+thismonthstring+'/'
  
  ;file_path = '/jill_hdd/Batch_process/2019_09/old/'
  hsfilepath = file_path+'hs/'
  interpfilepath = file_path+'interpol/'
  mizfindpath = file_path+'mizfind/'

  simonpath = mizfindpath+'csv/'
 

  spawn, 'find ' + hsfilepath + ' -name "*.sav" -maxdepth 1 |sort', hsfilenames
  n_lines = n_elements(hsfilenames)
  
  hsfilenames = hsfilenames[where(strmatch(hsfilenames,'*20190523232549_08510301*'))]

  for f=0, n_elements(hsfilenames)-1 do begin

    ; load in data for each file
    restore, hsfilenames[f]
    spawn, 'find ' + interpfilepath + ' -name "*' + lineID + '*" -maxdepth 1', thisinterpfile ; find matching interp data based on file id
    if n_elements(thisinterpfile) gt 1 then begin
      print, 'more than one matching interp file'
      stop
    endif
    restore, thisinterpfile

    ; extra criteria to filter out when first sic long is between 50 and 61W
    if (all_SIC_coords[1,0] lt -50) and (all_SIC_coords[1,0] gt -61) then continue


    ; calc variables - mean hs, std, tot_error, correct for snip problem (new sav files since 30 sep will be fixed)
    ; change distances to km

    all_SIC_dist_km = all_SIC_dist/1000.
    all_SIC_corrdist_km = all_SIC_corrdist/1000.
    mp_dist_from_ie = float(all_SIC_dist_km[good_mpinds])
    mp_dist_from_ie_SICcorr = float(all_SIC_corrdist_km[good_mpinds])
    mplat = reform(all_SIC_coords[0,good_mpinds])
    mplon = reform(all_SIC_coords[1,good_mpinds])
    mpSIC = all_SIC_vals[good_mpinds]

    ; combine spec and firf ests into one big array (easier from automation side) spec comes before FIRF
    allhs = fltarr(n_elements(good_mpinds), 14, 3)
    allhs[*,0:2,*] = spechs
    allhs[*,3:*,*] = firfhs

    mhs = mean(allhs, dimension = 3)
    stdhs = stddev(allhs, dimension = 3)
    height_error = mean(mean_error_sec, dimension = 2)
    
    if sec_num lt 10 then continue

    ; calculate total quadrature error
    toterr = fltarr(n_elements(good_mpinds), 14, /nozero)
    for hsm = 0, 13 do begin
      thishs = stdhs[*,hsm]
      toterr[*,hsm] = sqrt(thishs^2+height_error^2)
    endfor
    
    ; get raw heights snipped to the section inds (only the data that was used) 
    thisb1_height =b1_height[secind[0,0,0]:secind[780,0,sec_num-1]]
    thisb2_height = b2_height[secind[0,1,0]:secind[780,1,sec_num-1]]
    thisb3_height = b3_height[secind[0,2,0]:secind[780,2,sec_num-1]]

    disttoiceedgeb2 = ll2dist(b2_lat[0],b2_lon[0], ie_start_coord[0],ie_start_coord[1])
    extradist = b2_newx[secind[0,1,0]]

    thisb1_dist = b1_newx[secind[0,0,0]:secind[780,0,sec_num-1]]
    thisb1_dist = thisb1_dist + disttoiceedgeb2+extradist
    thisb2_dist = b2_newx[secind[0,1,0]:secind[780,1,sec_num-1]]
    thisb2_dist = thisb2_dist  + disttoiceedgeb2+extradist
    thisb3_dist = b3_newx[secind[0,2,0]:secind[780,2,sec_num-1]]
    thisb3_dist = thisb3_dist  + disttoiceedgeb2+extradist
    ;
    ;    pl = plot(thisb1_dist/1000., thisb1_height, layout=[1,3,1], yrange = [-4,4])
    ;    pl = plot(thisb2_dist/1000., thisb2_height, layout=[1,3,2], /current, yrange = [-4,4])
    ;    pl = plot(thisb3_dist/1000., thisb3_height, layout=[1,3,3], /current, yrange = [-4,4])

    length = n_elements(thisb2_dist)



    csvname = simonpath+lineID+'-mizfind.csv'
    
    bigarr = fltarr(n_elements(all_SIC_dist),39)
    bigarr[*] = !values.f_nan ; initialise with nan
    
    mp_num = sec_num-1
    bigarr[0:mp_num,0] = mp_dist_from_ie
    bigarr[0:mp_num,1] = mp_dist_from_Ie_SICcorr
    bigarr[0:mp_num,2:15] = mhs
    bigarr[0:mp_num,16:29] = toterr
    bigarr[0:mp_num,30] = mplat
    bigarr[0:mp_num,31] = mplon
    bigarr[0:mp_num,32] = mpsic
    bigarr[0:mp_num,33] = good_mpinds
    bigarr[*,34] = all_SIC_dist_km
    bigarr[*,35] = all_SIC_corrdist_km
    bigarr[*,36] = all_SIC_vals
    bigarr[*,37] = all_SIC_coords[0,*]
    bigarr[*,38] = all_SIC_coords[1,*]
    
    Firfnames = 'F'+strtrim(round(centre_wl),2)
    firferr = firfnames+'err'
    
    bigarr = transpose(bigarr)
    
    write_csv, csvname, bigarr, $
      header = ['Dist','CDist','hm0bc','hm0hann','std',firfnames,'hm0bcerr','hm0hannerr','stderr',firferr,'mplat','mplon','mpsic', 'good_mpinds', 'all_SIC_dist','all_SIC_corrdist', 'all_SIC_vals', 'all_lat', 'all_lon']

  endfor

endfor


stop
end