pro sic_processing_2day_combined

  ; get lat and lon arrays ready

  ;paths for lat and lon of grid
  areapath = '/jill_hdd/sept2019_SIC/grid_coordinates/pss06area_v3.dat'
  latspath = '/jill_hdd/sept2019_SIC/grid_coordinates/pss06lats_v3.dat'
  lonspath = '/jill_hdd/sept2019_SIC/grid_coordinates/pss06lons_v3.dat

  ; create empty arrays for variables
  area = lonarr(1264, 1328) 
  lats = lonarr(1264, 1328)
  lons = lonarr(1264, 1328)

  OPENR, lun, areapath, /GET_LUN
  readu, lun, area
  free_lun, lun

  OPENR, lun, latspath, /GET_LUN
  readu, lun, lats
  free_lun, lun

  OPENR, lun, lonspath, /GET_LUN
  readu, lun, lons
  free_lun, lun

  area = area/1000.0
  lats = lats/100000.0
  lons = lons/100000.0

  ; check alignment with lons but everything is the same
  ;cgimage, bytscl(lons)

  ; flip vertically

  area = rotate(area, 7)
  lats = rotate(lats, 7)
  lons = rotate(lons, 7)

  ; create easting and northing arrays
  ; make space

  SICeasting = fltarr(1264)
  SICnorthing = fltarr(1328)

  for i = 0, 1263 do begin
    temp = ps2xy_chadstyle(lats[i,0], lons[i,0], -70.0, 0.0)
    SICeasting[i] = temp[0]
  endfor

  for i = 0, 1327 do begin
    temp = ps2xy_chadstyle(lats[0,i], lons[0,i], -70.0, 0.0)
    SICnorthing[i] = temp[1]
  endfor

  SIClats = lats & SIClons = lons & SICarea = area
  save, SICeasting, SICnorthing, SIClats, SIClons, SICarea, filename = '/home/brouwerj/sav_files/sic_E&N.sav'
  
  ; read in data  ; NEED TO CHANGE AND SORT FILES IN HERE WITH 
  ;SICpath = '/jill_hdd/sept2019_SIC/hd5/'
  SICpath = '/jill_hdd/all_SIC'
  SIClist = FILE_SEARCH(SICpath,'*hdf*')
  ndays = n_elements(SIClist) ; specify number of days we have SIC to hard code arrays
  SICdate = strarr(ndays)
  SICmaps = fltarr(1264, 1328, ndays)
  
  ; save date for matching with IS-2 (later)
  for n=0, ndays-1 do begin
    SICfilesplit = strsplit(SIClist[n], '-', /EXTRACT)
    SICdate[n] = SICfilesplit[3]
    SICmaps[*,*,n] = HDFreturnSDS(SIClist[n], "ASI Ice Concentration")
    print, SICdate[n]
  endfor
  
  ; combine 2 day together by selecting max of each 
  SICdate2day = strarr(ndays)
  SICmaps2day = fltarr(1264,1328,ndays)
  npixels = n_elements(SICmaps[*,*,0])
  SICmax = fltarr(1264,1328)

  for n=0, ndays/2-1 do begin
    d1ind = 2.*n & d2ind = d1ind+1.
    SICdate2day[d1ind] = SICdate[d1ind]
    SICdate2day[d2ind] = SICdate[d2ind]
    SICmapday1 = reform(SICmaps[*,*,d1ind])
    SICmapday2 = reform(SICmaps[*,*,d2ind])
    for p=0, npixels-1 do begin
      SICmax[p] = max([SICmapday1[p],SICmapday2[p]])
    endfor
    SICmaps2day[*,*,d1ind] = SICmax
    SICmaps2day[*,*,d2ind] = SICmax
    print, d1ind, d2ind
  endfor
  
  ; calculate ice edge for each map
  mizedge_sic_2day=bytarr(1264,1328,ndays)
  increment=0.1

  for n=0, ndays/2-1 do begin
    d1ind = 2.*n & d2ind = d1ind+1.
    sicmap = reform(sicmaps2day[*,*,d1ind])
    for i=-180.0, 180.0, increment do begin
      print, i
      lonslice = lons gt i and lons lt (i+increment)
      ;now basically find the northernmost one of these with ge 15% conc....
      ge15perc_thisslice=lonslice*(sicmap ge 15.)
      ;get the lat of all these and find max..
      lats_thisslice=lats[where(ge15perc_thisslice)]
      maxlat=max(lats_thisslice)
      ;now the index of that one....
      thisind=where(lats eq maxlat and lons gt i and lons lt (i+increment))
      ;decompose into x and y... How does this work?? thisind is a 1d subscript, need to figure out where it is in a 1264x1328 array
      thisx=thisind mod 1264 ; divides by 1264 and calculates the remainder (very cool)
      thisy=thisind / 1264
      ;and fill in the map at this px....
      
      mizedge_sic_2day[thisx, thisy,d1ind]=1
      mizedge_sic_2day[thisx, thisy,d2ind]=1
      
    endfor
 ;  cgwindow, 'cgimage', bytscl(sicmaps2day[*,*,29]+120.*mizedge_sic_2day[*,*,29],min=0,max=120)
;    cgwindow, 'cgimage', bytscl(sicmaps2day[*,*,d2ind]+120.*mizedge_sic_2day[*,*,d2ind],min=0,max=120)
  endfor
  
  stop

save, sicmaps2day, sicdate2day, mizedge_sic_2day, SICeasting, SICnorthing, SIClats, SIClons, SICarea, filename = '/jill_hdd/sic_2day_processed_ALL.sav'

stop

end
