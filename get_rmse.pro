pro get_rmse

  monthstring = ['02','05','09','12']
  
  ; get total number of goodfiles for the four months
  goodcount = lonarr(4)
  for i = 0, n_elements(monthstring)-1 do begin
    file= '/jill_hdd/Batch_process/2019_'+monthstring[i]+'/mizfind/plots/2019_'+monthstring[i]+'-goodlog.txt'
    routpath = '/jill_hdd/Batch_process/2019_'+monthstring+'/mizfind/Rout/'
    data=read_csv(file)   
    goodbad=data.field2
    numgood = n_elements(where(strmatch(goodbad,'good')eq 1))
    goodcount[i] = numgood
  endfor
  
  Feb_inds = [0:goodcount[0]-1]
  May_inds = [goodcount[0]:total(goodcount[0:1])-1]
  Sep_inds = [total(goodcount[0:1]):total(goodcount[0:2])-1]
  Dec_inds = [total(goodcount[0:2]):total(goodcount[0:3])-1]
  
  tot = total(goodcount)

  RMSE = fltarr(14,2,2,tot) ; hsm,distm,lin/exp
  MIZ_err = fltarr(14,2,2,tot)
  mean_ice_noise = fltarr(14,tot) ; mean ice noise from miz end to 2*miz end
  std_ice_noise = fltarr(14,tot)
  MIZ_dist_est = fltarr(14,2,2,tot)
  lat0lon0 = dblarr(2,tot)
  SICdist = fltarr(tot)
  
  for i = 0, n_elements(monthstring)-1 do begin
    print, i
    if i eq 0 then savinds = Feb_inds
    if i eq 1 then savinds = May_inds
    if i eq 2 then savinds = Sep_inds
    if i eq 3 then savinds = Dec_inds
    
    ; get the relevant files
    rootpath = '/jill_hdd/Batch_process/2019_'+monthstring[i]+'/mizfind/'
    goodfile= rootpath+'plots/2019_'+monthstring[i]+'-goodlog.txt'
    routpath = rootpath+'Rout/'
    csvpath = rootpath+'csv/'

    ; read in good line IDs per month (again) - soz
    data=read_csv(goodfile)
    goodbad=data.field2
    filepath=data.field1
    trimmed_filepath=file_basename(filepath, '-mizcheck.pdf')
    wheregood = where(strmatch(goodbad,'good')eq 1)
    thisnumgood = n_elements(wheregood)
    goodfileids = trimmed_filepath[wheregood]

    ; mini arrays for saving a single month
    this_RMSE = fltarr(14,2,2,thisnumgood) ; hsm,distm,lin/exp
    this_MIZ_err = fltarr(14,2,2,thisnumgood)
    this_mean_ice_noise = fltarr(14,thisnumgood) ; mean ice noise from miz end to 2*miz end
    this_std_ice_noise = fltarr(14,thisnumgood)
    this_lat0lon0 = dblarr(2,thisnumgood)
    this_MIZ_dist_est = fltarr(14,2,2,thisnumgood)
    this_SICdist = fltarr(thisnumgood)

    
    for fnum = 0, thisnumgood-1 do begin
      
      lineID = goodfileids[fnum]
      spawn, 'find '+csvpath+' -name "*'+lineID+'*" -maxdepth 1', csvfile
      spawn, 'find '+routpath+' -name "*'+lineID+'*" -maxdepth 1', routfile

      rout= read_csv(routfile,table_header = routhead,types = ['string','float','float','float','float'])
      csvin = read_csv(csvfile,table_header = csvinhead)

      ; need to replace NAs with a number for it to be numeric - do 999
      strnames = tag_names(rout)

      for s = 0, n_elements(strnames)-1 do begin
        rout.(s)[where(strmatch(rout.(s),'NA')eq 1)]=9999999
      endfor
      
      ; headers for checking i've grabbed the right ones
      routheads = ["Hsm", "DC_slope_lin","DC_slope_lin_err","DC_int_lin","DC_int_lin_err","DC_bp_lin","DC_miz_lin","DC_RMSE_lin","DC_miz_lin_err","DC_slope_exp","DC_slope_exp_err","DC_exp_lin","DC_int_exp_err","DC_bp_exp","DC_miz_exp","DC_RMSE_exp","DC_miz_exp_err","D_slope_lin","D_slope_lin_err","D_int_lin","D_int_lin_err","D_bp_lin","D_miz_lin","D_RMSE_lin","D_miz_lin_err","D_slope_exp","D_slope_exp_err","D_exp_lin","D_int_exp_err","D_bp_exp","D_miz_exp","D_RMSE_exp","D_miz_exp_err","DC_segslop_lin","DC_segslop_lin_err","DC_segint_lin","DC_bp_lin_err","DC_segslop_exp","DC_segslop_exp_err","DC_segint_exp","DC_bp_exp_err","D_segslop_lin","D_segslop_lin_err","D_segint_lin","D_bp_lin_err","D_segslop_exp","D_segslop_exp_err","D_segint_exp","D_bp_exp_err","SIC_dist"]
      csvinheads = ["Dist","CDist","hm0bc","hm0hann","std","F1500","F1038","F719","F498","F345","F239","F165","F114","F79","F55","F38","hm0bcerr","hm0hannerr","stderr","F1500err","F1038err","F719err","F498err","F345err","F239err","F165err","F114err","F79err","F55err","F38err","mplat","mplon","mpsic","good_mpinds","all_SIC_dist","all_SIC_corrdist","all_SIC_vals","all_lat","all_lon"]

      this_RMSE[*,0,0,fnum] = float(rout.(7)) ; dc_RMSE_lin
      this_RMSE[*,0,1,fnum] = float(rout.(15)) ; dc_RMSE_exp
      this_RMSE[*,1,0,fnum] = float(rout.(23)) ; d_RMSE_lin
      this_RMSE[*,1,1,fnum] = float(rout.(31)) ; d_RMSE_exp

      this_miz_err[*,0,0,fnum] = float(rout.(8)) ;dc_miz_lin_err
      this_miz_err[*,0,1,fnum] = float(rout.(16)) ; dc_miz_exp_err
      this_miz_err[*,1,0,fnum] = float(rout.(24)) ; d_miz_lin_err
      this_miz_err[*,1,1,fnum] = float(rout.(32)) ; d_miz_exp_err
      
      this_sicdist[fnum] = float(rout.(49)[0]) ; miz dist
      this_lat0lon0[0,fnum] = double(csvin.(37)[0]) ;lat of track at ice edge
      this_lat0lon0[1,fnum] = double(csvin.(38)[0]) ; lon of track at ice edge
      
      this_miz_dist_est[*,0,0,fnum] = float(rout.(6)) ; dc_miz_lin
      this_miz_dist_est[*,0,1,fnum] = float(rout.(14)) ; dc_miz_exp
      this_miz_dist_est[*,1,0,fnum] = float(rout.(22)) ; d_miz_lin
      this_miz_dist_est[*,1,1,fnum] = float(rout.(30)) ; d_miz_exp
      
      ; figure out the mean ice noise after the miz break
      mpdist = csvin.(1)
      allmizdistest = this_miz_dist_est[6:9,*,*,fnum]
      allmizdistest[where(allmizdistest gt 9000000)] = !values.F_nan
      max_miz_est = max(allmizdistest, /nan)
      maxinds = where((mpdist ge max_miz_est) and (mpdist le 2*max_miz_est))
      
      if n_elements(maxinds) gt 2 then begin
        for hsm = 0, 13 do begin

          csvindex = hsm+2 ; find the right
          thisheight = csvin.(csvindex)
          mps = thisheight[maxinds]

          this_mean_ice_noise[hsm,fnum] = mean(mps, /nan)
          this_std_ice_noise[hsm,fnum] = stddev(mps, /nan)
        endfor
      endif else begin
        this_mean_ice_noise[*,fnum] = !values.f_nan
        this_std_ice_noise[*,fnum] = !values.f_nan
      endelse

    endfor
    
    RMSE[*,*,*,savinds] = this_RMSE
    MIZ_err[*,*,*,savinds] = this_miz_err
    mean_ice_noise[*,savinds] = this_mean_ice_noise
    std_ice_noise[*,savinds] = this_std_ice_noise
    MIZ_dist_est[*,*,*,savinds] = this_miz_dist_est
    lat0lon0[*,savinds] = this_lat0lon0
    SICdist[savinds] = this_SICdist
    
  endfor
  
stop
   
   ; remember to replace the 9999s
   RMSE[where(RMSE gt 9000000)] = !values.F_nan
   miz_err[where(miz_err gt 9000000)] = !values.f_nan
   mean_ice_noise[where(mean_ice_noise gt 9000000)] = !values.f_nan
   std_ice_noise[where(std_ice_noise gt 9000000)] = !values.f_nan
   miz_dist_est[where(miz_dist_est gt 9000000)] = !values.f_nan
   sicdist[where(sicdist gt 9000000)] = !values.f_nan

   Firfmiz_alex = median(miz_dist_est[6:9,*,*,*],/even,dimension = 1)
   
   
  save, RMSE, miz_err, mean_ice_noise, std_ice_noise, miz_dist_est, sicdist,lat0lon0,firfmiz_alex,Feb_inds,May_inds,Sep_inds,Dec_Inds, /compress, filename = '/home/brouwerj/sav_files/results.sav'
stop
  ; read all the info into the big arrays



    
    
    
   
    
    
    


  
 
  stop
  
  ; get average RMSE value and number of nans for each hsm
  hsmrmseave = fltarr(14)
  hsmrmsestd = fltarr(14)
  hsmnnans = fltarr(14)
  for hsm = 0, 13 do begin
      thishsrmse = RMSE[hsm,*,*,*]
      hsmrmseave[hsm] = mean(thishsrmse, /nan)
      hsmrmsestd[hsm] = stddev(thishsrmse, /nan)
      hsmnnans[hsm] = n_elements(where(~finite(thishsrmse)))      
  endfor
 
  rmse_dist1 = fltarr(2*tot*14)
  rmse_dist2 = rmse_dist1
  rmse_dist1 = reform(RMSE[*,*,0,*])
  rmse_dist2 = reform(RMSE[*,*,1,*],[2*tot*14])
  
  forplot1 = rmse_dist1[where(rmse_dist1 le 10000)]
  forplot2 = [rmse_dist2[where(rmse_dist1 le 10000)]]
  
  forplot1 = forplot1[sort(forplot1)]
  forplot2 = forplot2[sort(forplot2)]
  
  ;boxplot([[forplot1],[forplot2]])
  
stop
  
  RMSE = fltarr(14,2,2,numgood)

stop
end