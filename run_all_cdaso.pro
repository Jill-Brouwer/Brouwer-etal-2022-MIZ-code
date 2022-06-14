pro run_all_Cdaso

; read in year and month from textfile and copy files into folder....
; get this file list 
; do one month at a time - copy month 1 from cdaso
; run all the processing
; delete raw files from cdaso (but keep the sav files)


; info for downloading from cdaso
  monthstring_all = ['02']
  yearstring_all = ['2019','2019','2019']
  username = 'brouwerj'
  servername = '144.6.252.33'
  hdffolder = '/cdaso/PRIVATE/brouwerj/icesat2'
  workspacepath = '/home/brouwerj/batch/cpu1'
  prelimfiltpath = '/jill_hdd/Batch_process/'
  continue_from = 0

; specify some important variables before running
  sample_int = 8. ;distance between interpolated x samples (m)
  sec_dist = 6.25 ; section distance (km)
  sec_len = long(sec_dist*1000./sample_int) ; number of points in the array (add 1 because 0 isn't a distance)
  nan_thresh = 50. ; how many missing data points accepted in spec calcs (%)
  prelim_filter_thresh = 50. ; what level to prelim filter on 
  processing_date = '2019-10-23'
  good_line_count = 0 ; counter for checking all good files have been processed 
  
;  ; load in SIC data
  restore, '/home/brouwerj/sav_files/SIC_processed_small.sav' ; most recent combined .sav with all SIC
  sicfile = '/home/brouwerj/sav_files/sic_E&N.sav'
;  print, 'sicfiles loaded'
;  
;  ; loop over months 
;  
  for month = 0, n_elements(monthstring_all)-1 do begin
    
    monthstring = monthstring_all[month]
    yearstring = yearstring_all[month]
    
    ; where to save output .sav files and plots
      interppath = prelimfiltpath + yearstring + '_' + monthstring + '/interpol/' ; variables 1 from interpol and section selection
      hspath =  prelimfiltpath + yearstring + '_' + monthstring + '/hs/'
      plotpath = prelimfiltpath + yearstring + '_' + monthstring + '/plots/'
;        
;    ; names of savefiles 
      prelimfiltsavname = prelimfiltpath + 'prelimfilt_nanthresh_'+strtrim(prelim_filter_thresh, 2)+'_'+yearstring+monthstring+'.sav'
    
    ; get filenames from cdaso 
    spawn, "ssh "+username+"@"+servername+" 'find "+hdffolder+" -type f -name *"+yearstring+monthstring+"*.h5 -maxdepth 1 | sort '", hdflist
      print, 'found all good files'
      
      stop
      
      hdflist = hdflist[where(strmatch(hdflist,'*20190218101142_07940201*'))]
      
      stop
      
    ; set up a log
;    spawn, 'echo processing for:'+monthstring+'_'+yearstring+' completed on '+processing_date+' >> '+prelimfiltpath+YEARSTRING+'_'+MONTHSTRING+'/log.txt'
;    spawn, 'echo prelim filtthresh = '+string(prelim_filter_thresh)+'hs thresh ='+string(nan_thresh)+' >>'+prelimfiltpath+YEARSTRING+'_'+MONTHSTRING+'/log.txt'
;    spawn, 'echo total filenum ='+string(n_elements(Hdflist))+' >> '+prelimfiltpath+YEARSTRING+'_'+MONTHSTRING+'/log.txt'
;    
    ; copy file one by one to workspace and run interpol on them 
    for file = continue_from, n_elements(hdflist)-1 do begin
      print, 'reading file' + string(file) + 'out of' + string(n_elements(hdflist)-1)
    
      ; get atl filename (without the cdaso path stuff) for saving
      thishdf = hdflist[file]    
      
      namesplit = strsplit(thishdf, '/')
      atlfilename = thishdf.substring(namesplit[n_elements(namesplit)-1]) 
      
      ; copy file
      spawn, "scp "+username+"@"+servername+":"+thishdf+" "+workspacepath
      localfile = workspacepath+'/'+atlfilename
       
      ; get matching SIC data (do it here so don't have to restore inside multiple functions)
      atldate = atlfilename.substring(9,16)
      sicday1 = sicdate2day[*,0]
      sicday2 = sicdate2day[*,1]
      thissicday = where((strmatch(sicday1, atldate) eq 1) or (strmatch(sicday2, atldate) eq 1))
      SICmap = float(reform(sicmaps2day_new[*,*,thissicday]))
      SICmap[where(SICmap ge 250)] = !values.F_nan
      mizedge_SIC = mizedge_sic_2day[*,*,thissicday]
        
      ; get prelim filter check
      ; output = [good_start_flag, good_end_flag]
      goodflag = prelim_filter_cdaso(prelim_filter_thresh,localfile,sicmap,mizedge_sic,sicfile)  
    
      spawn, 'echo '+atlfilename+','+string(goodflag[0])+','+string(goodflag[1])+' >> '+prelimfiltpath+YEARSTRING+'_'+MONTHSTRING+'/log.txt'
    
      ;run interpol
      ; inputs are: filename, good_start_flag, good_end_flag, sicmap, mizedge_sic, sicpath, outpath,plotpath, sec_len, sample_int, nan_thresh
      ;interpol_and_segment_cdaso, localfile, goodflag[0], goodflag[1], sicmap, mizedge_sic, sicfile, interppath, plotpath, sec_len, sample_int, nan_thresh ;, do_plots = 1
      interpol_and_segment_cdaso, localfile, 0, 1, sicmap, mizedge_sic, sicfile, interppath, plotpath, sec_len, sample_int, nan_thresh ;, do_plots = 1

;      ; delete the file on the local drive
;      spawn, "rm "+localfile
;      
    endfor
;      
    ; restore interpolated section data from previous script
    spawn, 'find ' + interppath + ' -name "*.sav" -maxdepth 1 |sort', interpfilenames
    
    spectral_analysis_batch, interpfilenames, sample_int, sec_len, hspath, nan_thresh
  endfor

stop
end
