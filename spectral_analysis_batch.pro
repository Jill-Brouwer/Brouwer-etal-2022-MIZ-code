pro spectral_analysis_batch, filenames, sample_int, sec_len, outpath, nan_thresh


; create fft windows
hwin = hanning(sec_len)
bcwin = findgen(sec_len, start=1, increment=0)


; firf stuff
; create a string of filter filenames in order of descending wavelength

dirs = '/home/brouwerj/plots/filterbank'
filter = '*.sav'
fpaths = dirs + path_sep() + filter  ; Append the filter to each path
filtlist = file_search( fpaths )
order = []
for f =0, n_elements(filtlist)-1 do begin
  filtname = filtlist[f]
  fsplit =  strsplit(filtname, '_')
  fsplit2 = strsplit(filtname, 'm')
  split = filtname.Substring(fsplit[2], fsplit2[2]-2)
  order = [order, float(split)]
endfor

filtlists = filtlist[reverse(sort(order))] ; filenames in order from longest wl to shortest
filt_num = n_elements(filtlists)
thresh_perc = 50. ; choose what threshold percentage (here for now but can also be specified when making the filters)

; run analysis for all files
for file = 0, n_elements(filenames)-1 do begin
  
  restore, filenames[file]
  print, 'starting line' + lineid
  
  sec_num = n_elements(secind[0,0,*])

  ; firf processing starts here

  ; apply the convolution (to all beams)

  print,'start firf'

  restore, filtlists[0]
  kernel1 = fir_filt.fir_filt
  test=convol_gaps(b1_height,kernel1,dud_above=dud_above,n_thresh=n_thresh)
  centre_wn = fltarr(filt_num)

  print,'got firf data ready'

  ;filt_res = fltarr(n_elements(test), filt_num,3) ; [convolution result, filter #, beam #]

  b1firf = fltarr(n_elements(b1_height), filt_num)
  b2firf = fltarr(n_elements(b2_height), filt_num)
  b3firf = fltarr(n_elements(b3_height), filt_num)
  firfheightsec = fltarr(sec_len, sec_num, filt_num, 3)

  for i = 0, filt_num-1 do begin ;
    restore, filtlists[i]
    kernel = fir_filt.fir_filt
    centre_wn[i] = fir_filt.centre_freq
    n_thresh = thresh_perc/100.*fir_filt.n_span
    for beam = 0, 2 do begin
      if beam eq 0 then height_regular_snip = b1_height
      if beam eq 1 then height_regular_snip = b2_height
      if beam eq 2 then height_regular_snip = b3_height
      firf = convol_gaps(height_regular_snip, kernel, dud_above=dud_above, n_thresh=n_thresh)
      firf[where(firf gt 900)]=!Values.F_nan
      print, 'filter ' + strtrim(i,2) + 'beam ' + strtrim(beam+1) + ' convolved'
      ; calculate stdev inds of sections are in startind[beam#,section#] to get endind add sec_len-1 (780)
      for j = 0, sec_num-1 do begin
        starti = secind[0,beam,j]
        endi = starti + sec_len-1
        ; check that start ones are not above 50% threshold and if so set value of whole section to be 999 - not done because if est a std of only a few points the variability will be larger anyway
        firfsec = firf[starti:endi]
        firfheightsec[*,j,i,beam] = firf[starti:endi]
        ;nnans = float(n_elements(where(~finite(firfsec))))
        ;if nnans/sec_len*100. gt nan_thresh then firfheightsec[*,j,i,beam] = !values.F_nan else firfheightsec[*,j,i,beam] = firf[starti:endi]
      endfor
      ; save raw height data out of loop (just in case + for plotting)
      if beam eq 0 then b1firf[*,i] = firf
      if beam eq 1 then b2firf[*,i] = firf
      if beam eq 2 then b3firf[*,i] = firf
    endfor
  endfor

  centre_wl = (1./centre_wn)*8.

  ; fft processing starts here:
  print,'start fft'

  ; calculate missing data profiles
  heightsec = fltarr(sec_len, sec_num, 3)
  heightsec[*,*,0] = b1_height[reform(secind[*,0,*])]
  heightsec[*,*,1] = b2_height[reform(secind[*,1,*])]
  heightsec[*,*,2] = b3_height[reform(secind[*,2,*])]
  
  
  mdp = MAKE_ARRAY(sec_len,sec_num,3, /INTEGER, VALUE = 1)
  i_nan_test = where(~finite(heightsec), /null)
  mdp[i_nan_test] = 0

  secmeans = mean(heightsec, dimension = 1, /nan)

  ; use a test FFT to get the number of elements for each FFT (instead of hard coding the number)
  testheight = heightsec[*,0,0]
  nanind = where(~finite(testheight), n_nans, /null)
  testheight[nanind] = 0
  testFFT = FFT_powerspectrum(reform(testheight), sample_int, FREQ = wn)
  FFT_num = n_elements(testFFT)

  print,'get data fft'

  ; process data - zero mean, get rid of NANS, and apply the window function
  ; USE OUR OWN CALC OF STDEV TO GET SWH
  fhdata = fltarr(fft_num, sec_num, 3)
  fbdata = fltarr(fft_num, sec_num, 3)
  num_nans = fltarr(sec_num,3)

  for b = 0, 2 do begin
    for j = 0, sec_num-1 do begin
      this_mdp = mdp[*,j,b]
      mdp_hann = this_mdp*hwin ; combine missing data profile and hanning window for w_ss correction
      thisheightsec = heightsec[*,j,b] - secmeans[j,b]
      heightseccp = thisheightsec
      nanind = where(~finite(thisheightsec), n_nans, /null)
      heightseccp[nanind] = 0
      heightsecnonan = heightseccp
      hannheights = hwin*heightsecnonan
      num_nans[j,b] = float(n_nans)

      ; run the fft
      fallbc = FFT_powerspectrum(heightsecnonan, sample_int)
      fallhann = FFT_powerspectrum(hannheights, sample_int)

      ; calculate and apply wss corrections
      bc_wss = double(sec_len*total(this_mdp^2))
      hann_wss = double(sec_len*total(mdp_hann^2))

      bc_corrected = fallbc*(sec_len^2/bc_wss)
      hann_corrected = fallhann*(sec_len^2/hann_wss)

      fbdata[*,j,b] = bc_corrected
      fhdata[*,j,b] = hann_corrected

    endfor
  endfor

  percnanssec = num_nans/sec_len*100.
  print,'fft completed'

  ;calculate swh
  print,'start swh calcs'

  firfhs = 4*stddev(firfheightsec, dimension = 1, /nan)
  swh_std = 4*(stddev(heightsec, dimension = 1, /nan))
  
; bandpass filter the spectral data
  highwn = 1./16.
  lowwn = 1./1500.

  bc_m0_bpf = fltarr(sec_num,3)
  hann_m0_bpf = fltarr(sec_num,3)

  for b = 0, 2 do begin
    for j = 0, sec_num-1 do begin
      thisfbc = fbdata[*,j,b]
      thisfhann = fhdata[*,j,b]
      bpf_ind = where(wn lt highwn and wn gt lowwn)
      bc_m0_bpf[j,b] = total(thisfbc[bpf_ind])
      hann_m0_bpf[j,b] = total(thisfhann[bpf_ind])
    endfor
  endfor

  print,'swh spectrum loop fin'

  ; calc h_s from zeroeth moment of corrected spectrum
  bc_hm0_bpf = 4*sqrt(bc_m0_bpf)
  hann_hm0_bpf = 4*sqrt(hann_m0_bpf)
  
  ; combine spectral estimates into one array (sec_num, 0:bc_hm0/1:hann_hm0/2:hs_std, beam_num)
  spechs = fltarr(sec_num, 3, 3)
  spechs[*,0,*] = bc_hm0_bpf
  spechs[*,1,*] = hann_hm0_bpf
  spechs[*,2,*] = swh_std
  

  ; save variables
  save, b1firf, b2firf, b3firf, firfhs, centre_wl, spechs, fbdata, fhdata, wn, percnanssec, $
    lineID, filename = outpath+lineID+'_spec.sav', /compress
  
endfor


end
