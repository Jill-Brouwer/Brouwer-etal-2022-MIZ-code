function correct_dist, original_dist, SIC_fraction, step

  corrected_dist = dblarr(n_elements(original_dist), /nozero)
  startdist = 0
  corrected_dist[0] = startdist

  for i = 1, n_elements(corrected_dist)-1 do begin
    thissicfraction = sic_fraction[i]
    thisdist = step*thissicfraction
    totdist = startdist+thisdist
    corrected_dist[i] = totdist
    startdist = totdist
  endfor


return, corrected_dist

end