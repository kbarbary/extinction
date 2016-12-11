
pro run

  minwave = 1000
  maxwave = 50000
  n = 200
  
  ; create log-spaced wavelength vector
  minexp = alog10(minwave)
  maxexp = alog10(maxwave)
  step = (maxexp - minexp) / (n-1)
  exps = minexp + step * findgen(n)
  wave = 10.^exps

  ; create a "flat" flux vector
  flux = wave * 0.0 + 1.0

  r_vs = [2.3, 3.1, 4.0, 5.3]

  fnames = ['fm_unred_2.3.dat', 'fm_unred_3.1.dat', 'fm_unred_4.0.dat', 'fm_unred_5.3.dat']

  for i = 0, n_elements(r_vs) - 1 do begin
     a_v = 1.0
     r_v = r_vs(i)
     ebv = a_v / r_v
     fm_unred, wave, flux, -ebv, flux_new, R_V=r_v ;Redden (negative E(B-V))

     ; fm_unred does `result = flux * 10.^(0.4*curve)` where is scaled by input
     ; undo it to get the curve

     a_lambda = -2.5 * alog10(flux_new)

     data = fltarr(2, n_elements(wave))
     data(0,*) = wave
     data(1,*) = a_lambda

     openw, 1, fnames(i)
     printf, 1, data
     close, 1
  endfor
end
