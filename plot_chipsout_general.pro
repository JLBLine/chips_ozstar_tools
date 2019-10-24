;; IDL code to take the output binary files from CHIPS and plot them with Bryna's plotting code.
;; The code computes and applies the normalisations to take the units from Jy^2Hz^2sr^2 to mK^2 Mpc^3
;; and transposes and folds the data across the kpar = 0 line. The code will plot the output from xx and yy pols.

;; Inputs:
;; output_tag: Enter the output_tag that was used to run CHIPS (Entry required)
;; initials: Enter in your initials that were used to create your personal CHIPS_OUT directory (Default: empty)
;; FHD, RTS: Specify FHD or RTS inputs with those keywords, set to 1 in case of True.
;; oneD, twoD: Specify 1D and 2D plotting mode, set to 1 in case of True.
;; n_freq: Enter the number of frequency channels, defaults to 384
;; band, base_freq: Enter the band name (high, low, ultra-low, defaults to high) or the lowest frequency in the band in Hz.
;;   Setting band will default the base_freq to the expected value. Do not use if analyzing subbands
;; chan_width: Enter the frequency resolution in Hz, defaults to 80e3
;; lssa_num: Specify whether kriging has been performed on the input (0: kriging applied, 1: kriging not applied, defaults to 0)
;; band_point_weight: An array of 9 integers to indicate what pointings, from -4 to 4, were used. Default is just zenith.
;    For example, if pointings -2 to 2 were used, then band_point_weight = [0,0,1,1,1,1,1,0,0]
;; output_dir: Output directory for plots

;; Example call (using defaults):
;; idl
;; .compile <path>/<to>/<script>/plot_chipsout_general.pro
;; plot_chipsout_general, 'my_output_tag', initials='NB', RTS=1, twoD=1, n_freq=384, band='high', lssa_num=0, $
;;   beam_point_weight = [0.,0.,0.,0.,1.,0.,0.,0.,0.]

pro plot_chipsout_general, output_tag, initials=initials, FHD=FHD, RTS=RTS, oneD=oneD, twoD=twoD, $
  n_freq=n_freq, band=band, base_freq=base_freq, chan_width=chan_width, lssa_num=lssa_num, $
  beam_point_weight = beam_point_weight, output_dir=output_dir


  ;Set default timing resolution
  if keyword_set(FHD) then begin
    deltat = 2.
    print, "Using FHD inputs (Default: RTS)"
  endif else begin
    deltat = 8.
    print, "Using RTS inputs (Default: RTS)"
  endelse

  ;Default to making 2D plots
  if (~keyword_set(oneD) and ~keyword_set(twoD)) then twoD=1

  ;Default to kriging applied
  if ~keyword_set(lssa_num) then begin
    lssa_num = '0'
    print, "Kriging has been applied to inputs (Default: kriging applied, lssa_num=0)"
  endif else print, "Kriging has not been applied to inputs (Default: kriging applied, lssa_num=0)"

  ;Default number of frequencies
  if ~keyword_set(n_freq) then begin
    if keyword_set(band) then begin
      if (band EQ 'ultra-low') then n_freq = 320
    endif else n_freq = 384
  endif
  print, "Number of frequency channels is " + strtrim(n_freq,2) + " (Default: 384, 320 for ultra-low)"
  Nchan = n_freq
  Nchanall = n_freq


  ;***** Strings and output file names
  if ~keyword_set(initials) then basename = '/fred/oz048/MWA/CODE/CHIPS/CHIPS_OUT/' $
    else basename = '/fred/oz048/MWA/CODE/CHIPS/' + STRUPCASE(initials) + '_CHIPS_OUT/'
  if ~keyword_set(output_dir) then spawn,"echo $HOME",output_dir
  output_dir = output_dir + '/'
  print, "Output plot directory is " + output_dir

  outputstring=output_tag

  if keyword_set(FHD) then begin
    instring1 = 'xx_'+ strtrim(lssa_num,2) +'.iter.' + outputstring
    instring2 = 'yy_'+ strtrim(lssa_num,2) +'.iter.' + outputstring
  endif else if keyword_set(RTS) then begin
    instring1 = 'yy_'+ strtrim(lssa_num,2) +'.iter.' + outputstring
    instring2 = 'xx_'+ strtrim(lssa_num,2) +'.iter.' + outputstring
  endif
  print, 'Polarizations are presented in FHD-format (N-S = YY, E-W = XX)'
  ;*****

  ;***** General obs/integration details
  if ( ~keyword_set(band) and ~keyword_set(base_freq) ) then band = 'high'

  if ( keyword_set(band) and ~keyword_set(base_freq) ) then begin
    if (band NE 'high') AND (band NE 'low') AND (band NE 'ultra-low') then $
      message, "Band name not recognized. Options are high, low, and ultra-low. Please use base_freq instead."
    print, 'Frequency band is ' + band + ' (Default: high)'
  endif

  if ( ~keyword_set(band) and keyword_set(base_freq) ) then print, 'Base frequency is ' + strtrim(base_freq,2) + ' Hz'

  if ~keyword_set(chan_width) then chan_width = 80.e3
  print, 'Frequency channel resolution is ' + strtrim(chan_width,2) + ' Hz (Default: 80000 Hz)'
  ;*****

  ;***** Freq-dependent parameters
  if keyword_set(band) then begin
    case band of
      'high': if keyword_set(FHD) then base_freq = 167.075e6 else if keyword_set(RTS) then base_freq = 167.035e6
      'low': if keyword_set(FHD) then base_freq = 138.915e6 else if keyword_set(RTS) then base_freq = 138.875e6
      'ultra-low': if keyword_set(FHD) then base_freq = 75.92e6 else if keyword_set(RTS) then base_freq = 75.52e6
    endcase
  endif

  ;Scale expected Tsys to match the center of the frequency band
  Tsys_high = 220. ;K for 182 MHz (scales automatically for low)
  Tsys = Tsys_high*( 182./round((base_freq + (chan_width*n_freq)/2.)/1e6) )^2.6 ;K
  print, 'Expected Tsys is ' + strtrim(Tsys,2)

  ;*****

  ;***** Beam volume calculation
  beamxx = read_binary('beam_area_point_xx_srhz.dat',data_type=4)
  if ~keyword_set(beam_point_weight) then begin
    beam_point_weight = [0.,0.,0.,0.,1.,0.,0.,0.,0.]
  endif
  print, 'Pointings used out of -4,-3,-2,-1,0,1,2,3,4: ' + $
    strjoin(strtrim(beam_point_weight,2),",") + ' (Default: just zenith)' 

  beam_area = total(beam_point_weight*beamxx)
  beam_area_per_chan = beam_area/Nchanall

  ; beam area calculated from int (beam^2) for each coarse channel
  beam_area_steradian = beam_area_per_chan/chan_width
  ;beam_area_steradian = 0.07597 ; -2 to +2 for BH gridding kernel

  obs_volume = mean(beam_area_steradian)*chan_width*Nchan    ; sr. Hz

  ;*****

  ;***** Plotting options
  mk_units=1
  plot_max_1D = 1.e11
  plot_min_1D = 1.e3
  plot_max_2D = 1.e15
  plot_min_2D = 1.e3
  ;*****

  low_k_bin_1D = 5e-3

  ;***** Cutting options
  wedge = 3   ; [0 = no cut; 3.0 = horizon]
  ;; 103.7 is the horizon
  wedge_angle = 120.  ;alternatevely, provide sky angle to cut. Will overwrite wedge
  kperp_min_1D = 18. ;in lambda
  kperp_max_1D = 80. ;in lambda
  kpar_min_1D = 0.105  ;in Mpc^-1
  kpar_max_1D = 10.
  ;*****

  ;***** determine size of binary files
  size_test = read_binary(basename+'crosspower_'+instring1+'.dat',data_type=4)
  nbins = ((size(size_test))[1]) / Nchan
  if nbins NE 80 then print, "WARNING: nbins does not match expected. Check n_freq"
  ;*****

  ;***** cosmology
  hubble = 0.704*100.  ; km/s/Mpc
  hubble_param = hubble / 100.
  c = 299792458. ; m/s
  omega_matter = 0.272
  omega_baryon = 0.046
  omega_lambda = 0.7
  bw = Nchan*chan_width   ; Hz
  f21 = c/0.21
  ;*****

  freq = dindgen(Nchan)*chan_width + base_freq
  z0_freq = 1420.4057517667e6 ;; Hz
  z = mean(z0_freq/freq -1d)
  cosmology_measures, z, Ez=Ez, comoving_dist_los = DM, omega_matter=omega_matter, omega_lambda=omega_lambda, $
    wedge_factor=wedge_factor, hubble_param=hubble_param
  if keyword_set(wedge_angle) then wedge = wedge_factor * wedge_angle * !pi/180.
  ;*****

  umax = 300.
  Neta = Nchan/2+1
  Neta = Nchan+1

  eta_actual = dindgen(Nchan)-Nchan/2.+1.
  loc = where(abs(eta_actual) lt 0.01)

  tvlct,[0,255,0,0],[0,0,255,0],[0,0,0,255]
  device,decomposed=0

  ; **** read-in data ****

  crosspower = read_binary(basename+'crosspower_'+instring1+'.dat',data_type=4,data_dims=[Nchan,nbins])
  totpower = read_binary(basename+'totpower_'+instring1+'.dat',data_type=4,data_dims=[Nchan,nbins])
  flagpower = read_binary(basename+'flagpower_'+instring1+'.dat',data_type=4,data_dims=[Nchan,nbins])
  residpower = read_binary(basename+'residpower_'+instring1+'.dat',data_type=4,data_dims=[Nchan,nbins])
  weights = read_binary(basename+'outputweights_'+instring1+'.dat',data_type=4,data_dims=[Nchan,nbins])
  fg_num = read_binary(basename+'fg_num_'+instring1+'.dat',data_type=3,data_dims=[Nchan,nbins])

  crosspower_2 = read_binary(basename+'crosspower_'+instring2+'.dat',data_type=4,data_dims=[Nchan,nbins])
  totpower_2 = read_binary(basename+'totpower_'+instring2+'.dat',data_type=4,data_dims=[Nchan,nbins])
  flagpower_2 = read_binary(basename+'flagpower_'+instring2+'.dat',data_type=4,data_dims=[Nchan,nbins])
  residpower_2 = read_binary(basename+'residpower_'+instring2+'.dat',data_type=4,data_dims=[Nchan,nbins])
  weights_2 = read_binary(basename+'outputweights_'+instring2+'.dat',data_type=4,data_dims=[Nchan,nbins])
  fg_num_2 = read_binary(basename+'fg_num_'+instring2+'.dat',data_type=3,data_dims=[Nchan,nbins])

  ; **** remove zeroed bins ****

  inds = where(total(crosspower, 1) NE 0,n_count)
  if n_count GT 0 then begin
    nbins = n_count
    crosspower = crosspower[*,inds]
    totpower = totpower[*,inds]
    flagpower = flagpower[*,inds]
    residpower = residpower[*,inds]
    weights = weights[*,inds]
    fg_num = fg_num[*,inds]
    crosspower_2 = crosspower_2[*,inds]
    totpower_2 = totpower_2[*,inds]
    flagpower_2 = flagpower_2[*,inds]
    residpower_2 = residpower_2[*,inds]
    weights_2 = weights_2[*,inds]
    fg_num_2 = fg_num_2[*,inds]
  endif

  ; **** start transposing and folding data ****

  crosspower = transpose(crosspower)
  totpower = transpose(totpower)
  flagpower = transpose(flagpower)
  residpower = transpose(residpower)
  weights = transpose(weights)
  fg_num = transpose(fg_num)

  crosspower_2 = transpose(crosspower_2)
  totpower_2 = transpose(totpower_2)
  flagpower_2 = transpose(flagpower_2)
  residpower_2 = transpose(residpower_2)
  weights_2 = transpose(weights_2)
  fg_num_2 = transpose(fg_num_2)

  crosspower_temp = fltarr(nbins,Nchan/2)
  totpower_temp = fltarr(nbins,Nchan/2)
  flagpower_temp = fltarr(nbins,Nchan/2)
  residpower_temp = fltarr(nbins,Nchan/2)
  weights_temp = fltarr(nbins,Nchan/2)
  fg_num_temp = fltarr(nbins,Nchan/2)

  for i=0,nbins-1 do begin

    crosspower_temp[i,0] = crosspower[i,loc]
    flagpower_temp[i,0] = flagpower[i,loc]
    residpower_temp[i,0] = residpower[i,loc]
    totpower_temp[i,0] = totpower[i,loc]
    weights_temp[i,0] = weights[i,loc]
    fg_num_temp[i,0] = fg_num[i,loc]


    for j=1,Nchan/2-1 do begin

      crosspower_temp[i,j] = (crosspower[i,j+loc])
      totpower_temp[i,j] = (totpower[i,j+loc])
      flagpower_temp[i,j] = (flagpower[i,j+loc])
      residpower_temp[i,j] = (residpower[i,j+loc])
      weights_temp[i,j] = (weights[i,j+loc])
      fg_num_temp[i,j] = (fg_num[i,j+loc])


    endfor
  endfor

  crosspower = crosspower_temp
  residpower = residpower_temp
  totpower = totpower_temp
  weights = weights_temp
  flagpower = flagpower_temp
  fg_num = fg_num_temp


  for i=0,nbins-1 do begin

    crosspower_temp[i,0] = crosspower_2[i,loc]
    flagpower_temp[i,0] = flagpower_2[i,loc]
    residpower_temp[i,0] = residpower_2[i,loc]
    totpower_temp[i,0] = totpower_2[i,loc]
    weights_temp[i,0] = weights_2[i,loc]
    fg_num_temp[i,0] = fg_num_2[i,loc]

    for j=1,Nchan/2-1 do begin

      crosspower_temp[i,j] = (crosspower_2[i,j+loc])
      totpower_temp[i,j] = (totpower_2[i,j+loc])
      flagpower_temp[i,j] = (flagpower_2[i,j+loc])
      residpower_temp[i,j] = (residpower_2[i,j+loc])
      weights_temp[i,j] = (weights_2[i,j+loc])
      fg_num_temp[i,j] = (fg_num_2[i,j+loc])

    endfor
  endfor

  crosspower_2 = temporary(crosspower_temp)
  residpower_2 = temporary(residpower_temp)
  totpower_2 = temporary(totpower_temp)
  weights_2 = temporary(weights_temp)
  flagpower_2 = temporary(flagpower_temp)
  fg_num_2 = temporary(fg_num_temp)


  ; **** end transposing and folding data ****


  ; *****************************

  ; Define parameters

  ; define cosmological units and conversions

  boltz = 1.38e-23
  D = 4.6  ;m
  lambda = c/freq   ;m

  jy2_Nchan__to__jy2hz2 = chan_width^2*Nchanall
  jy2__to__K2sr2 = mean(lambda)^4/(2.*boltz*1.e26)^2

  ;Via Morales & Hewitt 2004, using equations 2,3, and 4
  mpc_conversion = DM^2*(c*1e-3)*(1.+z)^2/(100.*hubble_param)/f21/Ez    ; Mpc^3/sr.Hz

  ; mK^2 Mpc^3 / Jy^2 N_f (1e6 to go from K^2 to mK^2)
  normalisation = jy2_Nchan__to__jy2hz2*jy2__to__K2sr2/obs_volume*1.e6*mpc_conversion

  M = Nchan/2

  ;FHD calibrates to (XX+YY) and RTS calibrates to (XX+YY)/2,
  ;; RTS visibilities are 2x larger for XX and YY
  ;If presenting (XX+YY) limits and PS, need to scale down
  ;; RTS power by 2^2.
  ;If presenting (XX+YY)/2 limits and PS, need to multiply
  ;; by a factor of 2^2 for FHD
  ;For both FHD and scaled RTS, multiply by 2 (4 for variances) to get
  ;; to estimate of Stokes I rather than instrumental pol
  if keyword_set(FHD) then cal_scale_factor=1. else cal_scale_factor=2.
  stokes_I_factor=2.
  noise_factor=2.
  print, 'Received polarized power (XX+YY convention)'
  ;if keyword_set(FHD) then cal_scale_factor=(2.) else cal_scale_factor=1.
  ;stokes_I_factor=1.
  ;noise_factor=1.
  ;print, 'Single polarized limits on an unpolarized signal (XX+YY/2 convention)'

  ;CHIPS expects a certain maximum weights. 32 is assumed as the max
  ;; Thus it comes in as a factor of ^2 in the normalization of the weights.
  ;; Does not affect the power.
  if keyword_set(FHD) then weights_factor=(2.) else weights_factor=1.

  ;CHIPS visibility weights normalization (set to 0.5 max)
  vis_weights_norm = 2.

  ;Crosspower can be constructed from (even^2 - odd^2) / 4. This needs to also be applied to 
  ;; the weights
  evenodd_weights_factor = 4.

  ;Density correction currently not determined to be needed.
  print, 'Applying density correction'
  den_corr = 0.5 ;for 0.5 lambda resolution, see Barry et al. 2019a

  Tsys_correction = 1./2. * Tsys^2. / 110.^2.
  normalisation_power = normalisation * stokes_I_factor / cal_scale_factor^2. / den_corr
  normalisation_weights = normalisation * stokes_I_factor / weights_factor / Tsys_correction $
    / vis_weights_norm / evenodd_weights_factor / sqrt(Nchan)

  expected_noise = noise_factor*boltz*Tsys*1.e26/D^2/sqrt(bw/Nchan*deltat)   ; Jy
  expected_noise = (expected_noise)^2   ; square ->Jy^2


  ; ************************************************

  pratio = crosspower*0.0
  flagratio = crosspower*0.0
  ptotratio = crosspower*0.0
  residratio = crosspower*0.0
  kkk = where(weights ne 0.)

  ; mK^2 Mpc^3
  pratio[kkk] = crosspower[kkk]/weights[kkk]*normalisation_power
  flagratio[kkk] = flagpower[kkk]/weights[kkk]*normalisation_power*expected_noise
  ptotratio[kkk] = totpower[kkk]/weights[kkk]*normalisation_power
  residratio[kkk] = residpower[kkk]/weights[kkk]*normalisation_power

  weight_data = weights/double(normalisation_weights)^2

  crosspower = pratio
  residpower = residratio
  totpower = ptotratio
  weights = weight_data
  flagpower = flagratio
  ; ********

  pratio = crosspower*0.0
  flagratio = crosspower*0.0
  ptotratio = crosspower*0.0
  residratio = crosspower*0.0
  kkk = where(weights_2 ne 0.)

  ; mK^2 Mpc^3
  pratio[kkk] = crosspower_2[kkk]/weights_2[kkk]*normalisation_power
  flagratio[kkk] = flagpower_2[kkk]/weights_2[kkk]*normalisation_power*expected_noise
  ptotratio[kkk] = totpower_2[kkk]/weights_2[kkk]*normalisation_power
  residratio[kkk] = residpower_2[kkk]/weights_2[kkk]*normalisation_power

  weight_data = weights_2/double(normalisation_weights)^2

  crosspower_2 = temporary(pratio)
  residpower_2 = temporary(residratio)
  totpower_2 = temporary(ptotratio)
  weights_2 = temporary(weight_data)
  flagpower_2 = temporary(flagratio)


  ; setup ratios and differences between the two inputs

  kgood = where(crosspower_2 ne 0.)

  ratio_cross = crosspower*0.

  ratio_cross[kgood] =  crosspower[kgood]/crosspower_2[kgood]
  diff_cross =  crosspower - crosspower_2

  ; ********************************************************************

  eta = findgen(Neta)-0.5   ; Using Bryna's definition for zeroth bin

  lperp = dindgen(nbins)*umax*1.1/float(nbins)*2.*!pi
  lperp[0] = lperp[1]/2.

  ;Via Morales & Hewitt 2004, eq 4
  eta_kz_factor = (c*1e-3)*(1.+z)^2/(2.*!pi*(100.*hubble_param)*f21*Ez) 
  kpa = eta/bw/eta_kz_factor
  kper = lperp/DM

  conv_factor = DM / (2.*!pi)
  ns_conv = eta/bw*1.e9

  ;; assume 20 degrees from pointing center to first null, and 103.7 for horizon
  horizon = wedge_factor * 103.7 * !dpi / 180d
  fov =  wedge_factor * 20d * !dpi / 180d

  ; 1D cuts, convert from lambda to Mpc^-1
  kperp_min_1D = kperp_min_1D / conv_factor
  kperp_max_1D = kperp_max_1D / conv_factor

  mx = plot_max_2D
  mn = plot_min_2D

  Nchancut=Nchan

  delay_params = [ns_conv[2]-ns_conv[1] , ns_conv[Nchancut/2-1]]

;  delay_delta = 1e9/(n_freq*chan_width) ;; equivilent delay bin size for kparallel
;  delay_max = delay_delta * n_freq/2.    ;; factor of 2 b/c of neg/positive
;  ;delay_params = [delay_delta, delay_max]

  ; ****************************
  ; ****************************

  kperp_plot_range=[7.,200.]/conv_factor ;in h Mpc^-1

  if keyword_set(twoD) then begin

    struc = {ncol:2,nrow:1,ordering:'col'}
    output = output_dir + 'plots_'+outputstring+'.png'

    kpower_2d_plots,kperp_edges=kper[1:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=(crosspower[1:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,title_prefix='E-W',$
      start_multi_params=struc,plotfile=output,/png,data_range=[plot_min_2D,plot_max_2D],/baseline_axis,kperp_lambda_conv=conv_factor,$
      /delay_axis,delay_params=delay_params,kperp_plot_range=kperp_plot_range

    kpower_2d_plots,kperp_edges=kper[1:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=(crosspower_2[1:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,$
      multi_pos=[0.5,0.,1.0,1.0],title_prefix='N-S',plotfile=output,/png,data_range=[plot_min_2D,plot_max_2D],/baseline_axis,$
      kperp_lambda_conv=conv_factor,/delay_axis,delay_params=delay_params,note=outputstring,$
      kperp_plot_range=kperp_plot_range

    cgps_close, /png, /delete_ps, density = 600

    ; ********************************


    struc = {ncol:2,nrow:1,ordering:'col'}
    output = output_dir + 'plots_'+outputstring+'_ratiodiff.png'

    mxdiff = max(abs(diff_cross))*10.

    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=ratio_cross[2:nbins-10,0:Nchancut/2-1],noise_expval=abs(flagpower[2:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,$
      wedge_amp=[fov,horizon],window_num=0,full_title='E-W/N-S Ratio',plotfile=output,/png,data_range=[0.1,10.],$
      start_multi_params=struc,delay_axis,delay_params=delay_params,/no_units,$
      kperp_plot_range=kperp_plot_range
    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=diff_cross[2:nbins-10,0:Nchancut/2-1],/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,full_title='E-W - N-S Diff',$
      plotfile=output,/png,data_range=[-mxdiff,mxdiff],multi_pos=[0.5,0.,1.0,1.0],/delay_axis,$
      delay_params=delay_params,weights=weights[2:nbins-10,0:Nchancut/2-1],color_profile='sym_log',note=outputstring,$
      kperp_plot_range=kperp_plot_range

    cgps_close, /png, /delete_ps, density = 600

    ; ********************************


    struc = {ncol:3,nrow:2,ordering:'col'}
    output = output_dir + 'plots_'+outputstring+'_noise.png'

    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=crosspower[2:nbins-10,0:Nchancut/2-1],/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,title_prefix='E-W',$
      /delay_axis,delay_params=delay_params,start_multi_params=struc,weights=weights[2:nbins-10,0:Nchancut/2-1],/plot_sigma,$
      plotfile=output,/png,data_range=[1.e4,1.e13],kperp_plot_range=kperp_plot_range
    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=(crosspower[2:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,$
      noise_expval=abs(flagpower[2:nbins-10,0:Nchancut/2-1]),/plot_exp_noise,title_prefix='E-W',multi_pos=[0.33,0.5,0.67,1],$
      plotfile=output,/png,data_range=[1.e4,1.e13],/delay_axis,delay_params=delay_params,$
      kperp_plot_range=kperp_plot_range
    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=crosspower[2:nbins-10,0:Nchancut/2-1],/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,title_prefix='E-W',$
      multi_pos=[0.67,0.5,1.0,1.0],noise_meas=abs(residpower[2:nbins-10,0:Nchancut/2-1]),/plot_noise,plotfile=output,/png,$
      data_range=[1.e4,1.e13],/delay_axis,delay_params=delay_params,kperp_plot_range=kperp_plot_range

    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=crosspower_2[2:nbins-10,0:Nchancut/2-1],/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,title_prefix='N-S',$
      multi_pos=[0.,0.,0.33,0.5],weights=weights_2[2:nbins-10,0:Nchancut/2-1],/plot_sigma,plotfile=output,/png,data_range=[1.e4,1.e13],$
      /delay_axis,delay_params=delay_params,kperp_plot_range=kperp_plot_range
    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=(crosspower_2[2:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,$
      noise_expval=abs(flagpower_2[2:nbins-10,0:Nchancut/2-1]),/plot_exp_noise,title_prefix='N-S',multi_pos=[0.33,0.,0.67,0.5],$
      plotfile=output,/png,data_range=[1.e4,1.e13],/delay_axis,delay_params=delay_params,$
      kperp_plot_range=kperp_plot_range
    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=crosspower_2[2:nbins-10,0:Nchancut/2-1],/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,title_prefix='N-S',$
      multi_pos=[0.67,0.,1.0,0.5],noise_meas=abs(residpower_2[2:nbins-10,0:Nchancut/2-1]),/plot_noise,plotfile=output,/png,$
      data_range=[1.e4,1.e13],/delay_axis,delay_params=delay_params,note=outputstring,$
      kperp_plot_range=kperp_plot_range
    cgps_close, /png, /delete_ps, density = 600

    ; ****************************
    struc = {ncol:2,nrow:1,ordering:'col'}
    output = output_dir + 'plots_'+outputstring+'_errors.png'
    err_data_range=[1.e4,1.e13]

    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=(crosspower[2:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,title_prefix='E-W',$
      start_multi_params=struc,plotfile=output,/png,data_range=err_data_range,/baseline_axis,kperp_lambda_conv=conv_factor,/delay_axis,$
      delay_params=delay_params,/plot_sigma,weights=weights[2:nbins-10,0:Nchancut/2-1],kperp_plot_range=kperp_plot_range

    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=(crosspower_2[2:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,multi_pos=[0.5,0.,1.0,1.0],$
      title_prefix='N-S',plotfile=output,/png,data_range=err_data_range,/baseline_axis,kperp_lambda_conv=conv_factor,/delay_axis,$
      delay_params=delay_params,/plot_sigma,weights=weights_2[2:nbins-10,0:Nchancut/2-1],note=outputstring,$
      kperp_plot_range=kperp_plot_range

    cgps_close, /png, /delete_ps, density = 600

    ; ********************

    struc = {ncol:2,nrow:2,ordering:'col'}
    output = output_dir + 'plots_'+outputstring+'_snr.png'

    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=crosspower[2:nbins-10,0:Nchancut/2-1],noise_expval=abs(flagpower[2:nbins-10,0:Nchancut/2-1]),$
      noise_meas=abs(residpower[2:nbins-10,0:Nchancut/2-1]),/nnr,/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,$
      title_prefix='E-W',plotfile=output,/png,data_range=[0.01,100.],start_multi_params=struc,/delay_axis,$
      delay_params=delay_params,kperp_plot_range=kperp_plot_range
    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=crosspower[2:nbins-10,0:Nchancut/2-1],noise_expval=abs(flagpower[2:nbins-10,0:Nchancut/2-1]),/snr,/plot_wedge_line,$
      wedge_amp=[fov,horizon],window_num=0,title_prefix='E-W',plotfile=output,/png,data_range=[1.e-2,1.e6],$
      multi_pos=[0.,0.,0.5,0.5],/delay_axis,delay_params=delay_params,weights=weights[2:nbins-10,0:Nchancut/2-1],$
      note=outputstring,kperp_plot_range=kperp_plot_range


    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=crosspower_2[2:nbins-10,0:Nchancut/2-1],noise_expval=abs(flagpower_2[2:nbins-10,0:Nchancut/2-1]),/snr,/plot_wedge_line,$
      wedge_amp=[fov,horizon],window_num=0,title_prefix='N-S',plotfile=output,/png,multi_pos=[0.5,0.,1,0.5],data_range=[1.e-2,1.e6],$
      /delay_axis,delay_params=delay_params,weights=weights_2[2:nbins-10,0:Nchancut/2-1],kperp_plot_range=kperp_plot_range
    kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=hubble_param,$
      power=crosspower_2[2:nbins-10,0:Nchancut/2-1],noise_expval=abs(flagpower_2[2:nbins-10,0:Nchancut/2-1]),$
      noise_meas=abs(residpower_2[2:nbins-10,0:Nchancut/2-1]),/nnr,/plot_wedge_line,wedge_amp=[fov,horizon],window_num=0,$
      title_prefix='N-S',multi_pos=[0.5,0.5,1,1],plotfile=output,/png,data_range=[0.01,100.],/delay_axis,$
      delay_params=delay_params,kperp_plot_range=kperp_plot_range

    cgps_close, /png, /delete_ps, density = 600

  endif

  ; ******************************************************************************************

  if keyword_set(oneD) then begin

    ; form 1D power uncertainty

    Neta = Nchan/2.
    Netaa = 100.

    kmax = sqrt(kper[nbins-1]^2 + kpa[Neta-1]^2)

    ptot_xx = dblarr(Netaa)
    pmeas_xx = dblarr(Netaa)
    ptot_yy=dblarr(Netaa)
    pmeas_yy=dblarr(Netaa)

    num = intarr(Netaa)
    numedge = fltarr(Netaa)
    valedge = fltarr(Netaa)
    mask = weights*0.0

    ;ktot_bins = (dindgen(Netaa+1)/float(Netaa))^3.*kmax+0.05
    ;ktot_bins = (dindgen(Netaa+1)/float(Netaa))^1.55*kmax+0.03

    ktot_bins = (dindgen(Netaa+1)/float(Netaa))^1.4*kmax+low_k_bin_1D
    ;ktot_bins = findgen(Netaa+1)*0.0288 + low_k_bin_1D     ; nichole settings
    noise_obs_xx = sqrt(weights)
    noise_obs_yy = sqrt(weights_2)

    edge = dindgen(nbins)/nbins*5.    ; 2 sigma correlation length is 5 kpar bins

    for i = 1 , nbins-1 do begin
      for j = 0 , Neta-1 do begin

        for k = 0 , Netaa-1 do begin
          if ((sqrt(kper[i]^2 + kpa[j]^2) ge ktot_bins[k]) and (sqrt(kper[i]^2 + kpa[j]^2) lt ktot_bins[k+1])) then begin
            if ((kpa[j] gt kper[i]*wedge) and (kper[i] lt kperp_max_1D) and (kper[i] ge kperp_min_1D) and (kpa[j] gt kpar_min_1D)) then begin

              ptot_xx[k] = ptot_xx[k] + (noise_obs_xx[i,j])^2
              ptot_yy[k] = ptot_yy[k] + (noise_obs_yy[i,j])^2

              pmeas_xx[k] = pmeas_xx[k] + (crosspower[i,j])*(noise_obs_xx[i,j])^2
              pmeas_yy[k] = pmeas_yy[k] + (crosspower_2[i,j])*(noise_obs_yy[i,j])^2
              num[k] = num[k] + float(fg_num[i,j])/Nchan
              numedge[k] = numedge[k]+1.
              valedge[k] = valedge[k]+edge[i]
              mask[i,j] = k
            endif
          endif
        endfor

      endfor
    endfor

    ; find the average value of the +ve x error bar for the cells that enter each ktot_bin

    xerrhi = fltarr(Netaa)
    fin=where(numedge gt 0.)
    xerrhi[fin] = valedge[fin]/numedge[fin]*(kpa[5]-kpa[4])*2.

    xerrlo = dblarr(Netaa)+(kpa[5]-kpa[4])    ; for 2 sigma contained within that cell

    ; ptot here is the uncertainty on the dimensionless power in mK^2
    pmeas_xx = pmeas_xx/ptot_xx
    pmeas_yy=pmeas_yy/ptot_yy
    ptot_xx = 1./sqrt(ptot_xx)
    ptot_yy = 1./sqrt(ptot_yy)

    if keyword_set(mK_units) then begin
      pmeas_xx = pmeas_xx*ktot_bins^3/2./!pi^2
      pmeas_yy = pmeas_yy*ktot_bins^3/2./!pi^2
      ptot_xx = (ptot_xx)*(ktot_bins^3/2./!pi^2)
      ptot_yy = (ptot_yy)*(ktot_bins^3/2./!pi^2)
      y_title = 'log!D10!N $\Delta$$\exp2$ (k) (mK!E2!N)'
      full_title = '$\Delta$$\exp2$(k) = (k!E3!NP(k)/2!7p!3!E2!N) (mK!E2!N)'
      csv_units = 'Delta^2'
    endif else begin
      pmeas_xx = pmeas_xx * hubble_param^3.
      pmeas_yy = pmeas_yy * hubble_param^3.
      ptot_xx = ptot_xx * hubble_param^3.
      ptot_yy = ptot_yy * hubble_param^3.
      y_title = 'P (mK$\exp2$ Mpc$\exp3$ /  !8h!x$\exp3$)'
      full_title = 'P(k)'
      csv_units = 'mK^2 Mpc^3 / h^3'
    endelse

    ;for h Mpc^-1 units in k
    ktot_bins = ktot_bins / hubble_param

    sigma = 2. ;change for 1sigma or 2sigma errors
    ptot_xx = ptot_xx*sigma
    ptot_yy = ptot_yy*sigma

    output = output_dir + 'plots_'+outputstring+'_1D.png'
    cgPS_Open,output,/nomatch

    bin_start=1

    ;Write CSV file with all 1d binned values
    WRITE_CSV, output_dir + '1D_text_values_'+outputstring+'.csv', ktot_bins[bin_start:Netaa-1], pmeas_xx[bin_start:Netaa-1],$
      pmeas_yy[bin_start:Netaa-1],xerrhi[bin_start:Netaa-1],xerrlo[bin_start:Netaa-1],ptot_xx[bin_start:Netaa-1],$
      ptot_yy[bin_start:Netaa-1],header=['k [h Mpc^-1]','P XX ['+csv_units+']','P YY ['+csv_units+']','xerrhigh [h Mpc^-1]',$
      'xerrlow [h Mpc^-1]','thermal P XX ['+csv_units+']','thermal P YY ['+csv_units+']']

    ;fix error bars to edge of graph window for aesthetics
    xerrhi_plot = xerrhi
    xerrlo_plot = xerrlo
    ptot_yy_min_plot = ptot_yy
    ptot_xx_min_plot = ptot_xx
    for bin_i=bin_start,Netaa-1 do begin
      inds = where(xerrhi[bin_i] GT 4 - ktot_bins[bin_i],n_count)
      if n_count GT 0 then xerrhi_plot[bin_i] = 4 - ktot_bins[bin_i]
      inds = where(xerrlo[bin_i] GT ktot_bins[bin_i] - low_k_bin_1D,n_count)
      if n_count GT 0 then xerrlo_plot[bin_i] = ktot_bins[bin_i] - low_k_bin_1D
      inds = where(ptot_yy[bin_i] GT pmeas_yy[bin_i] - plot_min_1D,n_count)
      if n_count GT 0 then ptot_yy_min_plot[bin_i] = pmeas_yy[bin_i] - plot_min_1D
      inds = where(ptot_xx[bin_i] GT pmeas_xx[bin_i] - plot_min_1D,n_count)
      if n_count GT 0 then ptot_xx_min_plot[bin_i] = pmeas_xx[bin_i] - plot_min_1D
    endfor

    cgplot,[0,0],[0,0],charthick=1.75,thick=2.5,title='1D Power - ' + full_title, $
      xrange=[low_k_bin_1D,10],ystyle=1,yrange=[plot_min_1D,plot_max_1D],xtitle='log!D10!N k (h Mpc!E-1!N)',ytitle=y_title,$
      /ylog,/xlog,xstyle=1,charsize=1,/nodata
    cgoplot,ktot_bins[bin_start:Netaa-1],pmeas_xx[bin_start:Netaa-1],color='red',psym=1,thick=2,linestyle=6,$
      ERR_XLow=xerrlo_plot[bin_start:Netaa-1], ERR_XHigh=xerrhi_plot[bin_start:Netaa-1],ERR_YLow=ptot_xx_min_plot[bin_start:Netaa-1],$
      ERR_YHigh=ptot_xx[bin_start:Netaa-1], ERR_Color='red'
    cgoplot,ktot_bins[bin_start:Netaa-1],pmeas_yy[bin_start:Netaa-1],color='blue',thick=2,psym=1,linestyle=6,$
      ERR_XLow=xerrlo_plot[bin_start:Netaa-1], ERR_XHigh=xerrhi_plot[bin_start:Netaa-1], ERR_YLow=ptot_yy_min_plot[bin_start:Netaa-1],$
      ERR_YHigh=ptot_yy[bin_start:Netaa-1], ERR_Color='blue'

    cglegend, title=['E-W','N-S'], color=['red','blue'], location=[.2,.8], charsize=1, length=.05
    cgoplot,ktot_bins[bin_start:Netaa-1],ptot_xx[bin_start:Netaa-1]/2.,color='green',thick=2.5
    cgoplot,ktot_bins[bin_start:Netaa-1],ptot_yy[bin_start:Netaa-1]/2.,color='green',thick=2.5

    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage

  endif

end
