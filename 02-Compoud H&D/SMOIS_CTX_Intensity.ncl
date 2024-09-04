load "/lfs/home/killumy/functions/func_PJ.ncl"

begin

ispv = -999
fspv = 1e+36
year = ispan(1980,2019,1)
nyear = dimsizes(year)

year_po = (/1980, 1983, 1993, 1998, 1999, 2003, 2014, 2015/)
year_ne = (/1981, 1984, 1985, 1990, 1994, 1997, 2001, 2013, 2018/)

year1d = array_append_record(year_po,year_ne,0)
qsort(year1d)
;==================================================
;read HWV
;==================================================
print("*****XAXIS*****")
print("---read---")
/;
  diri = "/lfs/home/killumy/01-PJ_OBS/0data/"
  fn   = "CTX_variables_mask1000.nc"
  f = addfile(diri+fn,"r")

  JJA_DOY = ispan(152,243,1)
  temp = dim_avg_n_Wrap(f->CTX_pdf({year},JJA_DOY,:,:),1)

  x_po = dim_avg_n_Wrap(temp({year_po},:,:),(/1,2/))
  x_ne = dim_avg_n_Wrap(temp({year_ne},:,:),(/1,2/))

  temp({year1d},:,:) = fspv
  x_nu = dim_avg_n_Wrap(temp,(/1,2/))
  delete([/diri,fn,f/])
;/
;===================================================================================
;read temp data
;===================================================================================
print("*****XAXIS*****")
  print("---read CTX---")
  diri = "/lfs/home/killumy/01-PJ_OBS/0data/"
  ;fn   = "CTX_variables_mask1000.nc"
  fn = "00_CTXraw_umasked.nc"
  f_CTX = addfile(diri+fn,"r")

  CTX = f_CTX->CTXp90(152:243,:,:)
  CTX := dble2flt(CTX)
  CTX@_FillValue = fspv
  replace_ieeenan(CTX, fspv, 0)
  ;------------------------------------------
  ;mask height with 5km data
  ;------------------------------------------
  lat = dble2flt(f_CTX->lati)
  lon = dble2flt(f_CTX->long)
  H = 1000
  height = get_Taiwan_height(fspv,lat,lon)
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  scalarH = conform(CTX,height,(/1,2/))
  CTX = CTX * scalarH
  ;------------------------------------------

  CTX_mask = dim_avg_n(CTX,(/1,2/))
  delete([/diri,fn,f_CTX,CTX,scalarH/])
;===================================================================================
;read temp data
;===================================================================================
  f_tmax = addfile("/lfs/home/killumy/SYLIN_DATABASE/TCCIP_Grid/5km/Daily_tmax.nc","r")
  tmax = f_tmax->tmax({1980:2019},151:242,:,:)
  tmax@_FillValue = fspv

  scalarH = conform(tmax,height,(/2,3/))
  tmax = tmax*scalarH
  tmax_mask = dim_avg_n_Wrap(tmax,(/2,3/))
  delete([/lat,lon,height,scalarH,tmax/])

  intensity = tmax_mask


  do yy = 0, nyear-1
     data = tmax_mask(yy,:)
     intensity(yy,:) = where(data.gt.CTX_mask,data,fspv)
     delete(data)
  end do
;----------------------------------------------------------
  intensity!0 = "year"
  intensity&year = year

  x_po = dim_avg_n(intensity({year_po},:),1)
  x_ne = dim_avg_n(intensity({year_ne},:),1)

  ;nuetral
  intensity({year1d},:) = ispv
  x_nu = dim_avg_n(intensity,1)

  delete([/f_tmax,intensity/])
;===================================================================================
;read soil data
;===================================================================================
print("*****YAXIS*****")
  print("---read soil---")
  diri = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/OneF_Monthly/"
  fn   = "TReAD_2km_Monthly_SMOIS.nc"
  f = addfile(diri+fn,"r")

  ;soil@_Fillvalue = fspv
  ;soil@missing_value = fspv

  S1  = dim_avg_n_Wrap(f->Soil1(:,6:8,:,:),1)
  S2  = dim_avg_n_Wrap(f->Soil1(:,6:8,:,:),1)

  data = ((S1*0.07) + (S2*0.21)) / 0.28
  copy_VarCoords(S1,data)

  delete([/S1,S2/])

  lat = f->Lat
  lon = f->Lon
  nyear = dimsizes(f->Year)
  ;------------------------------------------
  ;mask height with 2km data
  ;------------------------------------------
  H = 1000
  height = get_Taiwan_height(fspv,lat,lon)
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  do yy = 0, nyear-1
     data(yy,:,:) = data(yy,:,:) * height
  end do

  y_po = dim_avg_n_Wrap(data({year_po},:,:),(/1,2/))
  y_ne = dim_avg_n_Wrap(data({year_ne},:,:),(/1,2/))

  data({year1d},:,:) = fspv
  y_nu = dim_avg_n_Wrap(data({1980:2019},:,:),(/1,2/))

  delete([/diri,fn,f/])
;==================================================
;plot
;frequency, HWV, HWMA, maxHWA, HWMD, HWMS, avgHWMS, HWMS
;==================================================
print("---plot---")
  fig_type = "X11"
  fig_type = "eps"
  
  outdir   = "./"
  fig_name = outdir + "240725_Scatter_CTXtmax_SMOIS_Neutral"  
  wks = gsn_open_wks(fig_type,fig_name)

  res  = True
  res@vpWidthF  = 0.6
  res@vpHeightF = 0.4

  res@gsnDraw  = False
  res@gsnFrame = False

  res@xyMarkLineMode  = "Markers"  ; choose which have markers
  res@xyMarker        = 1         ; choose type of marker 
  res@xyMarkerColor   = "red"      ; Marker color
  
  res@tiXAxisFont = 21
  res@tiXAxisFontHeightF = 0.025
  ;res@tiXAxisString = "HWV (10^~S~3~N~)"
  res@tiYAxisFont = 21
  res@tiYAxisFontHeightF = 0.025
  res@tiYAxisString = ""

  ;res@tmXBLabelFont = 22
  res@tmXBLabelFontHeightF = 0.025
  ;res@tmYLLabelFont = 22
  res@tmYLLabelFontHeightF = 0.025
  res@tmXBFormat = "0@;*.2f"

  res@tmYLMode = "Explicit"
  res@tmYLValues = fspan(30,35,11)
  res@tmYLLabels = "" + res@tmYLValues
  res@tmYLMinorValues = fspan(30,35,51)

  res@trXMinF = 0.28
  res@trXMaxF = 0.40
  res@trYMinF = 32
  res@trYMaxF = 34

  res@tiMainString  = ""

  res@xyMarkerSizeF = 0.08       ; Marker size (default 0.01)
  plot_po = gsn_csm_xy(wks,y_po,x_po,res)

  res@xyMarkerColor = "blue"                  ; Marker color
  plot_ne = gsn_csm_xy(wks,y_ne,x_ne,res)

  res@xyMarkerSizeF = 0.015       ; Marker size (default 0.01)
  res@xyMarkerThicknessF = 3.0
  res@xyMarker = 4         ; choose type of marker
  res@xyMarkerColor = "black"                  ; Marker color
  plot_nu = gsn_csm_xy(wks,y_nu,x_nu,res)

  pres = True
  pres@gsLineColor = "gray40"
  pres@gsLineDashPattern = 1
  pres@gsLineThicknessF = 2.0
  dum = gsn_add_polyline(wks,plot_nu,(/0.33,0.33/),(/0,100/), pres)

  overlay(plot_nu,plot_po)
  overlay(plot_nu,plot_ne)

  draw(plot_nu)
  frame(wks)
;---change figure type----------------------------------------------------------

  delete(wks)
  if (fig_type.eq."eps") then
  system("convert -trim -density 300 -crop 0x0 "+fig_name+".eps  "+fig_name+".png")
         ;==convert + quality(resolution) up + remove white boundary==;

  system("/bin/rm "+fig_name+".eps")
  end if
  print("Output Figure Name = "+fig_name+".png")
end
