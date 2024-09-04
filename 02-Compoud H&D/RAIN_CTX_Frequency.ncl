load "/lfs/home/killumy/functions/func_PJ.ncl"
;===================================================================================
;------------------------main code------------------------
;===================================================================================
begin

;begTime = get_cpu_time()
fspv = 1e+36
ispv = -999
year = ispan(1980,2019,1)
nyear = dimsizes(year)

year_po = (/1980, 1983, 1993, 1998, 1999, 2003, 2014, 2015/)
year_ne = (/1981, 1984, 1985, 1990, 1994, 1997, 2001, 2013, 2018/)

year_PJ = array_append_record(year_po,year_ne,0)
qsort(year_PJ)
;===================================================================================
;read temp data
;===================================================================================
print("*****XAXIS*****")
  print("---read CTX---")
  diri = "/lfs/home/killumy/01-PJ_OBS/0data/"
  fn   = "CTX_variables_mask1000.nc"
  f = addfile(diri+fn,"r")

  CTX = int2flt(f->CTX_month_number({year},5:7,:,:))
  CTX@_FillValue = fspv

  CTX_number = dim_avg_n_Wrap(CTX,1)
  CTX_number = mask(CTX_number,CTX_number.eq.0,False)

  ;JJA mean
  CTX_po = dim_avg_n_Wrap(CTX_number({year_po},:,:),(/1,2/))
  CTX_ne = dim_avg_n_Wrap(CTX_number({year_ne},:,:),(/1,2/))

  ;nuetral
  CTX_number({year_PJ},:,:) = fspv
  CTX_nu = dim_avg_n_Wrap(CTX_number,(/1,2/))

  delete([/CTX,CTX_number/])
  delete([/diri,fn,f/])
;===================================================================================
;read soil data
;===================================================================================
print("*****YAXIS*****")
  print("---read soil---")
  diri = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/Monthly/"
  fn   = "TReAD_2km_Monthly_PRECTmms.nc"
  f = addfile(diri+fn,"r")

  data = f->PRECTmms
  data@_FillValue = fspv
  data@missing_value = fspv

  data = data*86400    ;mm/day
  data_JJA = dim_avg_n_Wrap(data({year},5:7,:,:),1)
  data := data_JJA
  delete(data_JJA)

  lat = f->Lat
  lon = f->Lon
  ;------------------------------------------
  ;mask height with 2km data
  ;------------------------------------------
  H = 1000
  height = get_Taiwan_height(fspv,lat,lon)
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  do i = 0, nyear-1
     data(i,:,:) = data(i,:,:) * height
  end do

  soil_po = dim_avg_n_Wrap(data({year_po},:,:),(/1,2/))
  soil_ne = dim_avg_n_Wrap(data({year_ne},:,:),(/1,2/))

  data({year_PJ},:,:) = fspv
  soil_nu = dim_avg_n_Wrap(data,(/1,2/))
  ;printMinMax(soil_nu,0)

  delete([/diri,fn,f,data/])
;===================================================================================
;plot
;===================================================================================
print("---plot---")
  fig_type = "X11"
  fig_type = "eps"

  outdir   = "./"
  fig_name = outdir + "TReAD_Scatter_CTX_Rain_Neutral"

  wks = gsn_open_wks(fig_type,fig_name)

  res  = True
  res@vpWidthF  = 0.6
  res@vpHeightF = 0.4

  res@gsnDraw  = False
  res@gsnFrame = False

  res@xyMarkLineMode  = "Markers"  ; choose which have markers
  res@xyMarker        = 1         ; choose type of marker 
  res@xyMarkerColor   = "red"      ; Marker color
  ;res@xyDashPattern    = 1         ; solid line 
  ;res@xyLineThicknessF = 3.0       ; set second line to 2
  
  res@tiYAxisFont = 21
  res@tiYAxisFontHeightF = 0.025
  res@tiYAxisString = "day/month"

  res@tiXAxisFont = 21
  res@tiXAxisFontHeightF = 0.025
  res@tiXAxisString = "mm"

  ;res@tmXBLabelFont = 22
  res@tmXBLabelFontHeightF = 0.025
  ;res@tmYLLabelFont = 22
  res@tmYLLabelFontHeightF = 0.025

  ;res@tmYLMode = "Explicit"
  ;res@tmYLValues = fspan(0.31,0.51,11)
  ;res@tmYLLabels = "" + res@tmYLValues
  ;res@tmYLMinorValues = fspan(0.3,0.5,21)

  ;res@tmXBMode = "Explicit"
  ;res@tmXBValues = ispan(0,30,2)
  ;res@tmXBLabels = "" + res@tmXBValues
  ;res@tmXBMinorValues = ispan(0,30,1)

  res@trYMinF = 0
  res@trYMaxF = 12
  res@trXMinF = 4
  res@trXMaxF = 20

  res@xyMarkerSizeF = 0.08       ; Marker size (default 0.01)
  plot_po = gsn_csm_xy(wks,soil_po,CTX_po,res)

  res@xyMarkerColor = "blue"                  ; Marker color
  plot_ne = gsn_csm_xy(wks,soil_ne,CTX_ne,res)

  res@xyMarkerSizeF = 0.015       ; Marker size (default 0.01)
  res@xyMarkerThicknessF = 3.0
  res@xyMarker = 4         ; choose type of marker 
  res@xyMarkerColor = "black"                  ; Marker color
  plot_nu = gsn_csm_xy(wks,soil_nu,CTX_nu,res)
  
  overlay(plot_nu,plot_po)
  overlay(plot_nu,plot_ne)

  draw(plot_nu)
  frame(wks)

  delete([/CTX_po,CTX_ne,soil_po,soil_ne/])
;---change figure type----------------------------------------------------------

  delete(wks)
  if (fig_type.eq."eps") then
  system("convert -trim -density 300 -crop 0x0 "+fig_name+".eps  "+fig_name+".png")
         ;==convert + quality(resolution) up + remove white boundary==;

  system("/bin/rm "+fig_name+".eps")
  end if
  print("Output Figure Name = "+fig_name+".png")

end
