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
;===================================================================================
;read soil data
;===================================================================================
print("*****YAXIS*****")
  print("---read soil---")
  diri = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/OneF_Monthly/"
  fn   = "TReAD_2km_Monthly_PRECTmms.nc"
  f = addfile(diri+fn,"r")

  data = f->PRECTmms
  data@_FillValue = fspv
  data@missing_value = fspv

  data = data*86400    ;mm/day
  data_JJA = dim_avg_n_Wrap(data({year},5:7,:,:),1)
  data := data_JJA

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
  copy_VarCoords(data_JJA,data)

  y_po = dim_avg_n_Wrap(data({year_po},:,:),(/1,2/))
  y_ne = dim_avg_n_Wrap(data({year_ne},:,:),(/1,2/))

  data({year1d},:,:) = fspv
  y_nu = dim_avg_n_Wrap(data,(/1,2/))

  delete([/diri,fn,f,data_JJA/])
;==================================================
;regression
;==================================================
  x = dim_avg_n_Wrap(temp({year1d},:,:),(/1,2/))
  y = dim_avg_n_Wrap(data({year1d},:,:),(/1,2/))
  delete([/temp,data/])

  r = escorc(x,y)
  n = dimsizes(x)
  df = n-2
  t = r*sqrt((n-2)/(1-r^2))      
  p = student_t(t, df)
  ;print("R = "+r)
  ;print("P = "+p)
  delete([/n,df,t,p/])
  ;exit
;==================================================
;plot
;frequency, HWV, HWMA, maxHWA, HWMD, HWMS, avgHWMS, HWMS
;==================================================
print("---plot---")
  fig_type = "X11"
  fig_type = "eps"
  
  outdir   = "./"
  fig_name = "240703_Scatter_TReADRain_TCCIPCTXtmax"  
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

  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(5,20,5)
  res@tmXBLabels = "" + res@tmXBValues
  res@tmXBMinorValues = ispan(5,20,1)

  res@tmYLMode = "Explicit"
  res@tmYLValues = fspan(30,35,11)
  res@tmYLLabels = "" + res@tmYLValues
  res@tmYLMinorValues = fspan(30,35,51)

  res@trXMinF = 5
  res@trXMaxF = 20
  res@trYMinF = 32.8
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
