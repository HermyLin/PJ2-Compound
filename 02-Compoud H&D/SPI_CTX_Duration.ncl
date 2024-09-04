load "/lfs/home/killumy/functions/func_PJ.ncl"

begin

ispv = -999
fspv = 1e+36
year = ispan(1980,2019,1)
nyear = dimsizes(year)

year_po = (/1980, 1983, 1993, 1998, 1999, 2003, 2014, 2015/)
year_ne = (/1981, 1984, 1985, 1990, 1994, 1997, 2001, 2013, 2018/)
DOY_str = 152
DOY_end = 243

year1d = array_append_record(year_po,year_ne,0)
qsort(year1d)
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
  intensity@_FillValue = 0.0

  do yy = 0, nyear-1
     data = tmax_mask(yy,:)
     intensity(yy,:) = where(data.gt.CTX_mask,1.0,0.0)
     delete(data)
  end do
  ;print(intensity(0,:))
  ;print(num(intensity(0,:).eq.1))
  temp_data = toint(intensity)
;===================================================================================
;count CTX
;===================================================================================
  runs = dim_numrun_n(temp_data,0,1) ; calculate runs

  CTX_Duration = new(nyear,"integer",ispv)
  CTX_Duration!0 = "year"
  CTX_Duration&year = year

  do yy = 0, nyear-1
     run_yby = runs(yy,:)
     indmx = ind(run_yby.ne.0)
     mxrun = indmx(dimsizes(indmx)-1) + 1
     ;print(year(yy)+": "+mxrun)

     CTX_Duration(yy) = mxrun

     delete([/run_yby,indmx,mxrun/])
  end do
  delete(intensity)

  x_po = CTX_Duration({year_po})
  x_ne = CTX_Duration({year_ne})

  x_nu = CTX_Duration
  x_nu({year1d}) = ispv

  delete(CTX_Duration)
;==================================================
;rain
;==================================================   
print("*****YAXIS*****")
  print("---read SPI---")
  diri = "/lfs/home/killumy/SYLIN_DATABASE/TCCIP_Grid/5km/"
  fn   = "Monthly_rain.nc"
  f = addfile(diri+fn,"r")

  prep = f->prep({year},:,:,:)
  prep!0 = "time"
  prep@_FillValue = fspv
  prep@missing_value = fspv

  ndim = dimsizes(prep)
  ntim = ndim(0)*ndim(1)

  lat = f->lat
  lon = f->lon
  nlat = dimsizes(f->lat)
  nlon = dimsizes(f->lon)

  prep_re = reshape(prep,(/ntim,nlat,nlon/))
  copy_VarCoords(prep(0,0,:,:),prep_re(0,:,:))
  prep := prep_re
  delete(prep_re)

  ;mask mountain site
  H = 1000
  height = get_Taiwan_height(fspv,lat,lon)
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  do i = 0, ntim-1
     prep(i,:,:) = prep(i,:,:)*height
  end do

  delete(height)
  ;-----------------------------------------------
  ;calculate SPI
  ;-----------------------------------------------
  run = 3
  opt = False
  opt@spi_type = 3             ; 6.3.0 onward; only if opt=True

  spi = dim_spi_n(prep,run,opt,0)
  copy_VarCoords(prep,spi)

  ;JJA
  spi_year = reshape(spi,ndim)
  spi_JJA := dim_avg_n(dim_avg_n(spi_year(:,5:7,:,:),1),(/1,2/))

  spi_JJA!0 = "time"
  spi_JJA&time = year
  delete([/spi,spi_year/])

  spi_po = spi_JJA({year_po})
  spi_ne = spi_JJA({year_ne})

  ;nuetral
  spi_JJA({year1d}) = fspv
  spi_nu = spi_JJA

  delete(spi_JJA)

  y_po = spi_po
  y_ne = spi_ne
  y_nu = spi_nu
  delete([/spi_po,spi_ne,spi_nu/])
;===================================================================================
;regression
;===================================================================================
  ;x = array_append_record(ndtooned(temp_po),ndtooned(temp_ne),0)
  ;y = array_append_record(ndtooned(spi_po),ndtooned(spi_ne),0)
  
  ;a = regline_stats(x,y)
  ;r = escorc(x,y)
  ;print(r)
  ;n = 24
  ;df = n-2
  ;t = r*sqrt((n-2)/(1-r^2))      
  ;p = student_t(t, df)
  ;print(""+p)
;==================================================
;plot
;frequency, HWV, HWMA, maxHWA, HWMD, HWMS, avgHWMS, HWMS
;==================================================
print("---plot---")
  fig_type = "X11"
  fig_type = "eps"
  fig_name = "240725_Scatter_CTXduration_SPI_Neutral"

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

  ;res@tmYLMode = "Explicit"
  ;res@tmYLValues = fspan(30,35,11)
  ;res@tmYLLabels = "" + res@tmYLValues
  ;res@tmYLMinorValues = fspan(30,35,51)

  res@trXMinF = -2.2
  res@trXMaxF = 2.0
  res@trYMinF = 0
  res@trYMaxF = 12

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

  pres = True
  pres@gsLineColor = "gray40"
  pres@gsLineDashPattern = 1
  pres@gsLineThicknessF = 2.0
  dum = gsn_add_polyline(wks,plot_nu,(/0,0/),(/0,100/), pres)

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
