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

  CTX_days = new(nyear,"integer",ispv)

  do yy = 0, nyear-1
     CTX_days(yy) = num(tmax_mask(yy,:).gt.CTX_mask)
     ;print(year(yy)+": "+ CTX_days(yy) )
  end do
;----------------------------------------------------------
  CTX_days!0 = "year"
  CTX_days&year = year

  CTX_po = CTX_days({year_po})
  CTX_ne = CTX_days({year_ne})

  ;nuetral
  CTX_days({year_PJ}) = ispv
  CTX_nu = CTX_days

  delete([/f_tmax,CTX_days/])
;===================================================================================
;read soil data
;===================================================================================
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
  spi@_FillValue = fspv

  ;JJA
  spi_year = reshape(spi,ndim)
  spi_JJA := dim_avg_n(dim_avg_n(spi_year(:,5:7,:,:),1),(/1,2/))

  spi_JJA!0 = "time"
  spi_JJA&time = year
  delete([/spi,spi_year/])

  spi_po = spi_JJA({year_po})
  spi_ne = spi_JJA({year_ne})

  ;nuetral
  spi_JJA({year_PJ}) = fspv
  spi_nu = spi_JJA

  delete(spi_JJA)
;===================================================================================
;plot
;===================================================================================
print("---plot---")
  fig_type = "X11"
  fig_type = "eps"

  outdir   = "./"
  fig_name = outdir + "240725_Scatter_CTX_SPI_Neutral"
  plot = new(2,graphic)

  wks = gsn_open_wks(fig_type,fig_name)

  res  = True
  res@vpWidthF  = 0.6
  res@vpHeightF = 0.4

  res@gsnDraw  = False
  res@gsnFrame = False

  res@tiYAxisFont = 21
  res@tiYAxisFontHeightF = 0.025
  res@tiYAxisString = ""

  res@tiXAxisFont = 21
  res@tiXAxisFontHeightF = 0.025
  res@tiXAxisString = ""

  ;res@tmXBLabelFont = 22
  res@tmXBLabelFontHeightF = 0.025
  ;res@tmYLLabelFont = 22
  res@tmYLLabelFontHeightF = 0.025

  res@trYMinF = 0
  res@trYMaxF = 30
  res@trXMinF = -2.2
  res@trXMaxF = 2.0
  printMinMax(spi_nu,0)

  res@xyMarkLineMode  = "Markers"  ; choose which have markers
  res@xyMarker        = 1         ; choose type of marker 
  res@xyMarkerColor   = "red"      ; Marker color
  res@xyMarkerSizeF = 0.08       ; Marker size (default 0.01)
  plot_po = gsn_csm_xy(wks,spi_po,CTX_po,res)

  res@xyMarkerColor = "blue"                  ; Marker color
  plot_ne = gsn_csm_xy(wks,spi_ne,CTX_ne,res)

  res@xyMarkerSizeF = 0.015       ; Marker size (default 0.01)
  res@xyMarkerThicknessF = 3.0
  res@xyMarker = 4         ; choose type of marker 
  res@xyMarkerColor = "black"    
  plot_nu = gsn_csm_xy(wks,spi_nu,CTX_nu,res)

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
