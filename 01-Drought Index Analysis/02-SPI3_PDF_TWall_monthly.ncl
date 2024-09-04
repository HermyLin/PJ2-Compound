load "/lfs/home/killumy/functions//func_PJ.ncl"
;===================================================================================
;------------------------main code------------------------
;===================================================================================
begin

;begTime = get_cpu_time()
fspv = 1e+20
ispv = -999
year = ispan(1980,2019,1)
nyear = dimsizes(year)
month = ispan(1,12,1)
nmonth = dimsizes(month)

month_day = (/31,28,31,30,31,30,31,31,30,31,30,31/)

year_po = (/1980, 1983, 1993, 1998, 1999, 2003, 2014, 2015/)
year_ne = (/1981, 1984, 1985, 1990, 1994, 1997, 2001, 2013, 2018/)
year_PJ = array_append_record(year_po,year_ne,0)
qsort(year_PJ)
;===================================================================================
;data_list
;===================================================================================
  diri = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/Regrid_Monthly/RAINNC/"
  list_in = systemfunc("ls "+diri+"*{1980..2019}*.nc")
  f = addfiles(list_in,"r")
  ListSetType(f,"join")

  data = f[:]->RAINNC
  data = data*24   ;mm/hour > mm/day

  lat = data&Lat
  lon = data&Lon

  nlat = dimsizes(lat)
  nlon = dimsizes(lon)
;===================================================================================
;mask mountain site
;===================================================================================
  H = 1000
  height = get_Taiwan_height(fspv,lat,lon)
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  scalarH = conform(data,height,(/2,3/))
  data = data*scalarH

  ;multiply days to monthly value
  scalar_M = conform(data,month_day,1)
  data = data * scalar_M
  delete([/H,height,scalarH,scalar_M/])  

  data_re = reshape(data,(/nyear*12,nlat,nlon/))
;===================================================================================
;address data
;===================================================================================
;1.calculate SPI
  run = 3
  opt = False

  spi = dim_spi_n(data_re,run,opt,0)
  delete(data_re)

  spi_year = reshape(spi,dimsizes(data))
  copy_VarCoords(data,spi_year)
  spi_year!0 = "year"
  spi_year&year = year
  spi_JJA  = spi_year(:,5:7,:,:)

  delete([/spi,spi_year/])
;===================================================================================
;calculate PDF
;===================================================================================
;1.separate PJ
  spi_po = spi_JJA({year_po},:,:,:)
  spi_ne = spi_JJA({year_ne},:,:,:)
  spi_nu = spi_JJA
  spi_nu({year_PJ},:,:,:) = fspv
  delete(spi_JJA)

;2. calculate the events rate
  P_all = num(.not.ismissing(spi_po))
  N_all = num(.not.ismissing(spi_ne))
  A_all = num(.not.ismissing(spi_nu))

print("---Dry---")
  PD0 = num(spi_po.lt.0.0)
  NW0 = num(spi_ne.gt.0.0)

  ;print(100*tofloat(PD0)/tofloat(P_all))
  ;print(100*tofloat(NW0)/tofloat(N_all))
  
  print("")
  print("<-1.0")
  P_1 = num(spi_po.lt.-1.0.and.spi_po.ge.-1.5)
  N_1 = num(spi_ne.lt.-1.0.and.spi_ne.ge.-1.5)
  A_1 = num(spi_nu.lt.-1.0.and.spi_nu.ge.-1.5)

  print(""+P_1)
  print(""+N_1)
  print(""+A_1)

  ;print(100*tofloat(P_1)/tofloat(P_all))
  ;print(100*tofloat(N_1)/tofloat(N_all))
  ;print(100*tofloat(A_1)/tofloat(A_all))

  print("")
  print("<1.5")
  P_15 = num(spi_po.lt.-1.5.and.spi_po.ge.-2.0)
  N_15 = num(spi_ne.lt.-1.5.and.spi_ne.ge.-2.0)
  A_15 = num(spi_nu.lt.-1.5.and.spi_nu.ge.-2.0)

  print(""+P_15)
  print(""+N_15)
  print(""+A_15)

  ;print(100*tofloat(P_15)/tofloat(P_all))
  ;print(100*tofloat(N_15)/tofloat(N_all))
  ;print(100*tofloat(A_15)/tofloat(A_all))

  print("")
  print("<2.0")
  P_2  = num(spi_po.lt.-2.0)
  N_2  = num(spi_ne.lt.-2.0)
  A_2  = num(spi_nu.lt.-2.0)

  print(""+P_2)
  print(""+N_2)
  print(""+A_2)

  ;print(100*tofloat(P_2)/tofloat(P_all))
  ;print(100*tofloat(N_2)/tofloat(N_all))
  ;print(100*tofloat(A_2)/tofloat(A_all))

  print("")
  print("all dry")
  ;print(100*tofloat(num(spi_po.lt.-1.0))/tofloat(P_all))
  ;print(100*tofloat(num(spi_ne.lt.-1.0))/tofloat(N_all))
  ;print(100*tofloat(num(spi_nu.lt.-1.0))/tofloat(A_all))
  print(num(spi_po.lt.-1.0))
  print(num(spi_ne.lt.-1.0))
  print(num(spi_nu.lt.-1.0))
  exit

;3. calculate each probability
print("---calculate PDF---")

  opt = True
  opt@bin_min = -5.0
  opt@bin_max = 5.0
  bins = 40

  pdf_po = pdfx(spi_po,bins,opt)
  pdf_ne = pdfx(spi_ne,bins,opt)
  pdf_nu = pdfx(spi_nu,bins,opt)

  delete([/spi_po,spi_ne,spi_nu/])
;===================================================================================
;calculate PDF
;===================================================================================
print("---cal PDF---")
  nVar = 3
  nBin = bins          ; retrieve the number of bins

  xx = new((/nVar, nBin/), typeof(pdf_nu))

  xx(0,:) = pdf_po@bin_center
  xx(1,:) = pdf_ne@bin_center
  xx(2,:) = pdf_nu@bin_center

  yy = new((/nVar, nBin/), typeof(pdf_nu))
  yy(0,:) = (/pdf_po/)
  yy(1,:) = (/pdf_ne/)
  yy(2,:) = (/pdf_nu/)
  
  delete([/pdf_nu,pdf_po,pdf_ne/])
;===================================================================================
;plot
;===================================================================================
print("---plot---")
  fig_type = "X11"
  fig_type = "eps"

  fig_name = "PDF_TReAD_SPI3"

  wks = gsn_open_wks(fig_type,fig_name)

  res  = True
  res@vpWidthF  = 0.7
  res@vpHeightF = 0.3

  res@gsnDraw  = False
  res@gsnFrame = False
  
  res@xyLineThicknessF         = 3.0
  res@xyLineColors             = (/"red","blue","black"/)  
  res@xyMonoDashPattern        = True              ; all solid 

  ;--label----------------------------------------------------------------------
  label = ispan(nVar-1,0,1)
  res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Right"             ; Change location of 
  res@pmLegendParallelPosF   = 0.87                ; move units vertical
  res@pmLegendOrthogonalPosF = -0.3              ; move units horizontal
  res@pmLegendWidthF         = 0.08                ; Change width and
  res@pmLegendHeightF        = 0.08                ; height of legend.

  res@lgPerimOn              = False               ; turn off/on box around
  res@lgLabelFont            = 22                  ; label font  
  res@lgLabelFontHeightF     = 0.012               ; label font height
  res@lgItemOrder            = label
  ;res@xyExplicitLegendLabels = " "+(/"1960-2019","Positive","Negative"/)  ; create explicit labels
  res@xyExplicitLegendLabels = " "+(/"Positive Year","Negative Year","Neutral"/)  ; create explicit labels
  ;------------------------------------------------------------------------------
  
  ;x-y
  res@tmXTOn = False
  ;res@tmYROn = False

  ;tmax
  res@trXMinF = -5
  res@trXMaxF = 5
  ;res@trYMinF = 0
  ;res@trYMaxF = 10

  ;res@tmXBMode = "Explicit"
  ;res@tmXBValues = ispan(25,40,5)
  ;res@tmXBLabels = "" + res@tmXBValues
  ;res@tmXBMinorValues = "" + ispan(25,40,1)

  res@tiXAxisFont = 21
  res@tiXAxisFontHeightF = 0.016
  res@tiYAxisFont = 21
  res@tiYAxisFontHeightF = 0.016

  res@tiXAxisString = "SPI"
  res@tiYAxisString = "PDF (%)"

  ;string
  ;res@tiMainFont = 22               
  ;res@tiMainFontHeightF = 0.022               

  res@gsnStringFont = 21           
  res@gsnRightStringFontHeightF = 0.015           
  res@gsnRightStringOrthogonalPosF = 0.015
  ;res@gsnRightStringParallelPosF = 0.01
  res@gsnRightString = ""

  res@gsnLeftStringFontHeightF  = 0.015         
  res@gsnLeftString  = ""
  res@gsnLeftStringOrthogonalPosF = 0.02
 
  ;res@tiMainString  = varname + " PDF"
  plot = gsn_csm_xy (wks,xx,yy,res)

  respoly = True
  respoly@gsLineColor = "gray30"
  respoly@gsLineThicknessF = 1.0
  respoly@gsLineDashPattern = 1
  dum0 = gsn_add_polyline(wks,plot,(/0,0/),(/0,12/),respoly)
  ;dum1 = gsn_add_polyline(wks,plot,(/35,35/),(/0,10/),respoly)

  ;rest = True
  ;rest@txFont = 21
  ;rest@txFontColor = "gray50"
  ;rest@txFontHeightF = 0.01
  ;dum3 = gsn_add_text(wks,plot,"33~S~o~N~C",32.5,9.6,rest)
  ;dum4 = gsn_add_text(wks,plot,"35~S~o~N~C",35.5,9.6,rest)
  ;respoly@gsFillColor  = "gray70"
  ;respoly@gsFillOpacityF = 0.2
  ;dum = gsn_add_polygon(wks,plot,(/33,33,35,35,33/),(/0,6,6,0,0/),respoly)

  draw(plot)
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
