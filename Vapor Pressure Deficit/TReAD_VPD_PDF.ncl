load "/lfs/home/killumy/functions/func_PJ.ncl"
;VPD daily JJA PDF


begin

fspv = 1e+36
ispv = -999
year = ispan(1980,2019,1)
nyear = dimsizes(year)
year_po = (/1980, 1983, 1993, 1998, 1999, 2003, 2014, 2015/)
year_ne = (/1981, 1984, 1985, 1990, 1994, 1997, 2001, 2013, 2018/)
year_PJ = array_append_record(year_po,year_ne,0)
qsort(year_PJ)
;===================================================================================
;read: daily data need
;===================================================================================
  fn = "/lfs/home/killumy/02-PJ_LAI/00-data/TReAD_2km_JJAdaily_VPD.nc"

  f = addfile(fn,"r")
  data = f->VPD({year},:,:,:)
  data = data / 1000

  ndim = dimsizes(data)
  nday = ndim(1)
  nlat = ndim(2)
  nlon = ndim(3)

  lat = f->Lat
  lon = f->Lon

  data_3d = reshape(data,(/nyear*nday,nlat,nlon/))
  ntim = nyear*nday

  delete([/fn,f/])
;===================================================================================
;mask mountain site
;===================================================================================
  H = 1000
  height = get_Taiwan_height(fspv,lat,lon)
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  do t = 0, ntim-1
     data_3d(t,:,:) = data_3d(t,:,:) * height
  end do

  data_H = reshape(data_3d,(/nyear,nday,nlat,nlon/))
  copy_VarCoords(data,data_H)

  delete([/height,ntim/])
  delete([/data_3d,data/])
;===================================================================================
;calculate PDF
;===================================================================================
  mean = avg(data_H)
;1.separate PJ
  VPD_po = data_H({year_po},:,:,:)
  VPD_ne = data_H({year_ne},:,:,:)
  VPD_nu = data_H
  VPD_nu({year_PJ},:,:,:) = fspv

;cal percentage
  DRY_P = tofloat(num((VPD_po.gt.1.0)))
  DRY_N = tofloat(num((VPD_ne.gt.1.0)))

  ALL = tofloat(num(.not.ismissing(data_H)))
  print(DRY_P/ALL)
  print(DRY_N/ALL)
  exit

  delete(data_H)



;2. calculate each probability
print("---calculate PDF---")

  opt = True
  opt@bin_min = 0.
  opt@bin_max = 3.5
  bins = 35

  pdf_po = pdfx(VPD_po,bins,opt)
  pdf_ne = pdfx(VPD_ne,bins,opt)
  pdf_nu = pdfx(VPD_nu,bins,opt)

  delete([/VPD_po,VPD_ne,VPD_nu/])
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

  fig_name = "VPD_PDF_daily_V2"

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
  res@pmLegendParallelPosF   = 0.80                ; move units vertical
  res@pmLegendOrthogonalPosF = -0.38              ; move units horizontal
  res@pmLegendWidthF         = 0.08                ; Change width and
  res@pmLegendHeightF        = 0.12                ; height of legend.

  res@lgPerimOn              = False               ; turn off/on box around
  res@lgLabelFont            = 22                  ; label font  
  res@lgLabelFontHeightF     = 0.02               ; label font height
  res@lgItemOrder            = label
  ;res@xyExplicitLegendLabels = " "+(/"1960-2019","Positive","Negative"/)  ; create explicit labels
  res@xyExplicitLegendLabels = " "+(/"Positive Year","Negative Year","Neutral Year"/)  ; create explicit labels
  ;------------------------------------------------------------------------------
  
  ;x-y
  res@tmXTOn = False
  ;res@tmYROn = False

  ;tmax
  res@trXMinF = 0
  res@trXMaxF = 2
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

  res@tiXAxisString = "VPD (kPa)"
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
