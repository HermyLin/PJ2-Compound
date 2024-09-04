load "/lfs/home/killumy/functions/func_PJ.ncl"

;===================================================================================
;------------------------main code------------------------
;===================================================================================
begin

;begTime = get_cpu_time()
fspv = 1e+20
ispv = -999
year = ispan(1980,2019,1)
nyear = dimsizes(year)

year_po = (/1980, 1983, 1993, 1998, 1999, 2003, 2014, 2015/)
year_ne = (/1981, 1984, 1985, 1990, 1994, 1997, 2001, 2013, 2018/)

year_PJ = array_append_record(year_po,year_ne,0)
qsort(year_PJ)
;===================================================================================
;data_list
;===================================================================================
print("---read---")

  diri = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/Daily/" 
  fn = "TReAD_2km_Daily_SMOIS.nc"

  f = addfile(diri+fn,"r")

  lat  = f->Lat
  lon  = f->Lon
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)

  soil1 = (/f->Soil1({year},181:272,:,:)/)   ;layer, year, month, lat, lon
  soil1@_FillValue = fspv
  soil1@missing_value = fspv

  soil2 = (/f->Soil2({year},181:272,:,:)/)   ;layer, year, month, lat, lon
  soil2@_FillValue = fspv
  soil2@missing_value = fspv

  data = ((soil1*0.07) + (soil2*0.21)) / 0.28
  delete([/soil1,soil2/])
;===================================================================================
;mask mountain site
;===================================================================================
print("---height---")
  H = 1000
  height = get_Taiwan_height(fspv,lat,lon)
  height@_FillValue = fspv
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  scalarH = conform(data,height,(/2,3/))

  data = data*scalarH
  delete([/height,scalarH/])
  delete([/lat,lon,nlat,nlon/])

  data!0 = "year"
  data&year = year
;===================================================================================
;address data
;===================================================================================
  mean = avg(data)
  print(""+mean)
b = ndtooned(data({year_po},:,:,:))
bnum = num(.not.ismissing(b))

c = ndtooned(data({year_ne},:,:,:))
cnum = num(.not.ismissing(c))

b1 = num(b.lt.mean)
c2 = num(c.lt.mean)
print(tofloat(b1)/tofloat(bnum))
print(tofloat(c2)/tofloat(cnum))

;1.separate PJ
  data_po = data({year_po},:,:,:)
  data_ne = data({year_ne},:,:,:)
  data_nu = data
  data_nu({year_PJ},:,:,:) = fspv

;2. calculate each probability
  ;setting bins
  ;opt         = True
  opt = True
  opt@bin_min = 0.0
  opt@bin_max = 0.8
  bins = 40

  pdf_al = pdfx(data,bins,opt)
  pdf_po = pdfx(data_po,bins,opt)
  pdf_ne = pdfx(data_ne,bins,opt)
  pdf_nu = pdfx(data_nu,bins,opt)
;===================================================================================
;percentage of Temp > value
;===================================================================================
  ;Temp_std = 35

  ;a = tofloat(num(temp.gt.Temp_std))/tofloat(num(.not.ismissing(temp)))
  ;b = tofloat(num(data_po.gt.Temp_std))/tofloat(num(.not.ismissing(data_po)))
  ;c = tofloat(num(data_ne.gt.Temp_std))/tofloat(num(.not.ismissing(data_ne)))
  ;d = tofloat(num(data_nu.gt.Temp_std))/tofloat(num(.not.ismissing(data_nu)))

  ;print("60yr = "+100*a + "%")
  ;print("Threshold: "+Temp_std)
  ;print("POyr = "+100*b + "%")
  ;print("NEyr = "+100*c + "%")
  ;print("NUyr = "+100*d + "%")
;===================================================================================
;calculate PDF
;===================================================================================
print("---cal PDF---")
  nVar = 3
  nBin = pdf_al@nbins          ; retrieve the number of bins

  xx = new((/nVar, nBin/), typeof(pdf_al))

  xx(0,:) = pdf_po@bin_center
  xx(1,:) = pdf_ne@bin_center
  xx(2,:) = pdf_nu@bin_center

  yy = new ((/nVar, nBin/), typeof(pdf_al))
  yy(0,:) = (/pdf_po/)
  yy(1,:) = (/pdf_ne/)
  yy(2,:) = (/pdf_nu/)
  
  delete([/data_po,data_ne/])
  delete([/pdf_al,pdf_po,pdf_ne,pdf_nu/])
;===================================================================================
;plot
;===================================================================================
print("---plot---")
  fig_type = "X11"
  fig_type = "eps"

  ;out_dir  = "/lfs/home/killumy/PJ_pattern/paper/"
  fig_name = "PDF_SMOIS_JASdaily_TReAD"

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
  res@pmLegendSide           = "Left"             ; Change location of 
  res@pmLegendParallelPosF   = 0.80                ; move units vertical
  res@pmLegendOrthogonalPosF = -0.47              ; move units horizontal
  res@pmLegendWidthF         = 0.07                ; Change width and
  res@pmLegendHeightF        = 0.12                ; height of legend.

  res@lgPerimOn              = False               ; turn off/on box around
  res@lgLabelFont            = 22                  ; label font  
  res@lgLabelFontHeightF     = 0.020               ; label font height
  res@lgItemOrder            = label
  ;res@xyExplicitLegendLabels = " "+(/"1960-2019","Positive","Negative"/)  ; create explicit labels
  res@xyExplicitLegendLabels = " "+(/"Positive Year","Negative Year","Neutral Year"/)  ; create explicit labels
  ;------------------------------------------------------------------------------
  
  ;x-y
  res@tmXTOn = False
  ;res@tmYROn = False

  ;tmax
  res@trXMinF = 0.0
  res@trXMaxF = 0.6

  ;res@tmXBMode = "Explicit"
  ;res@tmXBValues = ispan(10,40,5)
  ;res@tmXBLabels = "" + res@tmXBValues

  res@tiXAxisFont = 21
  res@tiXAxisFontHeightF = 0.016
  res@tiYAxisFont = 21
  res@tiYAxisFontHeightF = 0.016

  res@tiXAxisString = "SMOIS (mm~S~3~N~/mm~S~3~N~)"
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
  dum0 = gsn_add_polyline(wks,plot,(/mean,mean/),(/0,20/),respoly)

  ;rest = True
  ;rest@txFont = 21
  ;rest@txFontColor = "gray50"
  ;rest@txFontHeightF = 0.01
  ;dum3 = gsn_add_text(wks,plot,"33~S~o~N~C",32,4.8,rest)
  ;dum4 = gsn_add_text(wks,plot,"35~S~o~N~C",36,4.8,rest)
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
