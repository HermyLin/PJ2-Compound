load "/lfs/home/killumy/functions//func_PJ.ncl"
;===================================================================================
;------------------------main code------------------------
;===================================================================================
begin

;begTime = get_cpu_time()
fspv = 1e+36
ispv = -999
year = ispan(1980,2019,1)
nyear = dimsizes(year)
month = ispan(1,12,1)
nmonth = dimsizes(month)

month_day = (/31,28,31,30,31,30,31,31,30,31,30,31/)

year_po = (/1980, 1983, 1993, 1998, 1999, 2003, 2014, 2015/)
year_ne = (/1981, 1984, 1985, 1990, 1994, 1997, 2001, 2013, 2018/)
;===================================================================================
;data_list
;===================================================================================
  diri = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/Regrid_Monthly/"
  fn   = "TReAD_2km_Monthly_PRECTmms.nc"
  f = addfile(diri+fn,"r")

  data = f->PRECTmms({year},:,:,:)
  data@_FillValue = fspv
  data@missing_value = fspv

  data = data*86400    ;mm/day

  lat = f->Lat
  lon = f->Lon
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)

  data = nameDim(data,(/"time","month","lat","lon"/),"","")

  ;multiply to monthly value
  scalar = conform(data,month_day,1)
  prep = data*scalar ;mm/mon
  copy_VarCoords(data,prep)

  delete([/data,scalar/])
;===================================================================================
;mask mountain site
;===================================================================================
  H = 1000
  height = get_Taiwan_height(fspv,lat,lon)
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  scalarH = conform(prep,height,(/2,3/))
  prep = prep*scalarH

  prep_time = reshape(prep,(/nyear*nmonth,nlat,nlon/))

  delete([/height,scalarH/])
;===================================================================================
;address data
;===================================================================================
;1.calculate SPI
  run = 3
  opt = False

  spi = dim_spi_n(prep_time,run,opt,0)
  spi_year = reshape(spi,(/nyear,nmonth,nlat,nlon/))
  spi_JJA  = dim_avg_n(spi_year(:,5:7,:,:),1)

  copy_VarCoords(prep(:,0,:,:),spi_JJA)

  spi_JJA!0 = "time"
  spi_JJA&time = year
  delete([/spi,spi_year/])

  spi_po = spi_JJA({year_po},:,:)
  spi_ne = spi_JJA({year_ne},:,:)
;===================================================================================
;t-test correlation
;ttest (ave1,var1,dims1,ave2,var2,dims2,iflag, tval_opt)
;===================================================================================
  prob_po = SigTest(fspv,lat,spi_JJA,spi_po)
  prob_po@_FillValue = fspv
  
  prob_ne = SigTest(fspv,lat,spi_JJA,spi_ne)
  prob_ne@_FillValue = fspv
  
  data_po = dim_avg_n_Wrap(spi_po,0)
  data_ne = dim_avg_n_Wrap(spi_ne,0)

  copy_VarCoords(data_po,prob_po)
  copy_VarCoords(data_ne,prob_ne)

  delete([/spi_JJA,spi_po,spi_ne/])
;===================================================================================
;plot
;===================================================================================
print("---plot---")

  fig_type = "X11"
  fig_type = "eps"
  fig_name = "./PNG_Composite/Composite_TReAD_SPI"
  plot = new(2,graphic)

  wks = gsn_open_wks(fig_type,fig_name)

  res  = True
  res@gsnDraw  = False
  res@gsnFrame = False
  
  res@cnFillOn  = True                ; turn on color
  res@cnLinesOn = False               ; no contour lines, lines off

  res@cnLineLabelsOn = False
  res@mpMinLatF = 21.5                          ;set lat interval
  res@mpMaxLatF = 25.5
  res@mpMinLonF = 119.2                         ;set lon interval
  res@mpMaxLonF = 122.3
  ;res@mpCenterLonF =  180.
  res@gsnMajorLatSpacing = 1
  res@gsnMajorLonSpacing = 1

  ;labelbar
  cmap = read_colormap_file("MPL_BrBG")
  res@cnFillPalette = cmap

  res@lbLabelBarOn = False

  res@cnLevelSelectionMode = "ManualLevels" ;color bar adjust mode
  res@cnLevelSpacingF = 0.2 ;color bar spacing
  res@cnMinLevelValF = -1.0
  res@cnMaxLevelValF = 1.0
  
  ;label
  res@tmXBLabelFont = 21         
  res@tmXBLabelFontHeightF = 0.017         
  res@tmXBMajorLengthF = 0.01
  res@tmXBMajorOutwardLengthF = 0.01
  res@tmXTOn = False

  res@tmYLLabelFont = 21       
  res@tmYLLabelFontHeightF = 0.017       
  res@tmYLMajorLengthF = 0.01
  res@tmYLMajorOutwardLengthF = 0.01

  ;map
  res@mpFillOn = True
  res@cnMissingValFillColor = "transparent" 
  res@mpOceanFillColor      = "white"
  res@mpLandFillColor       = "transparent"

  res@mpDataBaseVersion = "HighRes"
  res@mpOutlineBoundarySets = "AllBoundaries"
  ;res@mpFillDrawOrder    = "PostDraw"
  res@mpOutlineDrawOrder = "PostDraw"

  res@cnFillMode = "RasterFill"
  res@gsnAddCyclic = False
  plot(0) = gsn_csm_contour_map(wks,data_po,res)

  ;res@tmYROn = False
  plot(1) = gsn_csm_contour_map(wks,data_ne,res)
  delete(res)

;-----significant-------------------------------------------------------------
  resL  = True
  resL@gsnDraw  = False
  resL@gsnFrame = False

  resL@cnFillOn =True
  resL@cnLinesOn = False               ; no contour lines, lines off
  resL@cnLineLabelsOn = False
  resL@cnInfoLabelOn = False

  resL@lbLabelBarOn = False
  resL@cnMaxLevelValF = 0.1
  resL@cnMinLevelValF = 0.0

  resL@cnMonoFillPattern = True               ; no multiple patterns
  resL@cnMonoFillColor = True
  resL@cnFillScaleF = 0.7  
  resL@cnFillColor = "black"
  resL@cnLineThicknessF  = 2.0
  resL@cnFillPattern = 6
  ;resL@cnFillDotSizeF  = 0.004
  ;resL@cnFillPattern = 17
  resL@cnConstFEnableFill = False
  
  dum0 = gsn_csm_contour(wks,prob_po,resL)
  dum1 = gsn_csm_contour(wks,prob_ne,resL)
  delete(resL)
  ;-----significant-------------------------------------------------------------
  overlay(plot(0),dum0)
  overlay(plot(1),dum1)
  
  ;-----map---------------------------------------------------------------------------
  resH = True
  resH@gsnDraw  = False
  resH@gsnFrame = False
  
  resH@cnFillOn  = False                ; turn on color
  resH@cnLinesOn = True               ; no contour lines, lines off
  resH@cnLineLabelsOn = False
  resH@cnInfoLabelOn = False

  resH@cnLineThicknessF = 2
  resH@cnLevelSelectionMode = "ExplicitLevels" ;color bar adjust mode
  resH@cnLevels = (/1000,2000/)
  ;resH@cnLineDashPatterns = (/1,1,0/)
  resH@cnMonoLineColor = False
  resH@cnLineColors = (/"black","gray50"/)

  resH@tiMainString = ""
  resH@gsnLeftString = ""
  resH@gsnRightString = ""

  f = addfile("/lfs/home/killumy/01-PJ_OBS/0data/CTX_variables_mask1000.nc","r")
  lat5km = f->lat
  lon5km = f->lon
   
  height = get_Taiwan_height(fspv,lat5km,lon5km)
  map0 = gsn_csm_contour(wks,height,resH)
  map1 = gsn_csm_contour(wks,height,resH)
  
  delete(resH)
  delete(f)
;-----map---------------------------------------------------------------------------
  overlay(plot(0),map0)
  overlay(plot(1),map1)
  
;panel
  resP                            = True
  ;resP@gsnPanelXWhiteSpacePercent = 10.
  ;resP@gsnPanelYWhiteSpacePercent = 10.

  resP@gsnPanelLabelBar = True
  resP@lbOrientation = "Vertical"
  resP@lbBoxEndCapStyle = "TriangleBothEnds"
  resP@cnMissingValFillColor = "white"

  resP@lbLabelFont = 22
  resP@lbLabelFontHeightF = 0.015
  resP@pmLabelBarHeightF = 0.5
  resP@pmLabelBarWidthF = 0.08
  resP@pmLabelBarOrthogonalPosF = 0.01
  ;resP@pmLabelBarParallelPosF = 0.5
  gsn_panel(wks,plot,(/1,2/),resP)                      ; now draw as one plot
  delete(resP)
;---change figure type----------------------------------------------------------

  delete(wks)
  if (fig_type.eq."eps") then
  system("convert -trim -density 300 -crop 0x0 "+fig_name+".eps  "+fig_name+".png")
         ;==convert + quality(resolution) up + remove white boundary==;

  system("/bin/rm "+fig_name+".eps")
  end if
  print("Output Figure Name = "+fig_name+".png")

end
