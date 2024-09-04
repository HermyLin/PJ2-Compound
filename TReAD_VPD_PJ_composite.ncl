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
;-----------------------------------------------
;give date
;-----------------------------------------------
;date_split = str_split(systemfunc("date '+%D'"),"/")
;date = str_concat((/date_split(2),date_split(0),date_split(1),"_"/))
;===================================================================================
;read
;===================================================================================

  fn = "/lfs/home/killumy/02-PJ_LAI/0_data/TReAD_2km_JJAdaily_VPD.nc"

  f = addfile(fn,"r")
  data = dim_avg_n_Wrap(f->VPD({year},:,:,:),1)
  data = data / 100 

  lat  = f->Lat
  lon  = f->Lon
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)


  delete([/fn,f/])
;===================================================================================
;mask mountain site
;===================================================================================
  H = 1000
  height = get_Taiwan_height(fspv,lat,lon)
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  do t = 0, nyear-1
     data(t,:,:) = data(t,:,:) * height
  end do
  delete(height)

  data = nameDim(data,(/"time","lat","lon"/),"","hPa")
;===================================================================================
;address data
;===================================================================================
;;1.detrend
  data_dtr = dtrend_msg_n(data&time,data,True,False,0)
  copy_VarCoords(data,data_dtr)
  data := data_dtr   ;40yrs
  delete(data_dtr)
  
;2.anomaly on grids 
  data_clm = dim_avg_n_Wrap(data,0)  
  data_ano = data

  do yy = 0, nyear-1
     data_ano(yy,:,:) = data(yy,:,:) - data_clm
  end do

  data_po = data({year_po},:,:)  ;positive year
  data_ne = data({year_ne},:,:)  ;negative year
  data_ano_po = dim_avg_n_Wrap(data_ano({year_po},:,:),0)
  data_ano_ne = dim_avg_n_Wrap(data_ano({year_ne},:,:),0)
;===================================================================================
;t-test correlation
;ttest (ave1,var1,dims1,ave2,var2,dims2,iflag, tval_opt)
;===================================================================================
  prob_po = SigTest(fspv,lat,data,data_po)
  prob_po@_FillValue = fspv
  ;ind_po  = ind_resolve(ind(.not.ismissing(prob_po)),(/nlat,nlon/))

  prob_ne = SigTest(fspv,lat,data,data_ne)
  prob_ne@_FillValue = fspv
  ;ind_ne  = ind_resolve(ind(.not.ismissing(prob_ne)),(/nlat,nlon/))

  ;data_ratio = data_ano_po-data_ano_ne
  ;copy_VarCoords(data_ano_po,data_ratio)
  copy_VarCoords(data_ano_po,prob_po)
  copy_VarCoords(data_ano_ne,prob_ne)

  delete([/data_po,data_ne,data/])
;===================================================================================
;plot
;===================================================================================
print("---plot---")
  fig_type = "X11"
  fig_type = "eps"
  fig_name = "VPD_map"
  plot = new(2,graphic)

  wks = gsn_open_wks(fig_type,fig_name)

  res  = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@tmBorderThicknessF = 1.0
  
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
  res@lbLabelBarOn = False
  cmap = read_colormap_file("MPL_BrBG")
  res@cnFillPalette = cmap(::-1,:)
  res@cnMissingValFillColor = "white"

  res@cnLevelSelectionMode = "ManualLevels" ;color bar adjust mode
  res@cnMaxLevelValF = 1.5
  res@cnMinLevelValF = -1.5
  res@cnLevelSpacingF = 0.3 ;color bar spacing
  
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
  ;res@tmYROn = False

  ;map
  res@mpFillOn = True
  res@cnMissingValFillColor = "transparent" 
  res@mpOceanFillColor      = "white"
  res@mpLandFillColor       = "gray"

  res@mpDataBaseVersion = "HighRes"
  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpOutlineDrawOrder = "PostDraw"

  res@cnFillMode = "RasterFill"
  res@gsnAddCyclic = False
  
  plot(0) = gsn_csm_contour_map(wks,data_ano_po,res)
  plot(1) = gsn_csm_contour_map(wks,data_ano_ne,res)

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
  resL@cnFillScaleF = 0.8  ;增加形状填充的密度
  ;resL@cnFillDotSizeF  = 0.003
  resL@cnFillColor = "black"
  resL@cnFillPattern = 6
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
  
  delete([/f,lat,lon/])
;-----map---------------------------------------------------------------------------
  overlay(plot(0),map0)
  overlay(plot(1),map1)
  
;------panel-------------------------------------------------------------------------
  resP = True
  ;resP@gsnFrame = False
  resP@gsnMaximize = True
  resP@gsnPaperOrientation = "Portrait"

  resP@gsnPanelLabelBar = True
  resP@lbOrientation = "Vertical"
  resP@lbBoxEndCapStyle = "TriangleBothEnds"

  resP@lbLabelFont         =  22        
  resP@lbLabelFontHeightF  =  0.016        ; color bar的字大小
  resP@pmLabelBarOrthogonalPosF = 0.01     ; color bar垂直擺放位置
  ;resP@pmLabelBarParallelPosF = 0.01

  resP@pmLabelBarHeightF = 0.50
  resP@pmLabelBarWidthF = 0.08

  gsn_panel(wks,plot,(/1,2/),resP)
;---change figure type----------------------------------------------------------

  delete(wks)
  if (fig_type.eq."eps") then
  system("convert -trim -density 300 -crop 0x0 "+fig_name+".eps  "+fig_name+".png")
         ;==convert + quality(resolution) up + remove white boundary==;

  system("/bin/rm "+fig_name+".eps")
  end if
  print("Output Figure Name = "+fig_name+".png")

end
