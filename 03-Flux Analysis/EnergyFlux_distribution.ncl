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
JJA_DOY = ispan(5,7,1) ;JJAS: 5-8
;===================================================================================
;read
;vars: SSR,STR,SHF,LHF
;===================================================================================
;1. read
  diri = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/Regrid_Monthly/"

  ;SSR
  f1 = addfile(diri+"TReAD_2km_Monthly_FSDS.nc","r")
  f2 = addfile(diri+"TReAD_2km_Monthly_SWUP.nc","r")

  SSR = (f1->FSDS) - (f2->SWUP)
  SSR := dim_avg_n(SSR(:,JJA_DOY,:,:),1)
  delete([/f1,f2/])
  
  ;STR
  f1 = addfile(diri+"TReAD_2km_Monthly_FLDS.nc","r")
  f2 = addfile(diri+"TReAD_2km_Monthly_LWUP.nc","r")

  STR = (f1->FLDS) - (f2->LWUP)
  STR := dim_avg_n(STR(:,JJA_DOY,:,:),1)
  delete([/f1,f2/])

  ;SHF
  f = addfile(diri+"TReAD_2km_Monthly_SHF.nc","r")
  SHF = dim_avg_n_Wrap(f->SHF(:,JJA_DOY,:,:),1)
  delete(f)

  ;LHF
  f = addfile(diri+"TReAD_2km_Monthly_LHF.nc","r")
  LHF = dim_avg_n_Wrap(f->LHF(:,JJA_DOY,:,:),1)

  ;coordinate
  lat = LHF&Lat
  lon = LHF&Lon
;===================================================================================
;organize data
;===================================================================================
  data = (/SSR, STR, -1*SHF, -1*LHF/)   ;make downward positive
  data@_FillValue = fspv

  copy_VarCoords(LHF,data(0,:,:,:))
  delete([/SSR,STR,SHF,LHF/])

  data_JJA = data(:,{year},:,:)
  data_JJA!1 = "time"
  data_JJA!2 = "lat"
  data_JJA!3 = "lon"
  delete(data)
;===================================================================================
;mask height
;===================================================================================
  H = 1000
  height = get_Taiwan_height(fspv,lat,lon)
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  do n = 0, 3
  do t = 0, nyear-1
     data_JJA(n,t,:,:) = data_JJA(n,t,:,:) * height
  end do
  end do
  delete(height)
;===================================================================================
;address data
;===================================================================================
;;1.detrend
  data_dtr = dtrend_msg_n(data_JJA&time,data_JJA,True,False,1)
  copy_VarCoords(data_JJA,data_dtr)
  data := data_dtr   
  delete([/data_JJA,data_dtr/])
  
;2.anomaly on grids 
  ;coords -> var, year, lat, lon
  data_clm = dim_avg_n_Wrap(data,1)  
  data_ano = data

  do yy = 0, nyear-1
     data_ano(:,yy,:,:) = data(:,yy,:,:) - data_clm
  end do

  data_po = data(:,{year_po},:,:)  ;positive year
  data_ne = data(:,{year_ne},:,:)  ;negative year

  data_ano_po = dim_avg_n_Wrap(data_ano(:,{year_po},:,:),(/1,2,3/))
  data_ano_ne = dim_avg_n_Wrap(data_ano(:,{year_ne},:,:),(/1,2,3/))
  data_po := array_append_record(data_ano_po,sum(data_ano_po),0)
  data_ne := array_append_record(data_ano_ne,sum(data_ano_ne),0)
  delete([/data_ano_po,data_ano_ne/])

  print(data_po)
  print(data_ne)
;===================================================================================
;plot
;bar plot
;===================================================================================
print("---plot---")

  fig_type = "X11"
  fig_type = "eps"
  fig_name = "TReAD_JJA_Flux_distribution"

  plot = new(2,graphic)

  wks = gsn_open_wks (fig_type,fig_name)               ; send graphics to PNG file
  res = True                     ; plot mods desired

  res@gsnDraw = False
  res@gsnFrame = False
  
  res@vpHeightF = 0.4
  res@vpWidthF = 0.7
  ;res@gsnMaximize = True
  ;res@gsnPaperOrientation = "Portrait"

  res@gsnXYBarChart  = True                 ; turn on bar chart
  res@gsnXYBarChartBarWidth = 0.6
  res@gsnYRefLine    = 0.
  res@gsnAboveYRefLineColor = "black"
  res@gsnBelowYRefLineColor = "black"
  
  res@tmXBLabelFont = 21
  res@tmXBLabelFontHeightF = 0.02

  res@tmXBMode        = "Explicit"	
  res@tmXBValues      = (/0,1,2,3,4/) 
  res@tmXBLabels      = (/"SSR","STR","SHF","LHF","sum"/)
  ;res@tmXBMinorValues = ispan(1960,2019,2)

  res@tmYLMode = "Manual"
  res@tmYLTickSpacingF = 4.0
  res@tmYLTickStartF   = -8.0
  res@tmYLTickEndF     = 8.0
  ;res@tmYLMinorOn = True

  res@tmYLLabelFont = 21
  res@tmYLLabelFontHeightF = 0.02 
  res@tmYLMinorOn = False

  res@trXMinF = -0.5
  res@trXMaxF = 4.5
  res@trYMinF = -8.0
  res@trYMaxF = 8.0
  
  x = (/0,1,2,3,4/)
  plot(0) = gsn_csm_xy(wks,x,data_po,res) ; create plot
  plot(1) = gsn_csm_xy(wks,x,data_ne,res) ; create plot

;panel
  resP = True
  ;resP@gsnMaximize = True
  ;resP@gsnPaperOrientation = "Portrait"
  ;resP@gsnPanelMainString         = "EOF mode1 of 1950-2021 JJA Precipitation"
  ;resP@gsnPanelMainFontHeightF    = 0.02
  ;resP@gsnPanelXWhiteSpacePercent = 4.
  ;resP@gsnPanelYWhiteSpacePercent = 3.
  gsn_panel(wks,plot,(/2,1/),resP) 
;---change figure type----------------------------------------------------------

  delete(wks)
  if (fig_type.eq."eps") then
  system("convert -trim -density 300 -crop 0x0 "+fig_name+".eps  "+fig_name+".png")
         ;==convert + quality(resolution) up + remove white boundary==;

  system("/bin/rm "+fig_name+".eps")
  end if
  print("Output Figure Name = "+fig_name+".png")


end
