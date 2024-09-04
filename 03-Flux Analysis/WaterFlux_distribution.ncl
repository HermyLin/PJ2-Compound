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

JJA_DOY = ispan(5,8,1) ;JJAS: 5-8
;===================================================================================
;read
;vars: P, ET, R(SFR, UDR), SMSTOT
;===================================================================================
  diri = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/Regrid_Monthly/"

  VAR = (/"RAINNC","ACLHF"/)

  ;Rainfall
  ;------------------------------------------------------------------------
  list_in = diri+VAR(0)+"/"+"TReAD_monthly_"+year+"_"+VAR(0)+".nc"
  f       = addfiles(list_in,"r")
  ListSetType(f,"join")
  ;mm/hour
  
  RAIN = f[:]->$VAR(0)$(:,5:8,:,:)
  delete([/list_in,f/])
  ;------------------------------------------------------------------------
  ;LHF
  ;------------------------------------------------------------------------
  list_in = diri+VAR(1)+"/"+"TReAD_monthly_"+year+"_"+VAR(1)+".nc"
  f       = addfiles(list_in,"r")
  ListSetType(f,"join")
  
  LHF = f[:]->$VAR(1)$(:,5:8,:,:)
  delete([/list_in,f/])
  ;------------------------------------------------------------------------
  ;SMSTOT: need dt
  ;------------------------------------------------------------------------
  diri = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/Regrid_Daily/"
  list_in = diri+"SMSTOT/"+"TReAD_daily_"+year+"_SMSTOT.nc"
  f       = addfiles(list_in,"r")
  ListSetType(f,"join")
  
  SMSTOT = f[:]->SMSTOT

  SMSTOT_diff = RAIN
  ;Month: 6,7,8,9
  SMSTOT_diff(:,0,:,:) = SMSTOT(:,180,:,:) - SMSTOT(:,151,:,:)
  SMSTOT_diff(:,1,:,:) = SMSTOT(:,211,:,:) - SMSTOT(:,181,:,:)
  SMSTOT_diff(:,2,:,:) = SMSTOT(:,242,:,:) - SMSTOT(:,212,:,:)
  SMSTOT_diff(:,3,:,:) = SMSTOT(:,272,:,:) - SMSTOT(:,243,:,:)
  delete(SMSTOT)

  delete([/diri,list_in,f/])
  ;------------------------------------------------------------------------
  ;printVarSummary(RAIN)     ;mm/hr > *24
  ;printVarSummary(LHF)      ;J/m2 = 0.408*10^(-6) mm = 0.000000408 mm/hr > *24
  ;printVarSummary(SMSTOT)   ;kg/m2 >> mm
  
  ;MJ/m2/day = 0.408 mm /day
  ;28.94 W/m2 = 1 mm/day
  ;kg/m2 = mm
;===================================================================================
;organize data
;===================================================================================
  ;change units
  RAIN = RAIN*24 ;mm/day
  LHF  = LHF*0.000000408  ;mm/day?

  Residual = RAIN-LHF-SMSTOT_diff
  data = (/RAIN,LHF,SMSTOT_diff,Residual/)   
  data@_FillValue = fspv
  printVarSummary(data)
  exit
  copy_VarCoords(RAIN,data(0,:,:,:,:))
  delete([/RAIN,LHF,SMSTOT/])

  data!1 = "year"
  data&year = year

  data_JJA = dim_avg_n_Wrap(data(:,:,JJA_DOY,:,:),2)
  ;[4] x [year | 40] x [Lat | 173] x [Lon | 101]
  delete(data)
;===================================================================================
;mask height
;===================================================================================
  lat2km = data_JJA&Lat
  lon2km = data_JJA&Lon

  H = 1000
  height = get_Taiwan_height(fspv,lat2km,lon2km)
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  scalarH = conform(data_JJA,height,(/2,3/))

  data_JJA = data_JJA*scalarH
  data_all = dim_avg_n_Wrap(data_JJA,(/2,3/))

  delete([/H,height,scalarH/])
  delete([/lat2km,lon2km,data_JJA/])

  data_Mean = dim_avg_n_Wrap(data_all,1)
  printVarSummary(data_Mean)
;===================================================================================
;address data
;===================================================================================
;1.detrend
  data_dtr = dtrend_msg_n(data_all&year,data_all,True,False,1)
  copy_VarCoords(data_all,data_dtr)
  delete(data_all)
  
;2.anomaly on grids 
  data_clm = dim_avg_n_Wrap(data_dtr,1)  
  data_ano = data_dtr

  do yy = 0, nyear-1
     data_ano(:,yy) = data_dtr(:,yy) - data_clm
  end do

  data_all := data_ano
  delete([/data_clm,data_dtr,data_ano/])

  data_po  = dim_avg_n_Wrap(data_all(:,{year_po}),1)
  data_ne  = dim_avg_n_Wrap(data_all(:,{year_ne}),1)
;===================================================================================
;plot
;bar plot
;===================================================================================
print("---plot---")

  fig_type = "X11"
  fig_type = "eps"
  fig_name = "TReAD_JJAS_WaterBudget_distribution"
  ;fig_name = "TReAD_JJAS_WaterBudget_distribution_vertical"

  plot = new(4,graphic)

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
  res@tmXBValues      = (/0,1,2,3/) 
  res@tmXBLabels      = (/"PR","ET","SOIL","R+Ri"/)
  ;res@tmXBMinorValues = ispan(1960,2019,2)

  ;res@tmYLMode = "Manual"
  ;res@tmYLTickSpacingF = 1.0
  ;res@tmYLTickStartF   = -2.0
  ;res@tmYLTickEndF     = 2.0
  ;res@tmYLMinorOn = True

  res@tmYLLabelFont = 21
  res@tmYLLabelFontHeightF = 0.02 
  res@tmYLMinorOn = False

  res@tmXTOn = False
  res@tmYROn = False

  x = (/0,1,2,3/)
  res@trXMinF = -0.5
  res@trXMaxF = 3.5

  res@trYMinF = -1.0
  res@trYMaxF = 1.0
  plot(0) = gsn_csm_xy(wks,x,data_Mean,res)

  res@tmYLMode = "Manual"
  res@tmYLTickSpacingF = 0.1
  res@tmYLTickStartF   = -1.0
  res@tmYLTickEndF     = 1.0
  ;res@tmYLMinorOn = True

  res@trYMinF = -0.38
  res@trYMaxF = 0.38

  plot(2) = gsn_csm_xy(wks,x,data_po,res) ; create plot
  plot(3) = gsn_csm_xy(wks,x,data_ne,res) ; create plot

;panel
  resP = True
  ;resP@gsnMaximize = True
  ;resP@gsnPaperOrientation = "Portrait"
  ;resP@gsnPanelMainString         = "EOF mode1 of 1950-2021 JJA Precipitation"
  ;resP@gsnPanelMainFontHeightF    = 0.02
  ;resP@gsnPanelXWhiteSpacePercent = 4.
  ;resP@gsnPanelYWhiteSpacePercent = 3.

  gsn_panel(wks,plot,(/2,2/),resP) 
  ;gsn_panel(wks,plot(2:3),(/2,1/),resP) 
;---change figure type----------------------------------------------------------

  delete(wks)
  if (fig_type.eq."eps") then
  system("convert -trim -density 300 -crop 0x0 "+fig_name+".eps  "+fig_name+".png")
         ;==convert + quality(resolution) up + remove white boundary==;

  system("/bin/rm "+fig_name+".eps")
  end if
  print("Output Figure Name = "+fig_name+".png")


end
