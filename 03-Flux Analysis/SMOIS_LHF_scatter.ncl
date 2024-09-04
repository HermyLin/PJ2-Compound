load "/lfs/home/killumy/functions/func_PJ.ncl"

begin

ispv = -999
fspv = 1e+36

year_po = (/1980, 1983, 1993, 1998, 1999, 2003, 2014, 2015/)
year_ne = (/1981, 1984, 1985, 1990, 1994, 1997, 2001, 2013, 2018/)

year1d = array_append_record(year_po,year_ne,0)
qsort(year1d)

;5-8, JJAS
MM = ispan(5,8,1)
MM_name = "JJAS"
;==================================================
;read 
;==================================================
print("*****XAXIS*****")
  diri = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/OneF_Monthly/"
  fn   = "TReAD_2km_Monthly_LHF.nc"
  f = addfile(diri+fn,"r")

  ;JJAS
  LHF = f->LHF(:,MM,:,:)
  ;LHF = -1*LHF
  ;unit = W/m2

  delete([/fn,f/])

print("*****YAXIS*****")
  print("---read soil---")
  fn   = "TReAD_2km_Monthly_SMOIS.nc"
  f = addfile(diri+fn,"r")

  S1  = f->Soil1(:,MM,:,:)
  S2  = f->Soil1(:,MM,:,:)

  data = ((S1*0.07) + (S2*0.21)) / 0.28
  copy_VarCoords(S1,data)
  delete([/S1,S2/])

  lat = f->Lat
  lon = f->Lon
  nyear = dimsizes(f->Year)
  ;------------------------------------------
  ;mask height
  ;------------------------------------------
  H = 1000
  height = get_Taiwan_height(fspv,lat,lon)
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  scalarH = conform(data,height,(/2,3/))

  LHF  = LHF*scalarH
  soil = data*scalarH
  copy_VarCoords(LHF,soil)

  delete([/height,scalarH,data/])

  ;y    =  dim_avg_n_Wrap(LHF,(/1,2,3/))
  y    =  dim_avg_n_Wrap(LHF,(/2,3/))
  y_po = y({year_po},:)
  y_ne = y({year_ne},:)
  
  y_nu = y
  y_nu({year1d},:) = fspv
 
  ;x    = dim_avg_n_Wrap(soil,(/1,2,3/))
  x    = dim_avg_n_Wrap(soil,(/2,3/))
  x_po = x({year_po},:)
  x_ne = x({year_ne},:)

  x_nu = x
  x_nu({year1d},:) = fspv

  delete([/diri,fn,f/])
  delete([/LHF,soil/])
;==================================================
;regression
;==================================================
x0 = ndtooned(x_po)
x1 = ndtooned(x_ne)

y0 = ndtooned(y_po)
y1 = ndtooned(y_ne)


;positive
  rc_po = regline(x0,y0)
  line1 = rc_po*x0 + rc_po@yintercept
  print("Po:" + escorc(x0,y0))

;negative_<0.38
  x_ne1  = where(x1.le.0.38,x1,fspv)
  rc_ne1 = regline(x_ne1,y1)
  line2  = rc_ne1*x_ne1+ rc_ne1@yintercept
  print("Ne<=0.38:" + escorc(x_ne1,y1))
  delete(x_ne1)

;negative_>0.38
  x_ne2  = where(x1.gt.0.38,x1,fspv)
  rc_ne2 = regline(x_ne2,y1)
  line3  = rc_ne2*x_ne2+ rc_ne2@yintercept
  print("Ne>0.38:" + escorc(x_ne2,y1))
  delete(x_ne2)

  delete([/x_po,y_po,x_ne,y_ne/])
  ;exit
;==================================================
;plot
;frequency, HWV, HWMA, maxHWA, HWMD, HWMS, avgHWMS, HWMS
;==================================================
print("---plot---")
  fig_type = "X11"
  fig_type = "eps"
  fig_name = "TReADSMOIS_LHF_"+MM_name+"_Nu_038"

  wks = gsn_open_wks(fig_type,fig_name)

  res  = True
  res@vpWidthF  = 0.6
  res@vpHeightF = 0.4

  res@gsnDraw  = False
  res@gsnFrame = False

  res@xyMarkLineModes = (/"Markers","Lines"/)  ; choose which have markers
  res@xyMarker        = 1          ; choose type of marker 
  res@xyMarkerColor   = "red"      ; Marker color
  res@xyLineColor     = "red"      ; Marker color
  res@xyMarkerSizeF   = 0.055      ; Marker size (default 0.01)
  res@xyDashPattern   = 0         ; solid line 
  res@xyLineThicknesses = (/2,1/)  ; set second line to 2

  ;res@xyMarkLineMode  = "Markers"  ; choose which have markers
  ;res@xyDashPattern    = 1         ; solid line 
  ;res@xyLineThicknessF = 2.0      
  
  res@tiXAxisFont = 21
  res@tiXAxisFontHeightF = 0.025
  res@tiXAxisString = ""

  res@tiYAxisFont = 21
  res@tiYAxisFontHeightF = 0.025
  res@tiYAxisString = ""


  ;res@tmXBLabelFont = 22
  res@tmXBLabelFontHeightF = 0.023
  ;res@tmYLLabelFont = 22
  res@tmYLLabelFontHeightF = 0.023
  res@tmXBFormat = "0@;*.2f"

  ;res@tmXBMode = "Explicit"
  ;res@tmXBValues = fspan(30,32.5,6)
  ;res@tmXBLabels = "" + res@tmXBValues
  ;res@tmXBMinorValues = fspan(30,32.5,26)

  ;res@tmYLMode = "Explicit"
  ;res@tmYLValues = fspan(0,0.5,11)
  ;res@tmYLLabels = "" + res@tmYLValues
  ;res@tmYLMinorValues = fspan(0,0.5,51)

  res@trYMaxF = 130
  res@trYMinF = 75
  res@trXMinF = 0.28
  res@trXMaxF = 0.42

  res@tiMainString  = ""

  plot_po = gsn_csm_xy(wks,x0,(/y0,line1/),res)

  res@xyMarkLineModes := (/"Markers","Lines","Lines"/)  ; choose which have markers
  res@xyLineThicknesses := (/2,1,1/)                ; set second line to 2
  res@xyLineColor     = "blue"      ; Marker color
  res@xyMarkerColor   = "blue"      ; Marker color
  plot_ne = gsn_csm_xy(wks,x1,(/y1,line2,line3/),res)

  delete([/res@xyMarkLineModes,res@xyLineThicknesses/])

  res@xyMarkLineMode = "Markers"  ; choose which have markers
  res@xyMarker = 4         ; choose type of marker
  res@xyMarkerSizeF = 0.01       ; Marker size (default 0.01)
  res@xyMarkerThicknessF = 1.0
  res@xyMarkerColor = "gray"
  plot_nu = gsn_csm_xy(wks,x_nu,y_nu,res)

  overlay(plot_nu,plot_po)
  overlay(plot_nu,plot_ne)
  draw(plot_nu)

  ;draw(plot_po)
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
