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
nday = 92
ntim = nyear*nday

year_po = (/1980, 1983, 1993, 1998, 1999, 2003, 2014, 2015/)
year_ne = (/1981, 1984, 1985, 1990, 1994, 1997, 2001, 2013, 2018/)

year_PJ = array_append_record(year_po,year_ne,0)
qsort(year_PJ)

Po_all = 736.0
Ne_all = 828.0
Nu_all = 2116.0
;===================================================================================
;set rainfall threshold
;===================================================================================
  threshold = (/0,4,8,12,16,20,24,28,32,36,40,44,50,60,70,80,90,100,150,200,250,300/)

  bin_center = tofloat((/2,6,10,14,18,22,26,30,34,38,42,45,55,65,75,85,95,125,175,225,275,325/))
  Min_val = 1.9
;===================================================================================
;data_list
;===================================================================================
print("---read---")

  diri = "./"
  fn   = "TReAD_PJ_JJA_Rain_Freq_TWMean.nc"    ;Daily

  f = addfile(diri+fn,"r")

  pdf_po = tofloat(f->pdf_po) / Po_all
  pdf_ne = tofloat(f->pdf_ne) / Ne_all
  pdf_nu = tofloat(f->pdf_nu) / Nu_all

  ;pdf_po = tofloat(f->pdf_po) 
  ;pdf_ne = tofloat(f->pdf_ne)
  ;pdf_nu = tofloat(f->pdf_nu)

  ;pdf_po = 100*pdf_po
  ;pdf_ne = 100*pdf_ne
  ;pdf_nu = 100*pdf_nu
;===================================================================================
;PDF
;===================================================================================
print("---PDF setting---")
  nVar = 3
  nBin = dimsizes(pdf_po)

  xx = (/bin_center,bin_center,bin_center/)
  xx@_FillValue = fspv

  yy = (/pdf_po,pdf_ne,pdf_nu/)
  yy@_FillValue = fspv

  yy = where(yy.eq.0,fspv,yy)

;do i = nBin-9, nBin-2
;  print("Po: "+threshold(i)+"-"+threshold(i+1)+" : "+pdf_po(i))
;  print("Ne: "+threshold(i)+"-"+threshold(i+1)+" : "+pdf_ne(i))
;  print("Nu: "+threshold(i)+"-"+threshold(i+1)+" : "+pdf_nu(i))
;  print("----------------------------")
;end do
;===================================================================================
;plot
;===================================================================================
print("---plot---")
  fig_type = "X11"
  fig_type = "eps"

  fig_name = "PDF_Rainfall_TWMean"+"_Log"
  ;fig_name = "PDF_Rainfall_TWMean"+"_PDF"

  wks = gsn_open_wks(fig_type,fig_name)

  res  = True
  res@vpWidthF  = 0.7
  res@vpHeightF = 0.3

  res@gsnDraw  = False
  res@gsnFrame = False
  
  res@xyLineThicknessF         = 3.0        
  res@xyLineColors             = (/"red","blue","black"/)  
  res@xyMonoDashPattern        = True              ; all solid 

  res@trYLog = True
  res@trYAxisType = "LogAxis"
  ;res@trYMaxF = 1.0
  ;res@trYMinF = 0.0001

  res@trXLog = True
  res@trXAxisType = "LogAxis"

  ;res@tiYAxisString = "Frequency (%)"                ; y axis title
  res@tiYAxisString = "Frequency"                ; y axis title

  res@xyComputeYMin = True
  
  ;x-y
  res@tmXTOn = False
  res@tmYROn = False

  res@trXMinF = Min_val
  ;res@trXMaxF = 350

  ;res@tmXBMode = "Explicit"
  ;res@tmXBValues = threshold
  ;res@tmXBLabels = "" + res@tmXBValues

  res@tmXBMinorOn = False
  res@tmYLMinorOn = False

  res@tiXAxisFont = 21
  res@tiXAxisFontHeightF = 0.016
  res@tiYAxisFont = 21
  res@tiYAxisFontHeightF = 0.016

  res@gsnRightString = ""
  res@gsnLeftString  = ""
 
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
