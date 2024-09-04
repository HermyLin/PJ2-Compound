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
;===================================================================================
;data_list
;===================================================================================
print("---read---")
;1. read

  diri = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/Regrid_Daily/"
  fn   = "TReAD_2km_Daily_PRECTmms.nc"  ;monthly

  f = addfile(diri+fn,"r")

  ;upward postive
  data = f->PRECTmms
  data@_FillValue = fspv
  data@missing_value = fspv

  data = data*86400    ;mm/day
  ;[Year | 42] x [Day | 365] x [Lat | 173] x [Lon | 101]

  lat = f->Lat
  lon = f->Lon
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)

;2. get data
  data_JJA = data({year},151:242,:,:)
  delete(data)
;===================================================================================
;mask mountain site
;===================================================================================
print("---height---")
  H = 1000
  height = get_Taiwan_height(fspv,lat,lon)
  height@_FillValue = fspv
  height = where(height.gt.0.and.height.lt.H,1,fspv)

  scalarH = conform(data_JJA,height,(/2,3/))
  data = data_JJA*scalarH
  copy_VarCoords(data_JJA,data)
  data := dim_avg_n_Wrap(data,(/2,3/))

  delete([/data_JJA,height,scalarH/])
;===================================================================================
;address data
;===================================================================================
;1.separate PJ
  data_po = ndtooned(data({year_po},:))
  data_ne = ndtooned(data({year_ne},:))
  data_nu = data
  data_nu({year_PJ},:) = fspv

  data_nu := ndtooned(data_nu)
  delete(data)

  print("Po:" +num(.not.ismissing(data_po)))
  print("Ne:" +num(.not.ismissing(data_ne)))
  print("Nu:" +num(.not.ismissing(data_nu)))
  ;exit

;2. calculate frequency by bins
  ;set rainfall threshold
  threshold = (/0,4,8,12,16,20,24,28,32,36,40,44,50,60,70,80,90,100,150,200,250,300/)
  bins = dimsizes(threshold)

  pdf_po = new(bins,"integer",ispv)
  pdf_ne = new(bins,"integer",ispv)
  pdf_nu = new(bins,"integer",ispv)

  do i = 0, bins-2
     str_num = threshold(i)
     end_num = threshold(i+1)
     print(str_num+" <= counts < "+end_num)

     pdf_po(i) = num(data_po.ge.str_num.and.data_po.lt.end_num)
     pdf_ne(i) = num(data_ne.ge.str_num.and.data_ne.lt.end_num)
     pdf_nu(i) = num(data_nu.ge.str_num.and.data_nu.lt.end_num)

  end do

  pdf_po(bins-1) = num(data_po.gt.300)
  pdf_ne(bins-1) = num(data_ne.gt.300)
  pdf_nu(bins-1) = num(data_nu.gt.300)

;3. write out
  print("---output---")

  outfil = "./TReAD_PJ_JJA_Rain_Freq_TWMean.nc"
  system("rm -f "+outfil)
  fout = addfile(outfil,"c")

  fout->pdf_po = pdf_po
  fout->pdf_ne = pdf_ne
  fout->pdf_nu = pdf_nu

  fout@information = "not regular bins"
  print("output file: " + outfil)

end
