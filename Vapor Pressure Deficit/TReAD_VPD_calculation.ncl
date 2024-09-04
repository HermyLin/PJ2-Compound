begin

ispv = -999
fspv = 1e+36

;==========================================
;calculate JJA VPD
;
;SVP = 610.78 x e^(T / (T +237.3) x 17.2694)
;T in degC, SVP in Pa
;VPD = SVP x (1 – RH/100)
;==========================================

diri = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/Daily/"

;==========================================
;1. read Temeprature (DegC) & calculate SVP
;==========================================
print("---SVP---")
  f = addfile(diri+"TReAD_2km_Daily_TBOT.nc","r")

  T = (f->TBOT) - 273.15
  ;unit: K > DegC

  year = f->Year
  ;day  = f->Day
  lat  = f->Lat
  lon  = f->Lon
  ;----------------------------------------
  print("---SVP---")
  SVP =  610.78 * exp(17.2694*(T/(T+237.3)))
  SVP@_FillValue = fspv

  delete([/f,T/])
;==========================================
;2. read RH & calculate VPD
;==========================================
print("---VPD---")
  f = addfile(diri+"TReAD_2km_Daily_RH.nc","r")
  RH = (f->RH) / 100
  VPD = SVP * (1-RH)

  delete([/f,RH/])

  SVP = nameDim(SVP,(/"Year","Day","Lat","Lon"/),"","") 
  SVP&Year = year
  ;SVP&Day  = day
  SVP&Lat  = lat
  SVP&Lon  = lon
  SVP@_FillValue = fspv

  copy_VarCoords(SVP,VPD)
  delete([/year,lat,lon/])
;==========================================
;3. output
;==========================================
  dir_out = "/lfs/home/killumy/SYLIN_DATABASE/TReAD_2km/Daily/"

  outfil = dir_out + "TReAD_2km_Daily_VPD.nc"
  system("rm -f "+outfil)

  setfileoption("nc","Format","NetCDF4Classic")
  fout = addfile(outfil,"c")
  fout->VPD = VPD
  fout->SVP = SVP

  print("output file: " + outfil)

end
