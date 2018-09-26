def Td_from_spec_hum(spec_hum,p,t):
  # calculates Dewpoint temperature from specific humidity, temperature and pressure  
  Rd = 287.058 # gas constant for dry air [J * kg**-1 * K**-1]
  Rv = 461.5   # gas constant for water vapour [J * kg**-1 * K**-1]
  epsilon = Rd/Rv
  e = ((spec_hum/epsilon)*p)/(1+spec_hum/epsilon)  # vapor pressure in [hPa] ?
  T0 = 273.15 # [K]
  RvLv = 1.844*100**(-4) # [K**-1]
  e0 = 0.6113 # [KPa]
  Td = (1.0/T0 - RvLv * np.log((e*10.0)/e0))**(-1.0)
  return Td

def wet_bulb_potential_temperature(rec_t, rec_sp, rec_p):
  dp = 1 # [hPa]
  # find lifting condensation level
  found_LCL = False
  while(found_LCL == False):
    Td = Td_from_spec_hum(rec_sp,rec_p,rec_t)
    if(T <= TD):
      found = True
    else:
      rec_p = rec_p - dp
      dtdp = calc_dtdp_dry(rec_t, rec_p)
      rec_t = rec_t - (dp * dtdp)
  # bring down wet adiabatic to 1000 hPa
  rF = 0.0
  while(rec_p < 1000.):
    dtdp = calc_dtdp_wet(rec_t, rec_p, rF)
    rec_t = rec_t + (dp * dtdp)
    rec_p = rec_p + dp
  theta_wp = rec_t
  return theta_wp  
  
def Adedokun1(t850, sp850, t500):
  # calculates the Adedokun Index due to Haklander and VanDelden (2003)
  # input is t850: temperature [K] in 850 hPa
  # sp850: specific humidity[g/g] in 850 hPa
  # calculate wet bulb potential temperature in 850 hPa
  rec_t = t850
  rec_sp = sp850
  rec_p = 850.0
  theta_wp = wet_bulb_potential_temperature(t850, sp850, 850.0)
  # calculation of theta_sp
  rec_p = 500.0
  rec_t = t500
  while(rec_p < 1000.):
    dtdp = calc_dtdp_wet(rec_t, rec_p, rF)
    rec_t = rec_t + (dp * dtdp)
    rec_p = rec_p + dp
  theta_sp = rec_t
  res = theta_wp - theta_sp

#def Adedokun2()

def Boyden_Index(Z700, Z1000, T700):
  # calculates the Boyden Index
  # input variables are:
  # z700 and Z1000: geopotential heights of 700/1000 hPa level
  # T700 : temperature at 700 hPa in [K]
  res = 0.1 * (Z700-Z1000) - (T700 - 273.15) - 200
  return res    

def Bradbury_Index(t850, sp850, t500, sp500):
  # calculates the Bradbury Index due to Haklander and VanDelden (2003)
  # input is t850: temperature [K] in 850 hPa
  # sp850: specific humidity[g/g] in 850 hPa
  # t500: temperature [K] in 500 hPa
  # sp500: specific humidity[g/g] in 500 hPa
  # calculate wet bulb potential temperature in 850 hPa
  theta_w850 = wet_bulb_potential_temperature(t850, sp850, 850.0)
  theta_w500 = wet_bulb_potential_temperature(t500, sp500, 500.0)
  res = theta_w500 - theta_w850
  return res

def Cross_Totals_Index(t850, sp850, t500):
  # calculates the Cross Totals Index due to Haklander and VanDelden (2003)
  # input is t850: temperature [K] in 850 hPa
  # sp850: specific humidity[g/g] in 850 hPa
  # t500: temperature [K] in 500 hPa
  Td = Td_from_spec_hum(sp850,850.0,t850)
  res = Td - t500
  return res

 
#def Deep_Convective_Index()  
    
def Jefferson_Index(t850, sp850, t500, t700, sp700):
  # calculates the Jefferson Index due to Haklander and VanDelden (2003)
  # input is t850: temperature [K] in 850 hPa
  # sp850: specific humidity[g/g] in 850 hPa
  # t500: temperature [K] in 500 hPa
  # t700: temperature [K] in 700 hPa
  # sp700: specific humidity[g/g] in 700 hPa
  theta_w850 = wet_bulb_potential_temperature(t850, sp850, 850.0)
  Td700 =Td_from_spec_hum(sp700,700.0,t700)
  res = 1.6 * theta_w850 - t500 - 0.5 * (t700 - Td700) - 8.0
  return res
  
def K_Index(t850,sp850,t700,sp700,t500):
  # calculates the K Index due to Haklander and VanDelden (2003)
  # input is t850: temperature [K] in 850 hPa
  # sp850: specific humidity[g/g] in 850 hPa
  # t500: temperature [K] in 500 hPa
  # t700: temperature [K] in 700 hPa
  # sp700: specific humidity[g/g] in 700 hPa
  Td850 =Td_from_spec_hum(sp850,850.0,t850)
  Td700 =Td_from_spec_hum(sp700,700.0,t700)
  res = (t850 - t500) + Td850 - (t700 - Td700)
  return res
  
def Modified_K_Index(ts850,sps850,ps850,t700,sp700,t500):
  # calculates the modified K Index due to Haklander and VanDelden (2003)
  # input is ts850: temperature array [K] from surface to 850 hPa
  # sps850: specific humidityarray [g/g] from surface to 850 hPa
  # t500: temperature [K] in 500 hPa
  # t700: temperature [K] in 700 hPa
  # sp700: specific humidity[g/g] in 700 hPa
  import numpy as np
  t_mean = lnp_averaging(ps850, ts850)
  n = len(ts850)
  tds850 = []
  for i in range(0, n):
    td_temp =Td_from_spec_hum(sps850[i],ps850[i],ts850[i])
    tds850.append(td_temp)
  td_mean = lnp_averaging(ps850, tds850)  
  Td700 =Td_from_spec_hum(sp700,700.0,t700)
  res = (t_mean - t500) + td_mean - (t700 - Td700)    
  return res
  
def lnp_averaging(p, v):  
  # calculates the ln(p) averaging used in Haklander and VanDelden (2003)
  # input is:
  # p = pressure array [hPa]
  # v = variable to average
  import numpy as np    
  n = len(p)
  sum_v = 0.0
  sum_lnp = 0.0
  for i in range(0, n):
    sum_v = sum_v + np.log(p[i]) * v[i]
    sum_lnp = sum_lnp + np.log(p[i])
  v_mean = sum_v / sum_lnp
  return v_mean

  
def KO_Index(T500,sp500,T700,sp700,T850,sp850,T1000,sp1000):
  # calculates the KO Index due to Haklander and VanDelden (2003)
  # input is 
  # T500: temperature [K] in 500 hPa
  # sp500: specific humidity [g/g] in 500 hPa
  # T700: temperature [K] in 700 hPa
  # sp700: specific humidity [g/g] in 700 hPa
  # T850: temperature [K] in 850 hPa
  # sp850: specific humidity [g/g] in 850 hPa
  # T1000: temperature [K] in 1000 hPa
  # sp1000: specific humidity [g/g] in 1000 hPa
  
  thetae_500 = theta_e(500,T500,sp500)
  thetae_700 = theta_e(700,T700,sp700)
  thetae_850 = theta_e(850,T850,sp850)
  thetae_1000 = theta_e(1000,T1000,sp1000)
  res = 0.5 * (thetae_500 + thetae_700) - 0.5 * (thetae_850 + thetae_1000)
  return res
  
def parcel_ascent(p,t,sp,final_pressure):
  # calculates dry and subsequent wet adiabatic ascent up to the 'final_pressure' level 
  # with intial properties p, t, and sp
  from thermodynamic_routines import saturation_vapour_pressure
  epsilon = 0.622 #[g/g]  
  dp=1.0
  rec_p = p
  rec_t = t
  rec_sp = sp
  found_LCL = False
  while((found_LCL == False) and (rec_p > final_pressure)):
    Td = Td_from_spec_hum(rec_sp,rec_p,rec_t)
    if(rec_t <= Td):
      found_LCL = True
    else:
      rec_p = rec_p - dp
      dtdp = calc_dtdp_dry(rec_t, rec_p)
      rec_t = rec_t - (dp * dtdp)
  rF = 0.0
  while(rec_p > final_pressure):
    dtdp = calc_dtdp_wet(rec_t, rec_p, rF)
    rec_t = rec_t + (dp * dtdp)
    es = saturation_vapour_pressure(rec_T)
    rec_sp = (epsilon * es) / (rec_p - es * (1.0 - epsilon)) 
    rec_p = rec_p + dp
  return rec_t, rec_sp  

def Lifted_Index(p,t,sp,t500):
  # calculates the Lifted Index due to Haklander and VanDelden (2003)
  # input is
  # p: pressure of initial parcel [hPa]
  # t: temperature of initial parcel [K]
  # sp: specific humidity of initial parcel in [g/g]
  dp=1.0
  rec_p = p
  rec_t = t
  rec_sp = sp
  res_t, res_sp = parcel_ascent(rec_p,rec_t,rec_sp,500.0)
  res = t500 - res_t
  return res
  
def Lifted_Index_layer(p_array,t_array,sp_array,t500):
  # calculates the Lifted Index 50 hPa or 100 hPa due to Haklander and VanDelden (2003)
  # input is
  import numpy as np
  rec_t  = lnp_averaging(p_array, t_array)
  rec_sp  = lnp_averaging(p_array, sp_array)
  rec_p = np.min(p_array)
  res = Lifted_Index(rec_p,rec_t,rec_sp,t500)
  return res

def Lifted_Index_mu(p_array,t_array,sp_array,t500):
  # calculates the most unstable Lifted Index (ThetaE is max at 250 hPa above surface) due to Haklander and VanDelden (2003)
  # input is
  import numpy as np
  n = len(p_array)
  max_thetae = 0.0
  final_pressure = np.min(p_array) - 250.0
  for i in range(0, n):
    res_t, res_sp = parcel_ascent(p_array[i],t_array[i],sp_array[i],final_pressure)
    thetae_temp = theta_e(final_pressure,res_t,res_sp)
    if(thetae_temp > max_thetae):
      res = Lifted_Index(p_array[i],t_array[i],sp_array[i])
  return res

def Potential_Instability_Index(t925,sp925,Z925,t500,sp500,Z500):
  # calculates the Potential Instability Index due to Haklander and VanDelden (2003)
  # input is
  # Z925 and Z500: the geopotential height of the 925 hPa and the 500 hPa level
  thetae_925 = theta_e(925.0,t925,sp925)
  thetae_500 = theta_e(500.0,t500,sp500)
  res = (thetae_925 - thetae_500)/(Z500 - Z925)
  return res
  
def Rackliff_Index(t900,sp900,t500):
  # calculates the Rackliff Index due to Haklander and VanDelden (2003)
  theta_w_900 =  wet_bulb_potential_temperature(t900, sp900, 900.0)  
  res = theta_w_900 - t500
  return res

def Showalter_Index(t850,sp850,t500):
  # calculates the Showalter Index due to Haklander and VanDelden (2003)
  rec_t, rec_sp = parcel_ascent(850.0,t850,sp850,500.0)
  SSI = t500 - rec_t
  return SSI

def S_Index(t850,sp850,t700,sp700):
  # calculates the S-Index due to Haklander and VanDelden (2003)
  TT = Total_Totals_Index(t850,sp850,t500)
  VT = Vertical_Totals(t850,t500)
  if(VT > 25.9):
    A = 0.0
  else:
    if(VT >= 22.0):
      A = 2.0
    else:
      A = 6.0
  td700 = Td_from_spec_hum(sp700,700.0,t700)
  SI = TT - (t700 - td700) - A
  return SI
  
def SWEAT_Index(t850,sp850,t500,ws850,wd850,ws500,wd500):
  # calculates the SWEAT Index due to Haklander and VanDelden (2003)
  # be aware that the equation in Haklander and VanDelden (2003) is wrong ! 
  # input is
  # t850: temperature in 850 hPa [K]
  # sp850: specific humidity in 850 hPa [g/g]
  # t500: temperature in 500 hPa [K]
  # ws850: wind speed in 850 hPa [m/s]
  # wd850: wind direction in 850 hPa [0 - 360 deg] 
  # ws500: wind speed in 500 hPa [m/s]
  # wd500: wind direction in 500 hPa [0 - 360 deg]
  
  Td850 =Td_from_spec_hum(sp850,850.0,t850)
  Td850 = Td850 - 273.15
  TT = Total_Totals_Index(t850,sp850,t500)
  ws850_kn = ws850 * 1.94384
  ws500_kn = ws500 * 1.94384
  TT_term = 20 * (TT - 49.0)
  if(TT_term < 0.0):
    TT_term = 0.0
  wd_term = 0.0
  if(wd500 >= 210.0):
    if(wd500 <= 310.0):
      if(wd850 >= 130.0):
        if(wd850 <= 250.0):
          if(wd500 > wd850):
            if(ws500_kn >= 15.0):
              if(ws850_kn >= 15.0):
                wd_term =  125.0 * (math.sin(wd500-wd850)+0.2)  
  res = 12.0 * Td850  + TT_term + 2.0 * ws850_kn + ws500_kn + wd_term
  return res

def SWISS00_Index(t850,sp850,t600,sp600,t500,ws_3km,wd_3km,ws_6km,wd_6km):
  # calculates the SWISS00-Index  due to Haklander and VanDelden (2003)
  # input is:
  # t850: temperature in 850 hPa [K]
  # sp850: specific humidity in 850 hPa [g/g]
  # t600: temperature in 850 hPa [K]
  # sp600: specific humidity in 850 hPa [g/g]
  
  # t500: temperature in 500 hPa [K]
  # ws_3km: windspeed in 3 km above the ground [m/s]
  # wd_3km: wind direction in 3 km above the ground [0 - 360 deg]
  import math
  
  SHOW = Showalter_Index(t850,sp850,t500)  
  td600 =Td_from_spec_hum(sp600,600.0,t600)
  u_3km = -1.0 * ws_3km * math.sin(wd_3km)
  v_3km = -1.0 * ws_3km * math.cos(wd_3km)
  u_6km = -1.0 * ws_3km * math.sin(wd_3km)
  v_6km = -1.0 * ws_3km * math.cos(wd_3km)
  d_u = u_6km - u_3km
  d_v = v_6km - v_3km
  wsh = math.sqrt(d_u**2.0 + d_v**2.0)
  res = SHOW + 0.4 * wsh + 0.1 * (t600 - td600)
  return res

def SWISS12_Index(p_sfc,t_sfc,sp_sfc,t650,sp650,t500,ws_10m,wd_10m,ws_3km,wd_3km):
  # calculates the SWISS12-Index  due to Haklander and VanDelden (2003)
  # input is:
  # t850: temperature in 850 hPa [K]
  # sp850: specific humidity in 850 hPa [g/g]
  # t600: temperature in 850 hPa [K]
  # sp600: specific humidity in 850 hPa [g/g]
  
  # t500: temperature in 500 hPa [K]
  # ws_3km: windspeed in 3 km above the ground [m/s]
  # wd_3km: wind direction in 3 km above the ground [0 - 360 deg]
  import math
  
  LI_sfc =  Lifted_Index(p_sfc,t_sfc,sp_sfc,t500)  
  td650 =Td_from_spec_hum(sp650,650.0,t650)
  u_10m = -1.0 * ws_10m * math.sin(wd_10m)
  v_10m = -1.0 * ws_10m * math.cos(wd_10m)
  u_3km = -1.0 * ws_3km * math.sin(wd_3km)
  v_3km = -1.0 * ws_3km * math.cos(wd_3km)
  d_u = u_3km - u_10m
  d_v = v_3km - v_10m
  wsh = math.sqrt(d_u**2.0 + d_v**2.0)
  res = LI_sfc - 0.1 * wsh+ 0.1 * wsh + 0.1 * (t650 - td650)
  return res

def Thompson_Index(t50,p50,sp50,t850,sp850,t700,sp700,t500):
  # calculates the Thompson Index due to Haklander and VandDelden (2003)
  # t50: temperatures of the layer 50 hPa above surface 

  t_mean = lnp_averaging(p50, t50)
  sp_mean = lnp_averaging(p50, sp50)
  p_mean = max(p50)
  LI50 =  Lifted_Index(p_mean,t_mean,sp_mean,t500)
  KI = K_Index(t850,sp850,t700,sp700,t500)
  THOM = KI - LI50
  return THOM

def Total_Energy_Index(t500,sp500,z500,p_test,t_test,sp_test,z_test):
  # calculates the Total Energy Index due to Haklander and VandDelden (2003)
  # input is:
  # t_test: temperature in desired pressure level [K]
  # sp_test: specific humidity in desired pressure level [g/g]
  # z_test: geopotential height of desired pressure level [gpm]
  # t500: temperature in 500 hPa [K]
  # sp500: specific humidity in 500 hPa [g/g]
  # z500: geopotential height of 500 hPa level [gpm]

  r500 = mixing_ratio_from_specific_humidity(sp500,500.0)
  h500 = static_energy(r500, t500, z500)
  r_test = mixing_ratio_from_specific_humidity(sp_test,p_test)
  h_test = static_energy(r_test,t_test,z_test)
  TEI = h500 - h_test
  return TEI

def static_energy(r, T, Z):
  # calculates the static energy due to Haklander and VandDelden (2003)
  # input is:
  # r: mixing ratio [g/g]
  # T: temperature [K]
  # Z: geopotential height [m]
  cpd = 1004.0 # J / (kg * K) <- heat capacity of dry air at constant pressure
  cl = 4217.6  # J / (kg * K) <- heat capacity of liquid water
  Lv = 2.501 * 10.0**6  # [J/kg] <- latent heat of vaporization
  g = 9.81 # [m/s**2]
  h = (cpd + r * cl) * T + Lv * r + (1.0 + r)* g * Z 
  return h

def Total_Totals_Index(t850,sp850,t500):
  # calculates the Total Totals Index due to Haklander and VandDelden (2003)
  # original: Miller (1967)
  CT = Cross_Totals_Index(t850,sp850,t500)
  VT = Vertical_Totals(t850,t500)
  TT = VT + CT
  return TT

def Vertical_Totals(t850,t500):
  # calculates the Vertical Totals Index due to Haklander and VandDelden (2003)
  # input is:
  # t850: temperature in 850 hPa [K]
  # t500: temperature in 500 hPa [K]
  VT = t850 - t500
  return VT

def Yonetani_Index():
  # calculates the Yonetani Index due to Haklander and VandDelden (2003)

  YON = 0.966 * Gamma_L + 2.41 * (Gamma_U - Gamma_W)
  return YON

def vapour_pressure_from_specific_humidity(q,p):
  # derives vapour pressure from specific humidity
  # equation derived from Stull: Meteorology for Scientists and Engineers  
  # input is:
  # q: specific humidity [g_water_vapour / g_total]    
  # p: pressure in hPa
  # output is vapour pressure in [hPa]
  
  epsilon = 0.622 #[g/g]  
  e = (q*p) / (epsilon - q*(1.0 - epsilon))
  return e
  
def mixing_ratio_from_specific_humidity(q,p):
  # derives vapour pressure from specific humidity
  # equation derived from Stull: Meteorology for Scientists and Engineers    
  # input is:
  # q: specific humidity [g_water_vapour / g_total]    
  # p: pressure in hPa
  # output is mixing ration in [g/g]
  e = vapour_pressure_from_specific_humidity(q,p)
  epsilon = 0.622 #[g/g]
  r = (epsilon*e)/(p-e)
  return r

def specific_humidity_from_mixing_ratio(r, p):
  # derives specific humidity from mixing ratio
  # equation derived from Stull: Meteorology for Scientists and Engineers
  # input is:
  # r: mixing ratio [g_water_vapour / g_dry]    
  # p: pressure in hPa
  # output is specific humidity in [g/g]
  epsilon = 0.622 #[g/g]   
  e = (r*p)/(epsilon+r) 
  q = (epsilon*e)/(p-e*(1.0-epsilon))
  return q
  
def theta_e(p,T,sp):
  import numpy as np 
  import math
  # calculates theta e based on the equation of the 
  # Bolton (1980) Eq. 43 and Bryan (2008) Eq.6
  # input is:
  # p: pressure [hPa]
  # T: temperature [K]
  # sp: specific humidity [g/g]
    
  e = vapour_pressure_from_specific_humidity(sp,p)
  r = mixing_ration_from_specific_humidity(sp,p)
  TL = (2840.0 / (3.5 * np.log(T) - np.log(e) - 4.805)) + 55.0
  exp1 = 0.2854*(1.0-0.28 * r)    
  theta_e = T * (1000.0/p)**exp1 * math.exp((3376.0/TL - 2.54) * r * (1.0 + 0.81 *r))
  return theta_e

rec_p = 996.0
rec_T = 29.4 + 273.15
rec_r = 12.71/1000.0
rec_sp = specific_humidity_from_mixing_ratio(rec_r, rec_p)
thetaE = theta_e(rec_p,rec_T,rec_sp)
print(thetaE)
    
    
