# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 16:04:16 2018

@author: admin
"""
import numpy as np
import trackingpy as tp
import math
CSPEED = 299792000.


#%%
#Test:
from  imp import reload
reload(tp)
year = 2019#2003
month = 8#10
day = 17#17
hour = 18#19
minute = 16#30
second = 39#30
lon_mos = 115.2505#-105.1786 # in degrees
lat_mos = 42.2118#39.742476 # in degrees
ele_mos = 1365#1830.14 # in meters
dT = 68 #67 seconds used in 2017.
el_sat=38.11 # 卫星的俯仰角
az_sat=195.42 #卫星的方位角
delt_sat = 15 # 卫星的赤纬


jd = tp.cal_jd(year,month,day,hour,minute,second)
dxyzs = tp.museri_ant_pos
ha_sun,dec_sun = tp.sun_hrdec(jd,dT,lon_mos,lat_mos,ele_mos) # 计算太阳的时角和赤纬
LXYZs = tp.enu2ecef(lat_mos,lon_mos,dxyzs)

wdly_sun = tp.wdly_cal(ha_sun,dec_sun,LXYZs)# 
wdly_p_sun = wdly_sun - wdly_sun.min() # in m
tdly_p_sun = wdly_p_sun/CSPEED


lat = 52
az = 283.271028
alt = 19.334344
ha_sat,delt_sat = tp.azalt2hadec(az,alt,lat)
wdly_sat = tp.wdly_cal(ha_sat,delt_sat,LXYZs)# 
wdly_p_sat = wdly_sat - wdly_sat.min() # in m
tdly_p_sat = wdly_p_sat/CSPEED


# Format dlys into fixed-point word

#%%
year=2009
month = 6
day =19
hour = 18
minute = 0
second = 0

jd1 = tp.cal_jd(year, month, day, hour, minute, second)
tp.ut2gst(1980,4,22,14,36,51.67)
tp.gst2lst(4.668119,-64)
#deg=(tdly_p[0]-tdly_p[1])*1000/1024*189*360*1e6
#print(deg)
#print(np.remainder(deg,-360))


