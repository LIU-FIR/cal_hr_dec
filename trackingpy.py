# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 10:28:45 2021
包括太阳的(hour_angle,declination)计算
地心坐标系与地平面坐标系转换
空间频率分量[u,v,w]的计算
延时的计算

@author: LIU FEI
"""

import numpy as np
from numpy import pi
#from JD_cal import cal_jd
import math
CSPEED = 299792000.
#dT = 68 #67 seconds used in 2017.
#%%
## Load variables
L0_tab = np.load('L0.npy')
L1_tab = np.load('L1.npy')
L2_tab = np.load('L2.npy')
L3_tab = np.load('L3.npy')
L4_tab = np.load('L4.npy')
L5_tab = np.load('L5.npy')   
B0_tab = np.load('B0.npy')
B1_tab = np.load('B1.npy')
R0_tab = np.load('R0.npy')
R1_tab = np.load('R1.npy')
R2_tab = np.load('R2.npy')
R3_tab = np.load('R3.npy')    
R4_tab = np.load('R4.npy') 
TA43 = np.load('TA43.npy')
museri_ant_pos = np.load('museri_ant_pos.npy')

def sun_hrdec(JD,lon,lat,ele,dT=68):
    """
    计算太阳的指向：时角-赤纬
    输入：year,month, day, hour,minute,second.
    输出： hour-angle(小时)，declination(角度)
    """        
##   Step1: JD, JDE, JC, JCE, JME    
    JDE = JD + dT/86400
    JC = (JD - 2451545)/36525
    JCE = (JDE-2451545)/36525
    JME = JCE/10    
##   Step2: L, B, R (Earth heliocentric longitude,lattitude,radius vector)

    L0 = np.sum(L0_tab[0:,0:1]*np.cos(L0_tab[0:,1:2]+L0_tab[0:,2:3]*JME))
    L0t=0;
    for k in range(64):
        L0t=L0t+L0_tab[k][0]*math.cos(L0_tab[k][1]+L0_tab[k][2]*JME)
    L1 = np.sum(L1_tab[0:,0:1]*np.cos(L1_tab[0:,1:2]+L1_tab[0:,2:3]*JME))
    L2 = np.sum(L2_tab[0:,0:1]*np.cos(L2_tab[0:,1:2]+L2_tab[0:,2:3]*JME))
    L3 = np.sum(L3_tab[0:,0:1]*np.cos(L3_tab[0:,1:2]+L3_tab[0:,2:3]*JME))
    L4 = np.sum(L4_tab[0:,0:1]*np.cos(L4_tab[0:,1:2]+L4_tab[0:,2:3]*JME))
    L5 = np.sum(L5_tab[0:,0:1]*np.cos(L5_tab[0:,1:2]+L5_tab[0:,2:3]*JME))    
    Lrad = (L0+L1*JME+L2*JME**2+L3*JME**3+L4*JME**4+L5*JME**5)/1e8
    Ldeg = Lrad*180/pi # in degrees
    LF = math.modf(Ldeg/360)[0] # np.modf(v) returns decimal & integral part
    # if LF is always positive and L is negative, L = 360-360*F    
    L = 360*LF if Ldeg>=0 else 360+360*LF 

    B0 = np.sum(B0_tab[0:,0:1]*np.cos(B0_tab[0:,1:2]+B0_tab[0:,2:3]*JME))
    B1 = np.sum(B1_tab[0:,0:1]*np.cos(B1_tab[0:,1:2]+B1_tab[0:,2:3]*JME))
    Brad = (B0+B1*JME)/1e8
    B = Brad*180/pi # in degrees
  
    R0 = np.sum(R0_tab[0:,0:1]*np.cos(R0_tab[0:,1:2]+R0_tab[0:,2:3]*JME))
    R1 = np.sum(R1_tab[0:,0:1]*np.cos(R1_tab[0:,1:2]+R1_tab[0:,2:3]*JME))
    R2 = np.sum(R2_tab[0:,0:1]*np.cos(R2_tab[0:,1:2]+R2_tab[0:,2:3]*JME))
    R3 = np.sum(R3_tab[0:,0:1]*np.cos(R3_tab[0:,1:2]+R3_tab[0:,2:3]*JME))
    R4 = np.sum(R4_tab[0:,0:1]*np.cos(R4_tab[0:,1:2]+R4_tab[0:,2:3]*JME))    
    R = (R0+R1*JME+R2*JME**2+R3*JME**3+R4*JME**4)/1e8 #in AU
    ## Step 3: THETA & BETA(geocentric longitude and latitude)
    THETA_g = L+180
    # Limit THETA to the range from 0 to 360 degree.
    TF = math.modf(THETA_g/360)[0]
    THETA = 360*TF if THETA_g>=0 else 360+360*TF # in degrees
    BETA = -B # in degrees
    
    ## Step 4: dPSI & dEPS (the nutation in longitude & obliquity)

    # X0:The mean elongation of the moon from the sun (in degrees)
    X0 = 297.85036 + 445267.111480*JCE - 0.0019142*JCE**2 + JCE**3/189474
    # X1: the mean anomaly of the sun (in degrees)
    X1 = 357.52772 + 35999.050340*JCE - 0.0001603*JCE**2 - JCE**3/300000
    # X2: the mean anomaly of the moon (in degrees)
    X2 = 134.96298 + 477198.867398*JCE + 0.0086972*JCE**2 +JCE**3/56250
    # X3: the moon's argument of latitude (in degrees)
    X3 = 93.27191 + 483202.017538*JCE - 0.0036825*JCE**2 + JCE**3/327270
    # X4: the longitude of the ascending node of the moon's mean orbit on the ecliptic,measured from the mean equinox of the date (in degrees)
    X4 = 125.04452 - 1934.136261*JCE + 0.0020708*JCE**2 + JCE**3/450000
    Xv = np.array([X0,X1,X2,X3,X4])*pi/180
    dPSI = np.sum((TA43[0:63,5:6]+TA43[0:63,6:7]*JCE).reshape(63,)*np.sin(np.matmul(TA43[0:63,0:5],Xv)))/36000000 #in degrees
    dEPS = np.sum((TA43[0:63,7:8]+TA43[0:63,8:9]*JCE).reshape(63,)*np.cos(np.matmul(TA43[0:63,0:5],Xv)))/36000000 #in degrees
#    dPSI = 
    
    ## Step 5: EPS (the true obliquity of the ecliptic, in degrees)
    # EPS0 (the mean obliquity of the ecliptic, in arc seconds)
    U = JME/10
    EPS0 = 84381.448 - 4680.93*U - 1.55*U**2 + 1999.25*U**3 - \
        51.38*U**4 - 249.67*U**5 - 39.05*U**6 + 7.12*U**7 + \
        27.87*U**8 + 5.79*U**9 + 2.45*U**10
    EPS = EPS0/3600 + dEPS
    
    ## Step 6:dTAU (the aberration correction, in degrees)
    dTAU = -20.4898/(3600*R)
    
    ## Step 7: LAM (the apparent sun longitude, in degrees)
    LAM = THETA + dPSI + dTAU
    
    ## Step 8: VU (the apparent sidereal time at Greenwich at any given time: in degrees)
    VU0_deg = 280.46061837 + 360.98564736629*(JD-2451545) + \
            0.000387933*JC**2 - JC**3/38710000
    VU0F = math.modf(VU0_deg/360)[0]
    VU0 = 360*VU0F if VU0_deg>=0 else 360+360*VU0F # in degrees
    VU = VU0 + dPSI*math.cos(EPS*pi/180)
    
    ## Step 9: ALP (the geocentric sun right ascension, in degrees)
    ALP_deg = math.atan2((math.sin(LAM*pi/180)*math.cos(EPS*pi/180)-
                          math.tan(BETA*pi/180)*math.sin(EPS*pi/180)),
                            math.cos(LAM*pi/180))*180/pi
    ALPF = math.modf(ALP_deg/360)[0]
    ALP = 360*ALPF if ALP_deg>=0 else 360+360*ALPF
    
    ## Step 10: DELT (the geocentric sun declination, in degrees)
    DELT = math.asin(math.sin(BETA*pi/180)*math.cos(EPS*pi/180) +
                     math.cos(BETA*pi/180)*math.sin(EPS*pi/180)*math.sin(LAM*pi/180))*180/pi
                     
    ## Step 11: H (the ovserver local hour angle, in degrees)                     
    H_deg = VU + lon - ALP
    HF = math.modf(H_deg/360)[0]
    H = 360*HF if H_deg>=0 else 360+360*HF

    return (H,DELT)    
#    ## Step 12: ALP_pl (the topocentric sun right ascension, in degrees)
#    XI = 8.794/(3600*R) # in degrees
#    u = math.atan(0.99664719*math.tan(lat*pi/180)) # in radians
#    x = math.cos(u) + ele/6378140*math.cos(lat*pi/180)
#    y = 0.99664719*math.sin(u) + ele/6378140*math.sin(lat*pi/180)
#    dALP = math.atan2(-x*math.sin(XI*pi/180)*math.sin(H*pi/180),
#                      math.cos(DELT*pi/180)-x*math.sin(XI*pi/180)*math.cos(H*pi/180))*180/pi
#    ALP_pl = ALP + dALP           
#    # DELT_pl: topocentric sun declination, in degrees
#    DELT_pl = math.atan2((math.sin(DELT*pi/180)-y*math.sin(XI*pi/180))*math.cos(dALP*pi/180),
#                         math.cos(DELT*pi/180)-x*math.sin(XI*pi/180)*math.cos(H*pi/180))*180/pi           
#
#    ## Step 13: H_pl (the topocentric local hour angle, in degrees)                         
#    H_pl = H- dALP

#    return (H_pl, DELT_pl)    
#%%
def cal_jd(year, month, day, hour, minute, second):
    """
    计算Julian day
    输入：year,month, day, hour,minute,second.
    """        
    Y = year if month > 2 else year - 1
    M = month if month > 2 else month + 12
    D = day + hour/24 + minute/(24*60) + second/(24*3600)
    A = math.floor(Y/100)
    Bz = 2 - A + np.floor(A/4)
#    Bz = 0
    JD = math.floor(365.25*(Y + 4716)) + math.floor(30.6001*(M+1)) + D + Bz - 1524.5    
    return JD

def ut2gst(year, month, day, hour, minute, second,dT=68):
    """
    根据UT计算GST
    输入：year,month,day(Greenwich date), hour,minute,second(UT).
    输出：以hour为单位
    """    
    JD = cal_jd(year, month, day, hour, minute, second)
    JDE = JD + dT/86400
    JC = (JD - 2451545)/36525
    JCE = (JDE-2451545)/36525
    JME = JCE/10  
    # X0:The mean elongation of the moon from the sun (in degrees)
    X0 = 297.85036 + 445267.111480*JCE - 0.0019142*JCE**2 + JCE**3/189474
    # X1: the mean anomaly of the sun (in degrees)
    X1 = 357.52772 + 35999.050340*JCE - 0.0001603*JCE**2 - JCE**3/300000
    # X2: the mean anomaly of the moon (in degrees)
    X2 = 134.96298 + 477198.867398*JCE + 0.0086972*JCE**2 +JCE**3/56250
    # X3: the moon's argument of latitude (in degrees)
    X3 = 93.27191 + 483202.017538*JCE - 0.0036825*JCE**2 + JCE**3/327270
    # X4: the longitude of the ascending node of the moon's mean orbit on the ecliptic,measured from the mean equinox of the date (in degrees)
    X4 = 125.04452 - 1934.136261*JCE + 0.0020708*JCE**2 + JCE**3/450000
    Xv = np.array([X0,X1,X2,X3,X4])*pi/180
    dPSI = np.sum((TA43[0:63,5:6]+TA43[0:63,6:7]*JCE).reshape(63,)*np.sin(np.matmul(TA43[0:63,0:5],Xv)))/36000000 #in degrees
    dEPS = np.sum((TA43[0:63,7:8]+TA43[0:63,8:9]*JCE).reshape(63,)*np.cos(np.matmul(TA43[0:63,0:5],Xv)))/36000000 #in degrees

    U = JME/10
    EPS0 = 84381.448 - 4680.93*U - 1.55*U**2 + 1999.25*U**3 - \
        51.38*U**4 - 249.67*U**5 - 39.05*U**6 + 7.12*U**7 + \
        27.87*U**8 + 5.79*U**9 + 2.45*U**10
    EPS = EPS0/3600 + dEPS   
   
    ## Step 8: VU (the apparent sidereal time at Greenwich at any given time: in degrees)
    VU0_deg = 280.46061837 + 360.98564736629*(JD-2451545) + \
            0.000387933*JC**2 - JC**3/38710000
    VU0F = math.modf(VU0_deg/360)[0]
    VU0 = 360*VU0F if VU0_deg>=0 else 360+360*VU0F # in degrees
    VU = VU0 + dPSI*math.cos(EPS*pi/180)    
    gst = (VU/15) % 24
    return gst

#   以下计算版本：Practical Astronomy with your Calculator or Spreadsheet by Duffett-Smith P., Zwart J
#def ut2gst1(year, month, day, hour, minute, second):
#    """
#    根据UT计算GST
#    以下计算版本：Practical Astronomy with your Calculator or Spreadsheet
#    """        
#    JD = cal_jd(year, month, day, 0, 0, 0)
#    S = JD-2451545
#    T = S/36525
#    T0 = 6.697374558 + (2400.051336*T) + (0.000025862*T**2)
#    T0 = T0 % 24
#    UT = (second/60+minute)/60+hour
#    A = UT*1.002737909
#    GST = (A+T0) % 24
#    return GST

def gst2lst(gst,local_lon):
    """
    根据GST计算LST
    输入：GST(以hour为单位), local_lon(以度为单位)
    输出：LST(以hour为单位)
    """    
    LST = (gst+local_lon/15) % 24
    return LST

def ra2ha(ra,lst):
    """
    根据ra(right ascension)计算hour angle
    输入：ra(以hour为单位), lst(以度为单位)
    输出：LST(以hour为单位)
    """    
    ha = lst - ra
    if ha<0:
        ha = ha +24
    return ha

def azalt2hadec(az,alt,lat):
    """
    已知方位角、俯仰角和观测者纬度，计算时角和赤纬
    输入：az(方位角),alt(俯仰角)，lat(赤纬)，以角度为单位
    """           
    alt_rad = alt/180*np.pi
    az_rad = az/180*np.pi
    lat_rad = lat/180*np.pi
    delt_rad=math.asin(math.sin(alt_rad)*math.sin(lat_rad)+math.cos(az_rad)*math.cos(lat_rad)*math.cos(alt_rad))
    ha_rad = math.acos((math.sin(alt_rad)-math.sin(lat_rad)*math.sin(delt_rad))/math.cos(lat_rad)/math.cos(delt_rad))           
    ha = ha_rad*180/pi
    delt = delt_rad*180/pi
    if math.sin(az_rad)>0:
        ha = 360- ha 
    return (ha,delt)

def enu2ecef(lat,lon,dxyz):
    """
    从 enu(dx,dy,dz) 转换至 ecef(DX,DY,DZ)
    lat: 参考点的纬度，以角度为单位
    lon: 参考点的经度，以角度为单位
    dxyz: [dx, dy, dz]T, enu坐标系下的3xN array, 每一列是相对参考天线的位置.
    """    
#    DXYZ = np.zeros((dxyz.shape[0],dxyz.shape[1]))
    lat = pi/180*lat # in radians
    lon = pi/180*lon # in radians
    L_mat = np.array([[0,-math.sin(lat), math.cos(lat)],
                        [1,0,0],
                        [0,math.cos(lat),math.sin(lat)]])
    LXYZ = np.matmul(L_mat,dxyz)
    return LXYZ


def ecef2uvw(ha,delt,lam,LXYZ):
    """
    从 ECEF(LX,LY,LZ) 计算 uvw(u,v,w)
    ha: 参考点的本地时角, 以角度表示
    delt: 参考点的源赤纬, 以角度表示
    LXYZ: [dx, dy, dz]T, ecef坐标系下的3xN array,每一列是相对参考天线的位置.
    uvw: [u,v,w]T, 3xN array, 每一列是[u,v,w]T空间频率
    lam: 观测波长
    """    
    ha = pi/180*ha # in radians
    delt = pi/180*delt # in radians
    uvw_mat = np.array([math.sin(ha),math.cos(ha),0],
                        [-math.sin(delt)*math.cos(ha),math.sin(delt)*math.sin(ha),math.cos(delt)],
                        [math.cos(delt)*math.cos(ha),-math.cos(delt)*math.sin(ha),math.sin(delt)])/lam
    uvw = np.matmul(uvw_mat,LXYZ)                        
    return uvw  

def wdly_cal(ha,delt,LXYZ):
    """
    From ECEF(LX,LY,LZ) to caculate w
    ha: 参考点的本地时角, 以角度表示
    delt: 参考点的源赤纬, 以角度表示
    LXYZ: [dx, dy, dz]T, 3xN array, ecef坐标系下的3xN array,每一列是相对参考天线的位置.
    w: [u,v,w]的w方向路径长度。
    """      
    ha = pi/180*ha # in radians
    delt = pi/180*delt # in radians
    wdly = np.matmul(np.array([math.cos(delt)*math.cos(ha),-math.cos(delt)*math.sin(ha),math.sin(delt)]),LXYZ)   
    return wdly
