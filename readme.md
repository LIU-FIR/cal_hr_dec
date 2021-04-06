# readme

用于计算MUSER阵列天线的指向的c++程序：以时角-赤纬表示

cal_hr_dec.cpp包含了所需的函数sun_hr_dec()，利用结构体返回hour_angle和declination.

cal_jd()计算sun_hr_dec()所需要的JD。

dT根据不同年份给出。

main()测试cal_hr_dec的计算值。