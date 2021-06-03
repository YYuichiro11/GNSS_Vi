#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 19:11:59 2021

@author: yuichiro
"""

import csv


############## GNSSのデータを読み込みます ##############

#一つ目のGNSSファイルを開きます．$GNGGA
csv_file = open(r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524090750_gnss_gngga.csv", "r", encoding="ms932", errors="", newline="" )

#リスト形式
gnss_1 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

header = next(gnss_1)

#データを要素ごとに格納します
time_list = [] # UTC（協定世界時）での時刻．日本標準時は協定世界時より9時間進んでいる．hhmmss.ss
latitude = [] # 緯度．dddmm.mmmm．60分で1度なので，分数を60で割ると度数になります．ddd.dddd度表記は(度数 + 分数/60) で得る．
latitude_direction = [] # 	北緯か南緯か．
longitude = [] # 経度.dddmm.mmmm
longitude_direction = [] # 東経か西経か．
quality_posi = [] # 位置特定品質。0 = 位置特定できない、1 = SPS（標準測位サービス）モード、2 = differenctial GPS（干渉測位方式）モード
satellite = [] # 使用衛星数．
quality_horizontal = [] # 水平精度低下率
antenna_height = [] # アンテナの海抜高さ
geoid = []  #ジオイドの高さ
thread_time = [] # スレッドの書き込み時刻

for row in gnss_1:
    
    time_list.append([float(row[1])])
    latitude.append([float(row[2])])
    latitude_direction.append([row[3]])
    longitude.append([float(row[4])])
    longitude_direction.append([row[5]])
    quality_posi.append([float(row[6])])
    satellite.append([float(row[7])])
    quality_horizontal.append([float(row[8])])
    antenna_height.append([float(row[9])])
    geoid.append([float(row[11])])
    thread_time_cal  = float(row[18])+float(row[19])/60+float(row[20])/3600
    thread_time.append([thread_time_cal])
    

#二つ目のGNSSファイルを開きます．$GNRMC
csv_file = open(r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524090750_gnss_gnrmc.csv", "r", encoding="ms932", errors="", newline="" )

#リスト形式
gnss_2 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

header = next(gnss_2)


velocity = [] # 地表における移動の速度。000.0～999.9[knot]
true_direction = [] # 地表における移動の真方位。000.0～359.9度
for row in gnss_2:
    
    velocity.append([float(row[7])])
    true_direction.append([row[8]])

gnss_all = [time_list,latitude,latitude_direction,longitude,longitude_direction,quality_posi,satellite,quality_horizontal,antenna_height,geoid,thread_time,velocity,true_direction]


############## RealSenseのデータを読み込みます ##############

# １つ目のRealSenseのデータを読み込みます
csv_file = open(r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524090921_t265.csv", "r", encoding="ms932", errors="", newline="" )

#リスト形式
RS_1 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

header = next(RS_1)

camID_1 = []
captureID_1 = []
x_1 = []
y_1 = []
z_1 = []
omg_1 = []
phi_1 = []
kpp_1 = []
thread_time_vo1 = [] # スレッドの書き込み時刻

for row in RS_1:
    
    camID_1.append([float(row[0])])
    captureID_1.append([float(row[1])])
    x_1.append([float(row[2])])
    y_1.append([row[3]])
    z_1.append([float(row[4])])
    omg_1.append([row[5]])
    phi_1.append([float(row[6])])
    kpp_1.append([float(row[7])])
    thread_time_cal  = float(row[11])+float(row[12])/60+float(row[13])/3600
    thread_time_vo1.append([thread_time_cal])
    
    
# ２つ目のRealSenseのデータを読み込みます
csv_file = open(r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524091426_t265.csv", "r", encoding="ms932", errors="", newline="" )

#リスト形式
RS_2 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

header = next(RS_2)

camID_2 = []
captureID_2 = []
x_2 = []
y_2 = []
z_2 = []
omg_2 = []
phi_2 = []
kpp_2 = []
thread_time_vo2 = [] # スレッドの書き込み時刻

for row in RS_2:
    
    camID_2.append([float(row[0])])
    captureID_2.append([float(row[1])])
    x_2.append([float(row[2])])
    y_2.append([row[3]])
    z_2.append([float(row[4])])
    omg_2.append([row[5]])
    phi_2.append([float(row[6])])
    kpp_2.append([float(row[7])])
    thread_time_cal  = float(row[11])+float(row[12])/60+float(row[13])/3600
    thread_time_vo2.append([thread_time_cal])


# ３つ目のRealSenseのデータを読み込みます
csv_file = open(r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524091759_t265.csv", "r", encoding="ms932", errors="", newline="" )

#リスト形式
RS_3 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

header = next(RS_3)

camID_3 = []
captureID_3 = []
x_3 = []
y_3 = []
z_3 = []
omg_3 = []
phi_3 = []
kpp_3 = []
thread_time_vo3 = [] # スレッドの書き込み時刻

for row in RS_3:
    
    camID_3.append([float(row[0])])
    captureID_3.append([float(row[1])])
    x_3.append([float(row[2])])
    y_3.append([row[3]])
    z_3.append([float(row[4])])
    omg_3.append([row[5]])
    phi_3.append([float(row[6])])
    kpp_3.append([float(row[7])])
    thread_time_cal  = float(row[11])+float(row[12])/60+float(row[13])/3600
    thread_time_vo3.append([thread_time_cal])



############## RealSenseのデータをドッキングします ##############

camID = camID_1 + camID_2 + camID_3
captureID = captureID_1 + captureID_2 + captureID_3
x = x_1 + x_2 + x_3
y = y_1 + y_2 + y_3
z = z_1 + z_2 + z_3
omg = omg_1 + omg_2 + omg_3
phi = phi_1 + phi_2 + phi_3
kpp = kpp_1 + kpp_2 + kpp_3
thread_time_vo = thread_time_vo1 + thread_time_vo2 + thread_time_vo3


# (if文の復習，ifでやる必要はない)
# if thread_time_vo1[1] < thread_time_vo1[2]:
#     if thread_time_vo1[2] < thread_time_vo1[3]:
        
#         camID = camID_1 + camID_2 + camID_3
#         captureID = captureID_1 + captureID_2 + captureID_3
#         x = x_1 + x_2 + x_3
#         y = y_1 + y_2 + y_3
#         z = z_1 + z_2 + z_3
#         omg = omg_1 + omg_2 + omg_3
#         phi = phi_1 + phi_2 + phi_3
#         kpp = kpp_1 + kpp_2 + kpp_3
#         thread_time_vo = thread_time_vo1 + thread_time_vo2 + thread_time_vo3
        
        
#     else:
        
#         camID = camID_1 + camID_3 + camID_2
#         captureID = captureID_1 + captureID_3 + captureID_2
#         x = x_1 + x_3 + x_2
#         y = y_1 + y_3 + y_2
#         z = z_1 + z_3 + z_2
#         omg = omg_1 + omg_3 + omg_2
#         phi = phi_1 + phi_3 + phi_2
#         kpp = kpp_1 + kpp_3 + kpp_2
#         thread_time_vo = thread_time_vo1 + thread_time_vo3 + thread_time_vo2

# elif thread_time_vo1[2] < thread_time_vo1[1]:
#     if 








