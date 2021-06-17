#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 19:11:59 2021

@author: yuichiro

Read csv files and store data as list 
Visualize GNSS/RealSense data

"""

import csv
import matplotlib.pyplot as plt
from pyproj import Transformer

############## ファイルのパス ##############
gnss_file1 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524090750_gnss_gngga.csv"
gnss_file2 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524090750_gnss_gnrmc.csv"
rs_file1 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524090921_t265.csv"
rs_file2 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524091426_t265.csv"
rs_file3 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524091759_t265.csv"

############## 関数，クラスを定義します ##############

# GNSSの緯度経度を，度分表記（DMM）から度分秒(DMS)に変換します
def dmm2dms(dmm):
    s = dmm - int(dmm)
    s = s*60
    dms = int(dmm)*100 + s
    return dms

# GNSSの緯度経度を，WGS84系から直交座標系（数学系）へ変換
def wgs2rtg(target, origin):
    result = target-origin
    return result

# RealSenseを読みこむ関数を定義します
class RS_data:
    def __init__(self):
        self.camID = []
        self.captureID = []
        self.x = []
        self.y = []
        self.z = []
        self.omg = []
        self.phi = []
        self.kpp = []
        self.thread_time_vo = [] # スレッドの書き込み時刻
    
    def get_RS(self,RS_list):
        
        header = next(RS_list)
        
        camID = []
        captureID = []
        x = []
        y = []
        z = []
        omg = []
        phi = []
        kpp = []
        thread_time_vo = [] 
        
        for row in RS_list:
            
            camID.append(float(row[0]))
            captureID.append(float(row[1]))
            x.append(float(row[2]))
            y.append(float(row[3]))
            z.append(float(row[4]))
            omg.append(row[5])
            phi.append(float(row[6]))
            kpp.append(float(row[7]))
            thread_time_cal  = float(row[11])+float(row[12])/60+float(row[13])/3600
            thread_time_vo.append(thread_time_cal)
            
        self.camID = camID
        self.captureID = captureID
        self.x = x
        self.y = y
        self.z = z
        self.omg =omg
        self.phi =phi
        self.kpp = kpp
        self.thread_time_vo = thread_time_vo


############## GNSSのデータを読み込みます ##############

#一つ目のGNSSファイルを開きます．$GNGGA
csv_file = open(gnss_file1, "r", encoding="ms932", errors="", newline="" )

#リスト形式
gnss_1 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

header = next(gnss_1)

#データを要素ごとに格納します
time_list = [] # UTC（協定世界時）での時刻．日本標準時は協定世界時より9時間進んでいる．hhmmss.ss
latitude = [] # 緯度．dddmm.mmmm．60分で1度なので，分数を60で割ると度数になります．ddd.dddd度表記は(度数 + 分数/60) で得る．
lat_edit = [] # 直交座標系に変換した値
latitude_direction = [] # 	北緯か南緯か．
longitude = [] # 経度.dddmm.mmmm
lon_edit = [] # 直交座標系に変換した値
longitude_direction = [] # 東経か西経か．
quality_posi = [] # 位置特定品質。0 = 位置特定できない、1 = SPS（標準測位サービス）モード、2 = differenctial GPS（干渉測位方式）モード
satellite = [] # 使用衛星数．
quality_horizontal = [] # 水平精度低下率
antenna_height = [] # アンテナの海抜高さ
geoid = []  #ジオイドの高さ
thread_time = [] # スレッドの書き込み時刻

for row in gnss_1:
    
    time_list.append(float(row[1]))
    latitude_ori = float(row[2])
    latitude_ori = dmm2dms(latitude_ori)
    latitude.append(latitude_ori)
    lat_e = wgs2rtg(latitude_ori, latitude[0])
    lat_edit.append(lat_e)
    latitude_direction.append(row[3])
    longitude_ori = float(row[4])
    longitude_ori = dmm2dms(longitude_ori)
    longitude.append(longitude_ori)
    lon_e = wgs2rtg(longitude_ori, longitude[0])
    lon_edit.append(lon_e)
    longitude_direction.append(row[5])
    quality_posi.append(float(row[6]))
    satellite.append(float(row[7]))
    quality_horizontal.append(float(row[8]))
    antenna_height.append(float(row[9]))
    geoid.append(float(row[11]))
    thread_time_cal  = float(row[18])+float(row[19])/60+float(row[20])/3600
    thread_time.append(thread_time_cal)
    

#二つ目のGNSSファイルを開きます．$GNRMC
csv_file = open(gnss_file2, "r", encoding="ms932", errors="", newline="" )

#リスト形式
gnss_2 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

header = next(gnss_2)


velocity = [] # 地表における移動の速度。000.0～999.9[knot]
true_direction = [] # 地表における移動の真方位。000.0～359.9度
for row in gnss_2:
    
    velocity.append(float(row[7]))
    true_direction.append(row[8])

gnss_all = [time_list,latitude,latitude_direction,longitude,longitude_direction,quality_posi,satellite,quality_horizontal,antenna_height,geoid,thread_time,velocity,true_direction]


############## RealSenseのデータを読み込みます ##############

# １つ目のRealSenseのデータを読み込みます
csv_file = open(rs_file1, "r", encoding="ms932", errors="", newline="" )

#リスト形式
RS_1 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

# データの読み込み
rs1 = RS_data()
rs1.get_RS(RS_1)
    
    
# ２つ目のRealSenseのデータを読み込みます
csv_file = open(rs_file2, "r", encoding="ms932", errors="", newline="" )

#リスト形式
RS_2 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

# データの読み込み
rs2 = RS_data()
rs2.get_RS(RS_2)


# ３つ目のRealSenseのデータを読み込みます
csv_file = open(rs_file3, "r", encoding="ms932", errors="", newline="" )

#リスト形式
RS_3 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

# データの読み込み
rs3 = RS_data()
rs3.get_RS(RS_3)


# RealSenseのデータをドッキングします

camID = rs1.camID + rs2.camID + rs3.camID
captureID = rs1.captureID + rs2.captureID + rs3.captureID
x = rs1.x + rs2.x + rs3.x
y = rs1.y + rs2.y + rs3.y
z = rs1.z + rs2.z + rs3.z
omg = rs1.omg + rs2.omg + rs3.omg
phi = rs1.phi + rs2.phi + rs3.phi
kpp = rs1.kpp + rs2.kpp + rs3.kpp
thread_time_vo = rs1.thread_time_vo + rs2.thread_time_vo + rs3.thread_time_vo


########## GNSSデータの可視化を行います ##########

plt.plot(lon_edit, lat_edit);

plt.plot(rs1.x, rs1.y);
plt.plot(rs2.x, rs2.y);
plt.plot(rs3.x, rs3.y);









