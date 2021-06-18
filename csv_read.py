#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 19:11:59 2021

@author: yuichiro

Read csv files and store data as list 
Visualize GNSS/RealSense data
Convert GNSS data from wgs84 to rect jp9

"""

import csv
import matplotlib.pyplot as plt
import numpy as np

############## ファイルのパス ##############
gnss_file1 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524090750_gnss_gngga.csv"
gnss_file2 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524090750_gnss_gnrmc.csv"
rs_file1 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524090921_t265.csv"
rs_file2 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524091426_t265.csv"
rs_file3 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524091759_t265.csv"
rs_file4 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524091949_t265.csv"

############## 関数，クラスを定義します ##############
# GNSSの緯度経度を，度分表記（DMM）から度(DDD)に変換します
def dmm2ddd(dmm):
    d = dmm*0.01
    m = (d - int(d))*100
    ddd = int(d) + m/60
    return ddd

# GNSSの緯度経度を，平面直角座標系へ変換します
def calc_xy(phi_deg, lambda_deg, phi0_deg, lambda0_deg):
    """ 緯度経度を平面直角座標に変換する
    - input:
        (phi_deg, lambda_deg): 変換したい緯度・経度[度]（分・秒でなく小数であることに注意）
        (phi0_deg, lambda0_deg): 平面直角座標系原点の緯度・経度[度]（分・秒でなく小数であることに注意）
    - output:
        x: 変換後の平面直角座標[m]
        y: 変換後の平面直角座標[m]
    """
    # 緯度経度・平面直角座標系原点をラジアンに直す
    phi_rad = np.deg2rad(phi_deg)
    lambda_rad = np.deg2rad(lambda_deg)
    phi0_rad = np.deg2rad(phi0_deg)
    lambda0_rad = np.deg2rad(lambda0_deg)

    # 補助関数
    def A_array(n):
        A0 = 1 + (n**2)/4. + (n**4)/64.
        A1 = -     (3./2)*( n - (n**3)/8. - (n**5)/64. ) 
        A2 =     (15./16)*( n**2 - (n**4)/4. )
        A3 = -   (35./48)*( n**3 - (5./16)*(n**5) )
        A4 =   (315./512)*( n**4 )
        A5 = -(693./1280)*( n**5 )
        return np.array([A0, A1, A2, A3, A4, A5])

    def alpha_array(n):
        a0 = np.nan # dummy
        a1 = (1./2)*n - (2./3)*(n**2) + (5./16)*(n**3) + (41./180)*(n**4) - (127./288)*(n**5)
        a2 = (13./48)*(n**2) - (3./5)*(n**3) + (557./1440)*(n**4) + (281./630)*(n**5)
        a3 = (61./240)*(n**3) - (103./140)*(n**4) + (15061./26880)*(n**5)
        a4 = (49561./161280)*(n**4) - (179./168)*(n**5)
        a5 = (34729./80640)*(n**5)
        return np.array([a0, a1, a2, a3, a4, a5])

    # 定数 (a, F: 世界測地系-測地基準系1980（GRS80）楕円体)
    m0 = 0.9999 
    a = 6378137.
    F = 298.257222101

    # (1) n, A_i, alpha_iの計算
    n = 1. / (2*F - 1)
    A_array = A_array(n)
    alpha_array = alpha_array(n)

    # (2), S, Aの計算
    A_ = ( (m0*a)/(1.+n) )*A_array[0] # [m]
    S_ = ( (m0*a)/(1.+n) )*( A_array[0]*phi0_rad + np.dot(A_array[1:], np.sin(2*phi0_rad*np.arange(1,6))) ) # [m]

    # (3) lambda_c, lambda_sの計算
    lambda_c = np.cos(lambda_rad - lambda0_rad)
    lambda_s = np.sin(lambda_rad - lambda0_rad)

    # (4) t, t_の計算
    t = np.sinh( np.arctanh(np.sin(phi_rad)) - ((2*np.sqrt(n)) / (1+n))*np.arctanh(((2*np.sqrt(n)) / (1+n)) * np.sin(phi_rad)) )
    t_ = np.sqrt(1 + t*t)

    # (5) xi', eta'の計算
    xi2  = np.arctan(t / lambda_c) # [rad]
    eta2 = np.arctanh(lambda_s / t_)

    # (6) x, yの計算
    x = A_ * (xi2 + np.sum(np.multiply(alpha_array[1:],
                                       np.multiply(np.sin(2*xi2*np.arange(1,6)),
                                                   np.cosh(2*eta2*np.arange(1,6)))))) - S_ # [m]
    y = A_ * (eta2 + np.sum(np.multiply(alpha_array[1:],
                                        np.multiply(np.cos(2*xi2*np.arange(1,6)),
                                                    np.sinh(2*eta2*np.arange(1,6)))))) # [m]
    # return
    return x, y # [m]


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
        
        for row in RS_list:
            
            self.camID.append(float(row[0]))
            self.captureID.append(float(row[1]))
            self.x.append(float(row[2]))
            self.y.append(float(row[3]))
            self.z.append(float(row[4]))
            self.omg.append(row[5])
            self.phi.append(float(row[6]))
            self.kpp.append(float(row[7]))
            thread_time_cal  = float(row[11])+float(row[12])/60+float(row[13])/3600
            self.thread_time_vo.append(thread_time_cal)

############## GNSSのデータを読み込みます ##############
#一つ目のGNSSファイルを開きます．$GNGGA
csv_file = open(gnss_file1, "r", encoding="ms932", errors="", newline="" )
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
lat_edit = [] # 直交座標系に変換した値
lon_edit = [] # 直交座標系に変換した値

for row in gnss_1:
    time_list.append(float(row[1]))
    latitude_ori = float(row[2])
    latitude.append(latitude_ori)
        
    latitude_direction.append(row[3])
    longitude_ori = float(row[4])
    longitude.append(longitude_ori)
        
    longitude_direction.append(row[5])
    quality_posi.append(float(row[6]))
    satellite.append(float(row[7]))
    quality_horizontal.append(float(row[8]))
    antenna_height.append(float(row[9]))
    geoid.append(float(row[11]))
    thread_time_cal  = float(row[18])+float(row[19])/60+float(row[20])/3600
    thread_time.append(thread_time_cal)
    
    latitude_ori = dmm2ddd(latitude_ori)
    longitude_ori = dmm2ddd(longitude_ori)
    lat_e, lon_e = calc_xy(latitude_ori, longitude_ori, 36., 139+50./60)
    lat_edit.append(lat_e)
    lon_edit.append(lon_e)
    
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

# ４つ目のRealSenseのデータを読み込みます
csv_file = open(rs_file4, "r", encoding="ms932", errors="", newline="" )
#リスト形式
RS_4 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
# データの読み込み
rs4 = RS_data()
rs4.get_RS(RS_4)


# RealSenseのデータをドッキングします
# camID = rs1.camID + rs2.camID + rs3.camID + rs4.camID
# captureID = rs1.captureID + rs2.captureID + rs3.captureID + rs4.captureID
# x = rs1.x + rs2.x + rs3.x +rs4.x
# y = rs1.y + rs2.y + rs3.y + rs4.y
# z = rs1.z + rs2.z + rs3.z + rs4.z
# omg = rs1.omg + rs2.omg + rs3.omg + rs4.omg 
# phi = rs1.phi + rs2.phi + rs3.phi + rs4.phi
# kpp = rs1.kpp + rs2.kpp + rs3.kpp + rs4.kpp
# thread_time_vo = rs1.thread_time_vo + rs2.thread_time_vo + rs3.thread_time_vo + rs4.thread_time_vo


########## GNSSデータの可視化を行います ##########

# plt.plot(lon_edit, lat_edit);

plt.plot(rs1.x, rs1.y);
# plt.plot(rs2.x, rs2.y);
# plt.plot(rs3.x, rs3.y);
# plt.plot(rs4.x, rs4.y);









