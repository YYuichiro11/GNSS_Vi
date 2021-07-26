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
""" 
GNSSの緯度経度を，度分表記（DMM）から度(DDD)に変換します 
"""
def dmm2ddd(dmm):
    d = dmm*0.01
    m = (d - int(d))*100
    ddd = int(d) + m/60
    return ddd

""" 
GNSSの緯度経度を，平面直角座標系へ変換します 
"""
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


""" 
Real SenseのデータとGNSSのデータを使って，初期位置合わせを行います 
"""
def RS_GNSSposi(RS, gnss):
    rs_time = RS.thread_time_vo
    gnss_time = gnss.thread_time
        
    # gnssの時間とrealsenseの時間合わせを行います
    count = 0
    while gnss_time[count] - rs_time[0] < 0:
        count += 1
    else:
        if abs(gnss_time[count] - rs_time[0]) < abs(gnss_time[count-1] - rs_time[0]):
            ori_x = gnss.lon_edit[count]
            ori_y = gnss.lat_edit[count]
            
            # ここまでで時刻合わせが完了したから，gnssの初期値をrealsenseのデータに足し合わせる（並進）
            k = 0
            for rst in rs_time:
                RS.x[k] = RS.x[k] + ori_x
                RS.y[k] = RS.y[k] + ori_y
                k += 1
            
        else:
            ori_x = gnss.lon_edit[count-1]
            ori_y = gnss.lat_edit[count-1]
        
            # ここまでで時刻合わせが完了したから，gnssの初期値をrealsenseのデータに足し合わせる（並進）
            k = 0
            for rst in rs_time:
                RS.x[k] = RS.x[k] + ori_x
                RS.y[k] = RS.y[k] + ori_y
                k += 1
    

"""
Real SenseのデータとGNSSのデータを使って，realsenseの位置データに回転をかけて，正しい方向に直します
原点と二点目のgnssデータを使ってベクトルを求めて，真値とrealsenseのデータのすり合わせを行う二次元の回転行列を求めます 
"""


"""
Real SenseのデータとGNSSのデータの時刻合わせをします
"""
def RS_GNSSmatch(RS, gnss):
    gnss_time = gnss.thread_time
    rs_time = RS.thread_time_vo
    
    # gnssのデータを基軸に，rsとgnssの時刻合わせをします
    for gt in gnss_time:
        count = 0
        
        # gnssの時間がrealsenseの時間に近づくまで処理をスキップするためのif文
        if gt - rs_time[0] > 0:    
            
            # Real Senseの時間を合わせます
            for rt in rs_time:
                cal = abs(gt - rt)  # gnssとrealsenseの時刻の差を計算して，差の絶対値が最も小さいものが同時刻として，位置合わせを行う
                
                if count == 0:
                    minimum = cal
                    
                elif minimum > cal:
                    minimum = cal
                
                else:   # ここでgnssとrealsenseの時刻が同期されているから，この時のcountを用いてデータ合わせを行う
                    same_time = count
                    # ここに位置とかのすり合わせ処理を組み込む
                    # この関数において，realsenseのデータが最後まで読み込まれたにもかかわらず，時刻同期ができないときに処理を終わらせるようにする
                    # 次のrealsenseのデータを読み込むための，エラー回避用の処理を組む
                    
                    break
                
                count += 1
            

            
        
    
    

"""
 GNSSデータを読み込むクラスを定義します
"""
class GNSS_data:
    def __init__(self):
        #データを要素ごとに格納します
        self.time_list = [] # UTC（協定世界時）での時刻．日本標準時は協定世界時より9時間進んでいる．hhmmss.ss
        self.latitude = [] # 緯度．dddmm.mmmm．60分で1度なので，分数を60で割ると度数になります．ddd.dddd度表記は(度数 + 分数/60) で得る．
        self.latitude_direction = [] # 	北緯か南緯か．
        self.longitude = [] # 経度.dddmm.mmmm
        self.longitude_direction = [] # 東経か西経か．
        self.quality_posi = [] # 位置特定品質。0 = 位置特定できない、1 = SPS（標準測位サービス）モード、2 = differenctial GPS（干渉測位方式）モード
        self.satellite = [] # 使用衛星数．
        self.quality_horizontal = [] # 水平精度低下率
        self.antenna_height = [] # アンテナの海抜高さ
        self.geoid = []  #ジオイドの高さ
        self.thread_time = [] # スレッドの書き込み時刻
        self.lat_edit = [] # 直交座標系に変換した値
        self.lon_edit = [] # 直交座標系に変換した値
        self.velocity = [] # 地表における移動の速度。000.0～999.9[knot]
        self.true_direction = [] # 地表における移動の真方位。000.0～359.9度
    
    def get_GNGGA(self,GNGGA_list):
        header = next(GNGGA_list)
        for row in GNGGA_list:
            self.time_list.append(float(row[1]))
            latitude_ori = float(row[2])
            self.latitude.append(latitude_ori)
            
            self.latitude_direction.append(row[3])
            longitude_ori = float(row[4])
            self.longitude.append(longitude_ori)
        
            self.longitude_direction.append(row[5])
            self.quality_posi.append(float(row[6]))
            self.satellite.append(float(row[7]))
            self.quality_horizontal.append(float(row[8]))
            self.antenna_height.append(float(row[9]))
            self.geoid.append(float(row[11]))
            thread_time_cal  = float(row[18])+float(row[19])/60+float(row[20])/3600
            self.thread_time.append(thread_time_cal)
            
            latitude_ori = dmm2ddd(latitude_ori)
            longitude_ori = dmm2ddd(longitude_ori)
            lat_e, lon_e = calc_xy(latitude_ori, longitude_ori, 36., 139+50./60)
            self.lat_edit.append(lat_e)
            self.lon_edit.append(lon_e)
    
    def get_GNRMC(self,GNRMC_list):
        header = next(GNRMC_list)
        for row in GNRMC_list:
            self.velocity.append(float(row[7]))
            self.true_direction.append(row[8])
            

""" 
RealSenseを読みこむクラスを定義します
"""
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


"""############## GNSSのデータを読み込みます ##############"""
#一つ目のGNSSファイルを開きます．$GNGGA
csv_file = open(gnss_file1, "r", encoding="ms932", errors="", newline="" )
#リスト形式
gnss_1 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
# データの読み込み
gnss = GNSS_data()
gnss.get_GNGGA(gnss_1)

#二つ目のGNSSファイルを開きます．$GNRMC
csv_file = open(gnss_file2, "r", encoding="ms932", errors="", newline="" )
#リスト形式
gnss_2 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
gnss.get_GNRMC(gnss_2)

"""############## RealSenseのデータを読み込みます ##############"""
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


"""########## GNSSデータの可視化を行います ##########"""

plt.plot(gnss.lon_edit, gnss.lat_edit);
# plt.plot(rs1.x, rs1.y);
# plt.plot(rs2.x, rs2.y);
# plt.plot(rs3.x, rs3.y);
# plt.plot(rs4.x, rs4.y);

# 初期位置合わせを行います
RS_GNSSposi(rs1, gnss)

plt.plot(gnss.lon_edit, gnss.lat_edit);
plt.plot(rs1.x, rs1.y);

RS_GNSSposi(rs2, gnss)
plt.plot(rs2.x, rs2.y);

RS_GNSSposi(rs3, gnss)
plt.plot(rs3.x, rs3.y);

RS_GNSSposi(rs4, gnss)
plt.plot(rs4.x, rs4.y);

