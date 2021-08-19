#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 11:35:56 2021

@author: yuichiro

20210817 
リアルタイム処理をなるべく意識しました．
リアルタイムに座標を表示していくと処理が重ずぎるので，可視化の処理のみ最後にまとめて出力するようにしました
回転と累積誤差補正をしました
20210819
gnrmcフォーマットで出力します，出力は逐次的にすると重くなりそうなので．最後にまとめて行います
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import time

# 計測スタート
t0 = time.time() 

""" 
GNSSの緯度経度を，度分表記（DMM）から度(DDD)に変換します 
"""
def dmm2ddd(dmm):
    d = dmm*0.01
    m = (d - int(d))*100
    ddd = int(d) + m/60
    return ddd
""" 
GNSSの緯度経度を，度(DDD)から度分表記（DMM）に変換します 
"""
def ddd2dmm(ddd):
    d = int(ddd)
    m = (ddd - d)*60
    dmm = d*100 + m
    return dmm

""" 
GNSSの緯度経度を，平面直角座標系へ変換します 
平面直角座標系9系
'緯度経度と平面直角座標の相互変換をPythonで実装する'
https://qiita.com/sw1227/items/e7a590994ad7dcd0e8ab
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


def calc_lat_lon(x, y, phi0_deg, lambda0_deg):
    """ 平面直角座標を緯度経度に変換する
    - input:
        (x, y): 変換したいx, y座標[m]
        (phi0_deg, lambda0_deg): 平面直角座標系原点の緯度・経度[度]（分・秒でなく小数であることに注意）
    - output:
        latitude:  緯度[度]
        longitude: 経度[度]
        * 小数点以下は分・秒ではないことに注意
    """
    # 平面直角座標系原点をラジアンに直す
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

    def beta_array(n):
        b0 = np.nan # dummy
        b1 = (1./2)*n - (2./3)*(n**2) + (37./96)*(n**3) - (1./360)*(n**4) - (81./512)*(n**5)
        b2 = (1./48)*(n**2) + (1./15)*(n**3) - (437./1440)*(n**4) + (46./105)*(n**5)
        b3 = (17./480)*(n**3) - (37./840)*(n**4) - (209./4480)*(n**5)
        b4 = (4397./161280)*(n**4) - (11./504)*(n**5)
        b5 = (4583./161280)*(n**5)
        return np.array([b0, b1, b2, b3, b4, b5])

    def delta_array(n):
        d0 = np.nan # dummy
        d1 = 2.*n - (2./3)*(n**2) - 2.*(n**3) + (116./45)*(n**4) + (26./45)*(n**5) - (2854./675)*(n**6)
        d2 = (7./3)*(n**2) - (8./5)*(n**3) - (227./45)*(n**4) + (2704./315)*(n**5) + (2323./945)*(n**6)
        d3 = (56./15)*(n**3) - (136./35)*(n**4) - (1262./105)*(n**5) + (73814./2835)*(n**6)
        d4 = (4279./630)*(n**4) - (332./35)*(n**5) - (399572./14175)*(n**6)
        d5 = (4174./315)*(n**5) - (144838./6237)*(n**6)
        d6 = (601676./22275)*(n**6)
        return np.array([d0, d1, d2, d3, d4, d5, d6])

    # 定数 (a, F: 世界測地系-測地基準系1980（GRS80）楕円体)
    m0 = 0.9999 
    a = 6378137.
    F = 298.257222101

    # (1) n, A_i, beta_i, delta_iの計算
    n = 1. / (2*F - 1)
    A_array = A_array(n)
    beta_array = beta_array(n)
    delta_array = delta_array(n)

    # (2), S, Aの計算
    A_ = ( (m0*a)/(1.+n) )*A_array[0]
    S_ = ( (m0*a)/(1.+n) )*( A_array[0]*phi0_rad + np.dot(A_array[1:], np.sin(2*phi0_rad*np.arange(1,6))) )

    # (3) xi, etaの計算
    xi = (x + S_) / A_
    eta = y / A_

    # (4) xi', eta'の計算
    xi2 = xi - np.sum(np.multiply(beta_array[1:], 
                                  np.multiply(np.sin(2*xi*np.arange(1,6)),
                                              np.cosh(2*eta*np.arange(1,6)))))
    eta2 = eta - np.sum(np.multiply(beta_array[1:],
                                   np.multiply(np.cos(2*xi*np.arange(1,6)),
                                               np.sinh(2*eta*np.arange(1,6)))))

    # (5) chiの計算
    chi = np.arcsin( np.sin(xi2)/np.cosh(eta2) ) # [rad]
    latitude = chi + np.dot(delta_array[1:], np.sin(2*chi*np.arange(1, 7))) # [rad]

    # (6) 緯度(latitude), 経度(longitude)の計算
    longitude = lambda0_rad + np.arctan( np.sinh(eta2)/np.cos(xi2) ) # [rad]

    # ラジアンを度になおしてreturn
    return np.rad2deg(latitude), np.rad2deg(longitude) # [deg]


# windows環境
# gnss_file1 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524090750_gnss_gngga.csv"
# gnss_file2 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524090750_gnss_gnrmc.csv"
# rs_file1 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524090921_t265.csv"
# rs_file2 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524091426_t265.csv"
# rs_file3 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524091759_t265.csv"
# rs_file4 = r"C:\Users\jamis\Desktop\research\GNSS_visual\testdata_GNSS_VisualOdometry\dstr_20210524091949_t265.csv"

# mac環境
gnss_file1 = r"/Users/yuichiro/Research/code/GNSS_visual/20210529_benchmark_dataset/testdata_GNSS_VisualOdometry/dstr_20210524090750_gnss_gngga.csv"
gnss_file2 = r"/Users/yuichiro/Research/code/GNSS_visual/20210529_benchmark_dataset/testdata_GNSS_VisualOdometry/dstr_20210524090750_gnss_gnrmc.csv"
rs_file1 = r"/Users/yuichiro/Research/code/GNSS_visual/20210529_benchmark_dataset/testdata_GNSS_VisualOdometry/dstr_20210524090921_t265.csv"
rs_file2 = r"/Users/yuichiro/Research/code/GNSS_visual/20210529_benchmark_dataset/testdata_GNSS_VisualOdometry/dstr_20210524091426_t265.csv"
rs_file3 = r"/Users/yuichiro/Research/code/GNSS_visual/20210529_benchmark_dataset/testdata_GNSS_VisualOdometry/dstr_20210524091759_t265.csv"
rs_file4 = r"/Users/yuichiro/Research/code/GNSS_visual/20210529_benchmark_dataset/testdata_GNSS_VisualOdometry/dstr_20210524091949_t265.csv"

"""データの読み込み"""
csv_file = open(gnss_file1, "r", encoding="ms932", errors="", newline="" )
gnss_1 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
g_1 = [row for row in gnss_1]

csv_file = open(gnss_file2, "r", encoding="ms932", errors="", newline="" )
gnss_2 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
g_2 = [row for row in gnss_2]

csv_file = open(rs_file1, "r", encoding="ms932", errors="", newline="" )
rs1 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
r_1 = [row for row in rs1]

"""
データをひとつづつ読み込んで，仮想リアルタイム
リアルタイムにデータを読むやり方がよくわからないから，gnssとrealsenseのデータをひとつづつ読み込んでその都度処理を行うことで
未来のデータはわからない状態での処理を想定する
"""

"""
ここから先の処理は，リアルタイムのデータ受け取りを想定して，
まずは同じ時刻のデータを同じタイミングで受け取れるようにします（あとで関数化する）
"""
row_g = 1;
row_r = 1;
# まずそれぞれのデータにおける時間を引っ張り出す
g_time = float(g_1[row_g][18])+float(g_1[row_g][19])/60+float(g_1[row_g][20])/3600
rs_time = float(r_1[row_r][11])+float(r_1[row_r][12])/60+float(r_1[row_r][13])/3600

# 初期位置の位置合わせを行い，並進ベクトルを確定します
# gnssのデータの方が時刻が昔のデータである場合
if g_time - rs_time < 0:
    while g_time - rs_time < 0:
        row_g += 1
        g_time = float(g_1[row_g][18])+float(g_1[row_g][19])/60+float(g_1[row_g][20])/3600
    # realsenseのデータの方がデータが細かいため．realsenseのデータを使って時間調整
    while rs_time - g_time < 0:
        row_r += 1
        rs_time = float(r_1[row_r][11])+float(r_1[row_r][12])/60+float(r_1[row_r][13])/3600
   
    # ralsenseのデータの方が時刻が昔のデータである場合
else:
    while rs_time - g_time < 0:
        row_r += 1
        rs_time = float(r_1[row_r][11])+float(r_1[row_r][12])/60+float(r_1[row_r][13])/3600
            
# 初期位置合わせを行います
# gnssの緯度経度データを平面直角座標９系に変換
g_lat_0 = dmm2ddd(float(g_1[row_g][2]))
g_lon_0 = dmm2ddd(float(g_1[row_g][4]))
g_lat_0, g_lon_0 = calc_xy(g_lat_0, g_lon_0, 36., 139+50./60)
# realsenseのデータを読み込み
r_x_0 = float(r_1[row_r][2])
r_y_0 = float(r_1[row_r][3])
# 並進ベクトルtの初期値を求めます
t_x_0 = g_lon_0 - r_x_0
t_y_0 = g_lat_0 - r_y_0

# ここがrealsenseの初期地点
r_x_0 = r_x_0 + t_x_0
r_y_0 = r_y_0 + t_y_0

# 回転行列の初期値をセットします
cosR = 1
sinR = 0


"""
ここから実際にそれぞれのデータをひとつづつ読み込み，並進・回転・累積誤差補正を行います
リアルタイム処理を意識するために，時間の流れを意識した処理に寄せます
"""
# データをまとめる変数を用意します（今回はとりあえずappendで格納，本来は予め1万とか10万とかのlistを作成して，その中に格納する方が効率良い）
x_rs = []
y_rs = []

x_g = []
y_g = []

rs_gnrmcs = []

# gnssの現在地を読み込んでおきます
g_lat = dmm2ddd(float(g_1[row_g][2]))
g_lon = dmm2ddd(float(g_1[row_g][4]))
g_lat, g_lon = calc_xy(g_lat, g_lon, 36., 139+50./60)

# 角度の初期値を指定します
g_dig = ''
r_kpp = ''

# fig = plt.figure()

while True:
    g_time = float(g_1[row_g][18])+float(g_1[row_g][19])/60+float(g_1[row_g][20])/3600
    # gmssのデータ取得は1秒おきだから，次に取得できるのは今取得したデータの1秒後，これを使って，1秒後までrsのデータを読み込む処理をします
    next_g = float(g_1[row_g][18])+float(g_1[row_g][19])/60+(float(g_1[row_g][20])+1)/3600
    rs_time = float(r_1[row_r][11])+float(r_1[row_r][12])/60+float(r_1[row_r][13])/3600
    
    
    # gnssの地表に対する移動の真方位を読み込みます
    if g_2[row_g][8] != '':
        g_dig = float(g_2[row_g][8])
    # kppの値を読み込みます
    if r_1[row_r][7] != '':
        r_kpp = float(r_1[row_r][7])
    
    # gnssとrealsenseの方位のずれを計算します
    if g_dig != '' and r_kpp != '':
        if r_kpp < g_dig:
            theta = 360 + (r_kpp - g_dig)
        else:
            theta = r_kpp - g_dig
        
        # 回転行列を確定します
        cosR = np.cos(np.radians(theta))
        sinR = np.sin(np.radians(theta))
    
    
    # realsenseの現在地に，最初に計算した並進ベクトルを足し合わせます
    r_x = float(r_1[row_r][2]) + t_x_0
    r_y = float(r_1[row_r][3]) + t_y_0
    
    rs_vec = [r_x - r_x_0, r_y - r_y_0]
    rs_vector = [cosR*rs_vec[0] - sinR*rs_vec[1], sinR*rs_vec[0] + cosR*rs_vec[1]]
    r_x = rs_vector[0]+r_x_0
    r_y = rs_vector[1]+r_y_0
    
    
    # 次の更新タイミングへ
    row_r += 1
    
    # 1秒経過したことを疑似的に再現
    if rs_time > next_g:
        # 1秒経過したから，gnssのデータを更新
        row_g += 1
        g_lat = dmm2ddd(float(g_1[row_g][2]))
        g_lon = dmm2ddd(float(g_1[row_g][4]))
        g_lat, g_lon = calc_xy(g_lat, g_lon, 36., 139+50./60)
        x_g.append(g_lon)
        y_g.append(g_lat)
        
        # 累積誤差補正
        t_x_0 = g_lon - float(r_1[row_r][2])
        t_y_0 = g_lat - float(r_1[row_r][3])
        r_x_0 = float(r_1[row_r][2]) + t_x_0
        r_y_0 = float(r_1[row_r][3]) + t_y_0        
        
        
    
    x_rs.append(r_x)
    y_rs.append(r_y)
    
    # realsenseのx,yの結果を緯度軽度へ変換します
    r_lat, r_lon = calc_lat_lon(r_x, r_y, 33., 131.)
    r_lat = ddd2dmm(r_lat)
    r_lon = ddd2dmm(r_lon)
    
    # 7から12までgnssと同じ
    rs_gnrmc = ['$GNRMC', str(g_2[row_g][1]), str(g_2[row_g][2]), str(r_lat), str(g_2[row_g][4]), str(r_lon), str(g_2[row_g][6]),  str(g_2[row_g][7]), str(g_dig),  str(g_2[row_g][9]),  str(g_2[row_g][10]),  str(g_2[row_g][11]),  str(g_2[row_g][12]), r_1[row_r][8], r_1[row_r][9], r_1[row_r][10], r_1[row_r][11], r_1[row_r][12], r_1[row_r][13]]
    rs_gnrmcs.append(rs_gnrmc)
    
    
    # データ取得がうまくいかなくなったら一回処理を終了
    if g_1[row_g][2] == 'nan':
        break
    
    if r_1[row_r][2] == 'nan':
        break

# 可視化
plt.plot(x_rs, y_rs)
plt.plot(x_g, y_g)

# csvファイルへ出力
p_new = pathlib.Path('/Users/yuichiro/Research/code/GNSS_visual/rs_out.csv')
with p_new.open(mode='w') as f:
  f.write('')

header = g_2[0]
with open('/Users/yuichiro/Research/code/GNSS_visual/rs_out.csv', 'w', encoding='cp932') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(rs_gnrmcs)

t1 = time.time()                # 計測終了時間
elapsed_time = float(t1 - t0)   # 経過時間
print(elapsed_time)             # 経過時間を表示
