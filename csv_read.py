#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 19:11:59 2021

@author: yuichiro
"""

import csv

#一つ目のGNSSファイルを開きます．$GNGGA
csv_file = open("/Users/yuichiro/Research/20210529_benchmark_dataset/testdata_GNSS_VisualOdometry/dstr_20210524090750_gnss_gngga.csv", "r", encoding="ms932", errors="", newline="" )

#リスト形式
f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

header = next(f)
print(header)

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

for row in f:
    
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
csv_file = open("/Users/yuichiro/Research/20210529_benchmark_dataset/testdata_GNSS_VisualOdometry/dstr_20210524090750_gnss_gnrmc.csv", "r", encoding="ms932", errors="", newline="" )

#リスト形式
f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

header = next(f)
print(header)


velocity = [] # 地表における移動の速度。000.0～999.9[knot]
true_direction = [] # 地表における移動の真方位。000.0～359.9度
for row in f:
    
    velocity.append([float(row[7])])
    true_direction.append([row[8]])



