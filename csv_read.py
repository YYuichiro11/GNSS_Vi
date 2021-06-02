#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 19:11:59 2021

@author: yuichiro
"""

import csv

csv_file = open("/Users/yuichiro/Research/20210529_benchmark_dataset/testdata_GNSS_VisualOdometry/dstr_20210524090750_gnss_gngga.csv", "r", encoding="ms932", errors="", newline="" )

#リスト形式
f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

header = next(f)
print(header)

a = next(f)
print(a)


time_list = []
latitude = []
latitude_direction = []
longitude = []
longitude_direction = []
quality_posi = []
satellite = []
quality_horizontal = []
antenna_height = []
geoid = []
thread_time = []


for row in f:
    
    time_list.append([row[1]])
    latitude.append([row[2]])
    latitude_direction.append([row[3]])
    longitude.append([row[4]])
    longitude_direction.append([row[5]])
    quality_posi.append([row[6]])
    satellite.append([row[7]])
    quality_horizontal.append([row[8]])
    antenna_height.append([row[9]])
    geoid.append([row[11]])
    thread_time_cal  = float(row[18])+float(row[19])/60+float(row[20])/3600
    thread_time.append([thread_time_cal])
    
