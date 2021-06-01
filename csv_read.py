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

#辞書形式
g = csv.DictReader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)



