# csv2kml.py
#-*- coding: utf-8 -*-
#https://qiita.com/tomo001/items/b375e5fa578eb8880662
# usage csv2kml.py point_file.csv out.kml

# example of a csv file :
#ID1, 139.766389, 35.681340
#ID2, 139.763360, 35.675056
# ...
#ID10, 139.777195, 35.713714

import sys
import csv
import simplekml


n_argv=len(sys.argv)
if (n_argv != 3):
    print ("Incorrect number of arguments.\nUsage: csv2kml.py coordinate.csv output.kml")
    exit(1)

in_file = sys.argv[1]
out_file = sys.argv[2]

with open(in_file) as f:
    reader = csv.reader(f)
    sample_points = [row for row in reader]


kml = simplekml.Kml()
for point in sample_points:
    kml.newpoint(name=point[0], coords=[(point[1], point[2])])

kml.save(out_file)
