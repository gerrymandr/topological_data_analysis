# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 11:16:29 2018

@author: tug30201
"""

import os

# geospatial
import geopandas as gpd
import pysal as ps
import numpy as np
from tqdm import tqdm

# visualization
import matplotlib.pyplot as plt

# data retrieval
import os
from urllib import urlopen
from zipfile import ZipFile
from tqdm import tqdm

url = "https://www2.census.gov/geo/tiger/TIGER2012/VTD/"
name = "tl_2012_FIPS_vtd10.zip"

#this is missing KY and RI 

FIPS = {"02":"AK", "01":"AL", "05":"AR", "04":"AZ", "06":"CA",
        "08":"CO", "09":"CT", "10":"DE", "12":"FL", "13":"GA",
        "15":"HI", "19":"IA", "16":"ID", "17":"IL", "18":"IN",
        "20":"KS", "22":"LA", "25":"MA", "24":"MD",
        "23":"ME", "26":"MI", "27":"MN", "29":"MO", "28":"MS",
        "30":"MT", "37":"NC", "38":"ND", "31":"NE", "33":"NH",
        "34":"NJ", "35":"NM", "32":"NV", "36":"NY", "39":"OH",
        "40":"OK", "41":"OR", "42":"PA", "45":"SC",
        "46":"SD", "47":"TN", "48":"TX", "49":"UT", "51":"VA",
        "50":"VT", "53":"WA", "55":"WI", "54":"WV", "56":"WY"}

# data retrieval
def get_and_unzip(url, data_dir=os.getcwd()):
    basename = url.split("/")[-1]
    if not os.path.exists(os.path.join(data_dir, basename)):
        file_data = urlopen(url)
        data_to_write = file_data.read()
        with open(basename, "wb") as f:
            f.write(data_to_write)

        zip_obj = ZipFile(basename)
        zip_obj.extractall()
        del(zip_obj)


def get_centroids(fname, state):
    df_counties = gpd.read_file(fname)
    
    county_centroids = df_counties.centroid
    c_x = county_centroids.x
    c_y = county_centroids.y
        
    with open("%s.txt" % state, "w") as f:
        for z in zip(c_x, c_y):
            f.write("%.3f %.3f\n" % z)

for k in tqdm(FIPS.keys()):
    fname = name.replace("FIPS", k)
    get_and_unzip(url+fname, ".")
    get_centroids(fname.replace("zip", "shp"), FIPS[k])

