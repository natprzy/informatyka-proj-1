# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:06:05 2022

@author: admin
"""


import numpy as np
import kody_gw as k

el_grs84 = k.Transformacje(model = "wgs84")
plik = "wsp_int.txt"

tablica = np.genfromtxt(plik, delimiter = ",", skip_header = 4)
rows,cols = np.shape(tablica)

blh = np.zeros((rows,cols))
wsp_2000 = np.zeros((rows,3))
wsp_1992 = np.zeros((rows,2))
neu = np.zeros((rows,cols))


tablica_wsp = np.zeros((rows,8))

for i in range(rows):
    blh[i] = el_grs84.hirvonen(tablica[i,0], tablica[i,1], tablica[i,2])
    wsp_2000[i] = el_grs84.u2000(blh[i,0], blh[i,1], 0.999923, 21)
    wsp_1992[i] = el_grs84.u92(blh[i,0], blh[i,1], 0.9993)
    
    #neu[i] = el_grs84.neu(tablica[i,0], tablica[i,1], tablica[i,2], tablica[i,0]+1, tablica[i,1]+1)
    tablica_wsp[i,0:3] = blh[i]
    tablica_wsp[i,3:6] = wsp_2000[i]
    tablica_wsp[i,6:8] = wsp_1992[i]
    
np.savetxt("wsp_BLH.txt", tablica_wsp, delimiter=',')