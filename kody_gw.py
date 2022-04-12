# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 09:46:42 2022

@author: Natalia
"""
import math  
import numpy as np
import statistics as stt
class Transformacje:
    def __init__(self, model: str = "wgs84"):
        #https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.ecc2 = 2 * self.flattening - self.flattening ** 2




        def hirvonen(self,x, y, z, a, e2):
            eps1=0.000001 #sekundy
            eps=eps1/3600*math.pi/180 #rad
            
            
            r=math.sqrt(x**2 + y**2)
            
            
            fi=math.atan(z/(r*(1-e2))) #radiany
            fi2=2*fi
            while np.abs(fi2-fi)>eps:
                fi=fi2
                
            
                N=self.a/math.sqrt(1-self.ecc2*math.sin(fi2)**2)
                
            
                h=(r/np.cos(fi) - N)
                
            
                fi2=math.atan(z/(r*(1-self.ecc2* (N/(N+h))))) #radiany
                
        
            N=self.a/np.sqrt(1-self.ecc2*np.sin(fi2)**2)
            h=r/math.cos(fi) - N
            lam = math.atan(y/x)
            return fi2, lam, h , N
        
        
        def blh2xyz(self, fi, lam, h):
                N = self.a / math.sqrt(1 - self.ecc2 * (math.sin(fi)) ** 2)
                X = (N + h) * math.cos(fi) * math.cos(lam)
                Y = (N + h) * math.cos(fi)*math.sin(lam)
                Z = (N * (1 - self.ecc2) + h) * math.sin(fi)
                return X, Y, Z
        
        def s_A_z2neu(s, A, z):
            A = np.deg2rad(A)
            z = np.deg2rad(z)
            n = s*np.sin(z)*np.cos(A)
            e = s*np.sin(z)*np.sin(A)
            u = s*np.cos(z)
            return n, e, u
        
        def neu2dXYZ(n, e, u, F1, L1):
        
             R = np.array([[-np.sin(F1) * np.cos(L1), -np.sin(L1), np.cos(F1)*np.cos(L1)],
                  [-np.sin(F1)*np.sin(L1), np.cos(L1), np.cos(F1)*np.sin(L1)],
                  [np.cos(F1), 0, np.sin(F1)]])
             dx = np.linalg.inv(R.transpose()) @ np.array([n, e, u])
             return dx
    
        dx = neu2dXYZ(n, e, u, fi, lam)
        
        
        
       
        e_2 = self.ecc2 / (1-self.ecc2)
  
        
        def u1992(self,fi, lam, a, e2, m_0):
            N = self.a/(math.sqrt(1-self.ecc2 * np.sin(fi)**2))
            t = np.tan(fi)
            n2 = e_2 * np.cos(lam)**2
            lam_0 = math.radians(19) #poczatek ukladu w punkcie przeciecia poludnika L0 = 19st z obrazem równika 
            l = lam - lam_0
            
            A_0 = 1 - (self.ecc/4) - (3*(self.ecc2**2))/64 - (5*(self.ecc2**3))/256
            A_2 = 3/8 * (self.ecc2 + ((self.ecc2**2)/4) + ((15*self.ecc2**3)/128))
            A_4 = 15/256 * (self.ecc2**2 + (3*(self.ecc2**3))/4)
            A_6 = (35*(self.ecc2**3))/3072
            
            sigma = self.a* ((A_0*fi) - (A_2*np.sin(2*fi)) + (A_4*np.sin(4*fi)) - (A_6*np.sin(6*fi)))
            
            x = sigma + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
            y = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
        
            x92 = round(x * m_0 - 5300000, 3)
            y92 = round(y * m_0 + 500000, 3)   
            
            return x92, y92 
        
        
        
        
        
        def u2000(self,fi, lam, a, e2, m_0):
            N = self.a/(math.sqrt(1-self.ecc2 * np.sin(fi)**2))
            t = np.tan(fi)
            n2 = self.ecc2 * np.cos(lam)**2
            lam = math.degrees(lam)
            
            if lam > 13.5 and lam < 16.5:
                s = 5
                lam_0 = 15
            elif lam > 16.5 and lam < 19.5:
                s = 6
                lam_0 = 18
            elif lam > 19.5 and lam < 22.5:
                s = 7
                lam_0 = 21
            elif lam > 22.5 and lam < 25.5:
                s = 8
                lam_0 = 24
                
            lam = math.radians(lam)
            lam_0 = math.radians(lam_0)
            l = lam - lam_0
            
            A_0 = 1 - (self.ecc2/4) - (3*(self.ecc2**2))/64 - (5*(self.ecc2**3))/256
            A_2 = 3/8 * (self.ecc2 + ((self.ecc2**2)/4) + ((15*self.ecc2**3)/128))
            A_4 = 15/256 * (self.ecc2**2 + (3*(self.ecc2**3))/4)
            A_6 = (35*(self.ecc2**3))/3072
            
            
            sigma = self.a* ((A_0*fi) - (A_2*np.sin(2*fi)) + (A_4*np.sin(4*fi)) - (A_6*np.sin(6*fi)))
            
            x = sigma + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
            y = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
        
            x00 = round(x * m_0, 3)
            y00 = round(y * m_0 + (s*1000000) + 500000, 3)   
            
            return x00, y00 
    
