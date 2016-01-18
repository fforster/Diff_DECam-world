#!/usr/local/bin/python
# coding: latin-1

#import math as m
#import astroTools as astro

# Convert deg to arcsec:
def conv ( Decin ):
   
   if(Decin<0):
      sign = -1
      dec  = -Decin
   else:
      sign = 1
      dec  = Decin

   d = int( dec )
   dec -= d
   dec *= 100.
   m = int( dec*3./5. )
   dec -= m*5./3.
   s = dec*180./5.
   
   sec = float(d)*3600. + float(m)*60. + s 

   if(sign == -1):
      out = -sec
   else: out = sec
   
   print out
    
   return out
   
if __name__ == "__main__" :
	import sys
	conv (float (sys.argv[1]))   

