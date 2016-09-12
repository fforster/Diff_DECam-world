#!/usr/bin/python2.7

'''
This module (may) contains the following functions:
- jack_res
accepts a sample and performs a left-one out resampling
- jackknife 
call jack_res for the arrays to resample
- jackCI
accepts the flux computed for the original sample and all the fluxes computed for all the sub-samples 
to compute the confidence interval at 95%
- bootstrap
accepts a sample and performs a reampling with repetition
- bootCI 
accepts the flux computed for the original sample and all the fluxes computed for all the random samples 
to compute the confidence interval at 95%
'''


import os
import sys
import numpy as np
import numpy.random as npr
from astropy.stats import bootstrap


def jack_res (v) :
	
    #print np.shape(v)
	
    if v.ndim == 1 :
        for i in range(len(v)) :
            if i==0 :
                matv = np.delete (v, i)
            else :
                vtemp = np.delete (v, i)
                matv = np.vstack ((matv, vtemp))
        matvsplit = np.split (matv, len(v))
        return [arr.flatten() for arr in matvsplit]
    
    elif v.ndim == 3 :
        for i in range(np.shape(v)[2]) :
            if i==0 :
                matv = np.delete (v, i, axis=2)
            else :
                vtemp = np.delete (v, i, axis=2)
                matv = np.vstack ((matv, vtemp))
        return np.split (matv, np.shape(v)[2])


def jackknife (psf1s, psf2s, f1sel, e_f1sel, f2sel, e_f2sel, r1sel, r2sel, list_isource) :
	
	res_psf1s = jack_res (psf1s)
	res_psf2s = jack_res (psf2s)
	res_f1sel = jack_res (f1sel)
	res_e_f1sel = jack_res (e_f1sel)
	res_f2sel = jack_res (f2sel)
	res_e_f2sel = jack_res (e_f2sel)
	res_r1sel = jack_res (r1sel)
	res_r2sel = jack_res (r2sel)
	res_isourcelist = jack_res (list_isource)
	
	return res_psf1s, res_psf2s, res_f1sel, res_e_f1sel, res_f2sel, res_e_f2sel, res_r1sel, res_r2sel, res_isourcelist


if __name__ == "__main__" :
	
	jackknife ()
