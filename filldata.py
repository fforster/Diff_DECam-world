#!/usr/bin/python2.7

import os
import re # use regular patterns
import sys, getopt # system commands
import subprocess

address = 'development.dim.uchile.cl'
outdir = '/media/fforster_fondecyt_2Tb_2014/DATA/DATA_CMMPIPE'
indir = "/mnt/astro-arc/DECAM-14A-P"
field = sys.argv[1]
CCD = sys.argv[2]
ref = int(sys.argv[3])

# check available dates
ls = subprocess.Popen(['ssh','development.dim.uchile.cl', 'ls %s/ARCHIVE/%s/' % (indir, field)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err =  ls.communicate()
if err != '':
    print "Cannot list remote directory"
    sys.exit()

# find matching MJD according to .mapped files
MJDs = re.findall("(.*?)\n", out)
for iMJD in range(len(MJDs)):
    ls = subprocess.Popen(['ssh', address, 'cat %s/ARCHIVE/%s/%s/.mapped' % (indir, field, MJDs[iMJD])], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = ls.communicate()
    if int(out) == ref:
        print "Found match for ref %s in date %s" % (ref, MJDs[iMJD])
        break

# download data
savedir = "%s/%s/%s" % (outdir, field, CCD)
if not os.path.exists(savedir):
    os.makedirs(savedir)
command = "rsync -trvz %s:%s/ARCHIVE/%s/%s/%s/%s_%s_%s_image_crblaster.fits %s/%s_%s_%02i_image_crblaster.fits" % (address, indir, field, MJDs[iMJD], CCD, field, CCD, MJDs[iMJD], savedir, field, CCD, ref)
print command
os.system(command)

command = "rsync -trvz %s:%s/SHARED/%s/%s/%s_%s_%02i_image_crblaster.fits-catalogue_wtmap.dat %s/" % (address, indir, field, CCD, field, CCD, ref, savedir)
print command
os.system(command)
