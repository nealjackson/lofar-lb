import os,sys
fo=open('doit','w')
fo.write('#!/bin/csh\nsource /home/njj/.cshrc\nsource /aips/LOGIN.CSH\n')
fo.write('wget https://raw.githubusercontent.com/nealjackson/lofar-lb/master/proc_casaimg.py\n')
fo.write('parseltongue proc_casaimg.py')
fo.close()
os.system('chmod 755 doit')
os.system('./doit')
os.system('rm doit')
os.system('rm doit_frastra.py*')

