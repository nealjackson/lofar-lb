import os,sys
fo=open('doit','w')
fo.write('#!/bin/csh\nsource /home/njj/.cshrc\nsource /aips/LOGIN.CSH\ncp /mirror3/scratch/scratch/njj/rqq/scripts/proc_casaimg.py /home/njj/public_html\n')
fo.close()
os.system('chmod 755 doit')
os.system('./doit')
os.system('rm doit')
os.system('rm doit_frastra.py*')

