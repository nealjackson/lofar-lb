import os,sys,numpy as np
os.system('source /home/njj/.cshrc')
os.system('source /aips/LOGIN.CSH')
os.system('sed "s/source[i]/sources[i]/g" proc_pmkmap.py >proc_pmkmap1.py')
os.system('grep sources proc_pmkmap1.py')
os.system('which parseltongue')
os.system('parseltongue proc_pmkmap1.py ../newdata/44059320F.fits')

