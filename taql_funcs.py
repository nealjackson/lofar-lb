import os
def taql_calc (vis, vistable, qtab, qtype):
    os.system('taql \'CALC '+qtype+' ([select '+qtab+' from '+vis+\
              '/'+vistable+'])\' >taql_out')
    f=open('taql_out')
    for v in f:
        try:
            val = float(v.rstrip('\n'))
        except:
            pass
    f.close(); os.system('rm taql_out')
    return val

def taql_num (vis,vistable,qtab):
    os.system('taql \'select '+qtab+' from '+vis+'/'+vistable+'\' >taql_out')
    f=open('taql_out')
    for v in f:
        if 'select result of' in v:
            n = int(v.split('of')[1].split('row')[0])
            break
    f.close(); os.system('rm taql_out')
    return n

def taql_from (vis,vistable,qtab):
    os.system('taql \'select '+qtab+' from '+vis+'/'+vistable+'\' >taql_out')
    f = open('taql_out')
    for v in f:
        pass
    f.close(); os.system('rm taql_out')
    return v.rstrip('\n').rstrip(']').lstrip('[').split(',')

