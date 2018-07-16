import array
import math
import os
import shutil 
import sys
import re
sys.path.append(os.environ['HOME']+'/bin')
from cospy import get_file_len,get_files_from_dir, \
                  search_word_replace_line

def replace(oldfile, newfile, pattern, subst):
    #Create temp file
    new_file = open(newfile,'w')
    old_file = open(oldfile)
    for line in old_file:
        new_file.write(line.replace(pattern, subst))
    #close temp file
    new_file.close()
    old_file.close()
    #Remove original file
    #Move new file


def head(fname,nline):

   f=open(fname,'r')

   head=[next(f) for x in xrange(nline) ]
   f.close()

   return ''.join(data).splitlines() 
   
def headn(file_name, n):
    result = []
    nlines = 0
    assert n >= 1
    for line in open(file_name):
        result.append(line)
        nlines += 1
        if nlines >= n:
            break
    return result



def tail(f, window):
    """
    Returns the last `window` lines of file `f` as a list.
    """
    if window == 0:
        return []
    BUFSIZ = 1024
    f.seek(0, 2)
    bytes = f.tell()
    size = window + 1
    block = -1
    data = []
    while size > 0 and bytes > 0:
        if bytes - BUFSIZ > 0:
            # Seek back one whole BUFSIZ
            f.seek(block * BUFSIZ, 2)
            # read BUFFER
            data.insert(0, f.read(BUFSIZ))
        else:
            # file too small, start from begining
            f.seek(0,0)
            # only read what was not read
            data.insert(0, f.read(bytes))
        linesFound = data[0].count('\n')
        size -= linesFound
        bytes -= BUFSIZ
        block -= 1
    return ''.join(data).splitlines()[-window:]


def get_natom_from_vasp_incar(file_poscar):

  fil=open(file_poscar,'r')
  natom=0
  for line in fil:
   tmpa=line
   ln=line.strip()
   if ( (ln[0] == 'd') or (ln[0]=='D') or (ln[0]=='c') or (ln[0]=='C') or (ln[0]=='k') or (ln[0]=='K') ) :
     nat_line=tmpb.strip().split()
     nspecies=len(nat_line)
     for i in range(nspecies):
      natom = natom + int(nat_line[i])
     break

   tmpb=line
  return natom,nat_line
 

usage="""
      python prepare_vasp.py  dir_with_gin  dir_run 
  
      dir_with_gin - where the script is launch
      dir_run      - where the run directories are created
      
"""


if (len(sys.argv) > 3):
   print usage
   exit(0)
if (len(sys.argv)!= 3):
   print usage
   exit(0)

nbulk=500
base=os.getcwd()
dir_with_gin=sys.argv[1]
dir_out=sys.argv[2]
list_of_gin=get_files_from_dir(dir_with_gin)


vasp_input=base+'/INP_VASP'
if  not (os.path.isdir(vasp_input) or os.path.islink(vasp_input)):
  print '%s is is not present in the path . Put the correct path for base'%vasp_input
  exit(0)

files_to_test=['KPOINTS.vasp','INCAR.vasp','jsub','POTCAR.vasp']

for tmp in files_to_test:
  if  not (os.path.isfile(vasp_input+'/'+tmp) or os.path.islink(vasp_input+'/'+tmp)):
    print '%s is is not present in the path . Put the correct path for base'%tmp
    exit(0)

exe_to_test=['gin2poscar','gin2selectiveposcar','read_write_scale_gin.pl']

for tmp in files_to_test:
  if  not (os.path.isfile(vasp_input+'/'+tmp) or os.path.islink(vasp_input+'/'+tmp)):
    print '%s is is not present in the path . Put the correct path for base'%tmp
    exit(0)


os.system('mkdir -p %s'%dir_out)

for tmp in list_of_gin:
   os.chdir(base+'/'+dir_out)
   os.system('mkdir -p %s'%(tmp))
   os.system('cp -r %s/INCAR.vasp   %s/INCAR'%(vasp_input,tmp))
   os.system('cp -r %s/KPOINTS.vasp %s/KPOINTS'%(vasp_input,tmp))
   os.system('cp -r %s/POTCAR.vasp  %s/POTCAR'%(vasp_input,tmp))
   os.system('cp -r %s/jsub  %s/jtmp'%(vasp_input,tmp))
   os.system('%s/read_write_scale_gin.pl  %s  2.8553 4.0397375 > %s/l.gin'%(vasp_input,base+'/'+dir_with_gin+'/'+tmp, tmp))
   os.system('%s/gin2poscar %s/l.gin %s/POSCAR'%(vasp_input,tmp,tmp))
   os.chdir(tmp)
   ttmp=str(tmp).split("_")[:]
   replace('jtmp','jtmp01','SSUFF',dir_with_gin)
   replace('jtmp01','jsub','NNAMEE','j'+ttmp[0]+'_'+ttmp[1])
   ntot,nat_line=get_natom_from_vasp_incar('POSCAR')
   nsia=ntot-nbulk
   if nsia != 0:
      os.system('echo %s > tt'%nsia)
      os.system('%s/gin2selectiveposcar l.gin POSCAR < tt'%(vasp_input))
#   os.system('ccc_msub jsub')



exit(0)

# get natom_matrix ... number of atoms in matrix ....


