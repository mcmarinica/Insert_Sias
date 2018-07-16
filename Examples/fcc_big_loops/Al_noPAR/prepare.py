import array
import math
import os
import shutil 
import sys
import re
sys.path.append(os.environ['HOME']+'/bin')
from cospy import get_file_len,get_files_from_dir, \
                  search_word_replace_line

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
      python prepare_vasp.py  base 
  
      base - where the script is launched
      default - the current directory 
      
"""

if (len(sys.argv) > 2):
   print usage
   exit(0)
if (len(sys.argv)==2):
   base=sys.argv[1]
if (len(sys.argv)==1):
   base=os.getcwd()


list_of_poscar_with_h=get_files_from_dir(base+'/'+'POSCAR_ur')


vasp_input=base+'/INP_VASP'
if  not (os.path.isdir(vasp_input) or os.path.islink(vasp_input)):
  print '% is is not present in the path . Put the correct path for base'%vasp_input
  exit(0)

files_to_test=['KPOINTS','INCAR','jsub','POSCAR']

for tmp in files_to_test:
  if  not (os.path.isfile(vasp_input+'/'+tmp) or os.path.islink(vasp_input+'/'+tmpn)):
    print '% is is not present in the path . Put the correct path for base'%tmp
    exit(0)

# get natom_matrix ... number of atoms in matrix ....

matrix_poscar=base+'/INP_VASP'+'/'+'POSCAR'
natom_matrix,nnline_matrix=get_natom_from_vasp_incar(matrix_poscar)

# get the position of the matrix ...
fmatrix=open(matrix_poscar,'r')
idirect=0
itest=0
iline=0
positions_matrix=[]
for line in fmatrix:
#     print line
     tmpa=line
     ln=line.strip()
     if (len(ln) > 0):
      if ( (ln[0] == 'd') or (ln[0]=='D') or (ln[0]=='c') or (ln[0]=='C') or (ln[0]=='k') or (ln[0]=='K') ) :
       itest=1
       iline=idirect
      if((itest==1) and (idirect > iline) ):
#       data="".join(line.rstrip())
       positions_matrix.append(line.rstrip())
      idirect += 1
if (itest==0):
     print 'The Direct POSCAR matrix is not accepted in the input'
     print('Change the file %s'%(matrix_poscar))
     exit(0)
#print positions_matrix
fmatrix.close()



for tmp in list_of_poscar_with_h:

  os.chdir(base)

  f=open(base+'/'+'POSCAR_ur'+'/'+tmp)
# get the last line ... H coordinates.
  ttmp=str(tmp).split(".")[:]
#  print tmp, ttmp
  os.makedirs('all'+'/'+ttmp[0])
  os.chdir('all'+'/'+ttmp[0])

  last_line=tail(f,2)
#  print last_line[0]
  f.seek(0)

 
  first=headn(base+'/'+'POSCAR_ur'+'/'+tmp,5)
  fout=open('POSCAR','w')
  for i in range(5):
   fout.write(first[i])
  fout.write('H   W \n')
  fout.write('1  %s \n'%(natom_matrix))
  fout.write('Direct \n')
  fout.write(last_line[0]+' \n')
  for i in range(len(positions_matrix)):
   fout.write('%s \n'%(positions_matrix[i]))  
  
  f.close()
  fout.close()


# get the INCAR, POTCAR, KPOINTS and jsub
  os.symlink(vasp_input+'/'+'INCAR','INCAR')
  os.symlink(vasp_input+'/'+'POTCAR','POTCAR')
  os.symlink(vasp_input+'/'+'KPOINTS','KPOINTS')
  shutil.copy(vasp_input+'/'+'jsub','jsub')

  os.chdir(base) 
