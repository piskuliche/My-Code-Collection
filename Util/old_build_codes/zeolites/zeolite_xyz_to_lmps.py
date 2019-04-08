#!/usr/bin/env python
import numpy as np
import sys 

def enterbounds(name):
    try:
        val = float(raw_input('Enter %s:\n' % name))
    except ValueError:
        print("Not a floating point number")
        exit(1)
    return val

arg = sys.argv
inp_file=""


if len(arg) == 3:
    inp_file=str(arg[1])
    out_file=str(arg[2])
else:
    print("Usage: python zeolite_xyz_to_lmps.py inp_file out_file")
    exit(1)

with open(inp_file) as f:
    lines=f.readlines()


natms=int(lines[0].rstrip())
type=[]
rx=[]
ry=[]
rz=[]

xlo=enterbounds('xlo')
xhi=enterbounds('xhi')
ylo=enterbounds('ylo')
yhi=enterbounds('yhi')
zlo=enterbounds('zlo')
zhi=enterbounds('zhi')

print("There are %s atoms" % natms)
for i in range(2,int(natms)+2):
    type.append(lines[i].rstrip().split()[0])
    rx.append(lines[i].rstrip().split()[1])
    ry.append(lines[i].rstrip().split()[2])
    rz.append(lines[i].rstrip().split()[3])

out = open(out_file,'w')
natypes=len(np.unique(type))
print np.unique(type)
out.write('Lammps Description\n\n')

out.write('%s atoms\n' % natms)
out.write('0 bonds\n' )
out.write('0 angles\n')
out.write('0 dihedrals\n')
out.write('0 impropers\n\n')


out.write('%s atom types\n' % natypes )
out.write('0 bond types\n')
out.write('0 angle types\n')
out.write('0 dihedral types\n')
out.write('0 improper types\n\n')

out.write('%s %s xlo xhi\n' % (xlo,xhi))
out.write('%s %s ylo yhi\n' % (ylo,yhi))
out.write('%s %s zlo zhi\n\n' % (zlo,zhi))

out.write('Masses\n\n')
count=0
q=[]
for elem in np.unique(type):
    count += 1
    if elem == 'Si':
        q.append(2.050)
        out.write('%s 28.0855\n' % count)
    elif elem =='O':
        q.append(-1.025)
        out.write('%s 15.999\n' % count)
    elif elem =='C4':
        q.append(0.0)
        out.write('%s 16.043\n' % count)
out.write('\n')

out.write('Atoms\n\n')

for i in range(natms):
    aindex = i+1
    molindex = 1
    atype = np.unique(type).tolist().index(type[i])+1
    Q = q[atype-1]
    out.write('%s %s %s %2.4f   %3.6f   %3.6f   %3.6f\n' % (aindex,molindex,atype,float(Q),float(rx[i]),float(ry[i]),float(rz[i])))

