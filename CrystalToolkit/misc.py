
for i in `cat re_calc_dirs_list | cut -d ' ' -f 12 | cut -d '/' -f1` ; do rm hec5_20sqs_relax/$i/OUTCAR ; echo `pwd`/hec5_20sqs_relax/$i >> /workspace/RZ/RunJobs/CalcJobList ; done

for i in `head -n 10 re_calc_dirs_list | cut -d ' ' -f 12 | cut -d '/' -f1` ; do echo `pwd`/hec5_20sqs_relax/$i ; done

from pymatgen.transformations.advanced_transformations import EnumerateStructureTransformation
import pymatgen as mg
nacl = mg.Structure.from_spacegroup(
     'Fm-3m', mg.Lattice.cubic(5.6),['Na', 'Cl'],[[0.5, 0.5, 0.5], [0, 0, 0]] )

s = nacl.get_primitive_structure()
s['Na']='Na0.5K0.5'

e = EnumerateStructureTransformation(max_cell_size=2)
e.apply_transformation(s)


sc = StructureMatcher()
for i,d1 in enumerate(d):
    s1 = d1['structure']
    for j,s2 in enumerate(s_49):
        if sc.fit(s1, s2):
            print(i,j)
            break






>>> os.chdir('..')
>>> os.getcwd()
'/media/nfs_1/workspace/RZ'
>>> os.mkdir('hec5')
>>> os.chdir('hec5')
>>> os.getcwd()
'/media/nfs_1/workspace/RZ/hec5'
>>> for i,struc in enumerate(c):
...     struc[1].to(filename=str(i)+'.vasp',fmt='poscar')
... 
Traceback (most recent call last):
  File "<stdin>", line 2, in <module>
KeyError: 1
>>> for i,struc in enumerate(c):
...     struc['structure'].to(filename=str(i)+'.vasp',fmt='poscar')
... 
>>> os.getcwd()
'/media/nfs_1/workspace/RZ/hec5'
>>> os.chdir('../hec5_49')
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
FileNotFoundError: [Errno 2] No such file or directory: '../hec5_49'
>>> os.chdir('../hec5_49')
>>> s_49 = []
>>> for i in os.listdir('.'):
...     s_49.append(
... 
... 
... )
... 
Traceback (most recent call last):
  File "<stdin>", line 2, in <module>
TypeError: append() takes exactly one argument (0 given)
>>> from pymatgen import Structure
>>> for i in os.listdir('.'):
...     s_49.append(Structure.from_file(i))
... 
>>> s_49[0]
Structure Summary
Lattice
    abc : 10.58190958195514 10.58190958195514 10.58190958195514
 angles : 164.05762953971183 164.05762953971183 22.619867597024125
 volume : 89.37832582062589
      A : -1.467447 1.467447 10.376416
      B : 1.467447 -1.467447 10.376416
      C : 1.467447 1.467447 -10.376416
PeriodicSite: Ti (0.0000, 0.0000, 6.2258) [0.3000, 0.3000, -0.0000]
PeriodicSite: V (0.0000, 0.0000, 18.6775) [0.9000, 0.9000, -0.0000]
PeriodicSite: Al (0.0000, 0.0000, 10.3764) [0.5000, 0.5000, -0.0000]
PeriodicSite: Si (0.0000, 0.0000, 2.0753) [0.1000, 0.1000, -0.0000]
PeriodicSite: Cr (0.0000, 0.0000, 14.5270) [0.7000, 0.7000, -0.0000]
PeriodicSite: N (0.0000, 0.0000, 8.3011) [0.4000, 0.4000, 0.0000]
PeriodicSite: N (0.0000, 0.0000, 12.4517) [0.6000, 0.6000, 0.0000]
PeriodicSite: N (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]
PeriodicSite: N (0.0000, 0.0000, 16.6023) [0.8000, 0.8000, 0.0000]
PeriodicSite: N (0.0000, 0.0000, 4.1506) [0.2000, 0.2000, 0.0000]
>>> import numpy as np
>>> 89.37 ** 0.33333
4.47085668556518
>>> s_0 = s_49[0]
>>> s_0
Structure Summary
Lattice
    abc : 10.58190958195514 10.58190958195514 10.58190958195514
 angles : 164.05762953971183 164.05762953971183 22.619867597024125
 volume : 89.37832582062589
      A : -1.467447 1.467447 10.376416
      B : 1.467447 -1.467447 10.376416
      C : 1.467447 1.467447 -10.376416
PeriodicSite: Ti (0.0000, 0.0000, 6.2258) [0.3000, 0.3000, -0.0000]
PeriodicSite: V (0.0000, 0.0000, 18.6775) [0.9000, 0.9000, -0.0000]
PeriodicSite: Al (0.0000, 0.0000, 10.3764) [0.5000, 0.5000, -0.0000]
PeriodicSite: Si (0.0000, 0.0000, 2.0753) [0.1000, 0.1000, -0.0000]
PeriodicSite: Cr (0.0000, 0.0000, 14.5270) [0.7000, 0.7000, -0.0000]
PeriodicSite: N (0.0000, 0.0000, 8.3011) [0.4000, 0.4000, 0.0000]
PeriodicSite: N (0.0000, 0.0000, 12.4517) [0.6000, 0.6000, 0.0000]
PeriodicSite: N (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]
PeriodicSite: N (0.0000, 0.0000, 16.6023) [0.8000, 0.8000, 0.0000]
PeriodicSite: N (0.0000, 0.0000, 4.1506) [0.2000, 0.2000, 0.0000]
>>> for i in ['V', 'Al', 'Si', 'Cr']:
...     s_0[i] = 'Ti'
... 
>>> s_0
Structure Summary
Lattice
    abc : 10.58190958195514 10.58190958195514 10.58190958195514
 angles : 164.05762953971183 164.05762953971183 22.619867597024125
 volume : 89.37832582062589
      A : -1.467447 1.467447 10.376416
      B : 1.467447 -1.467447 10.376416
      C : 1.467447 1.467447 -10.376416
PeriodicSite: Ti (0.0000, 0.0000, 6.2258) [0.3000, 0.3000, -0.0000]
PeriodicSite: Ti (0.0000, 0.0000, 18.6775) [0.9000, 0.9000, -0.0000]
PeriodicSite: Ti (0.0000, 0.0000, 10.3764) [0.5000, 0.5000, -0.0000]
PeriodicSite: Ti (0.0000, 0.0000, 2.0753) [0.1000, 0.1000, -0.0000]
PeriodicSite: Ti (0.0000, 0.0000, 14.5270) [0.7000, 0.7000, -0.0000]
PeriodicSite: N (0.0000, 0.0000, 8.3011) [0.4000, 0.4000, 0.0000]
PeriodicSite: N (0.0000, 0.0000, 12.4517) [0.6000, 0.6000, 0.0000]
PeriodicSite: N (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]
PeriodicSite: N (0.0000, 0.0000, 16.6023) [0.8000, 0.8000, 0.0000]
PeriodicSite: N (0.0000, 0.0000, 4.1506) [0.2000, 0.2000, 0.0000]
>>> s0 = s_0.get_primitive_structure()
>>> s0
Structure Summary
Lattice
    abc : 2.934894 2.934893823602523 2.934893823602522
 angles : 60.000003976417545 59.999998011791256 59.999998011791256
 volume : 17.875665164125174
      A : 0.0 -2.934894 1.7763568394002505e-15
      B : -1.467447 -1.467447 -2.0752831999999994
      C : 1.467447 -1.467447 -2.0752831999999977
PeriodicSite: N (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]
PeriodicSite: Ti (0.0000, -2.9349, -2.0753) [0.5000, 0.5000, 0.5000]
>>> s0['Ti'] = 'Ti0.2V0.2Al0.2Si0.2Cr0.2'
>>> s0
Structure Summary
Lattice
    abc : 2.934894 2.934893823602523 2.934893823602522
 angles : 60.000003976417545 59.999998011791256 59.999998011791256
 volume : 17.875665164125174
      A : 0.0 -2.934894 1.7763568394002505e-15
      B : -1.467447 -1.467447 -2.0752831999999994
      C : 1.467447 -1.467447 -2.0752831999999977
PeriodicSite: N (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]
PeriodicSite: Ti:0.200, Al:0.200, V:0.200, Cr:0.200, Si:0.200 (0.0000, -2.9349, -2.0753) [0.5000, 0.5000, 0.5000]
>>> d = e.apply_transformation(s0,return_ranked_list=100)
>>> len(d)
54
>>> from pymatgen.analysis.structure_matcher import StructureMatcher
>>> for d1 in d:
...     s1 = d1['structure']
...     for s2 in s_49:
...             if StructureMatcher(s1, s2)









