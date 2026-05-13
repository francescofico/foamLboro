import os
import numpy as np
import classy_blocks as cb

mesh = cb.Mesh()

x = 0.5
y = 0.05
z = 0.05
x_elem = 50

start_size_h = 1e-4 # first cell hartmann wall
start_size_s = 1e-4 # first cell side wall
count_h = 20 # hartmann
count_s =  20# side wall
count_c = 25 # core

coarse_box = cb.Box([0, 0, 0], [x, y, z])
    
coarse_box.chop(0,count=x_elem)

coarse_box.chop(1, count=count_h, length_ratio=0.5, start_size=start_size_h, invert=False)
#coarse_box.chop(1, count=count_c, c2c_expansion=1, length_ratio=0.8, invert=False)
coarse_box.chop(1, count=count_h, length_ratio=0.5, start_size=start_size_h, invert=True)

coarse_box.chop(2, count=count_s, length_ratio=0.5, start_size=start_size_s, invert=False)
#coarse_box.chop(2, count=count_c, c2c_expansion=1, length_ratio=0.8, invert=False)
coarse_box.chop(2, count=count_s, length_ratio=0.5, start_size=start_size_s, invert=True)



coarse_box.set_patch(['bottom','top','front','back'], 'fluid_surface')
coarse_box.set_patch('right', 'outlet')
coarse_box.set_patch('left', 'inlet')
"first duct"

mesh.add(coarse_box)
mesh.write(os.path.join (".","system", "blockMeshDict"))

