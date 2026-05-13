#%%
import os
import numpy as np
import classy_blocks as cb

mesh = cb.Mesh()

x = 0.15
y = 0.015
z = 0.015
x_elem = 50 # orig 60


start_size_x = 0.5e-4 
start_size_h = 1e-04 # first cell h - orig 5e-06
start_size_s = 1e-04 # first cell s - orig 5e-06
count_h = 14 # hartmann - orig 32
count_s = 14# side wall - orig 32
count_c = 16 # core - orig 52

coarse_box = cb.Box([-x, -y, -z], [x, y, z])
    
coarse_box.chop(0,count=x_elem)

coarse_box.chop(1, count=count_h, length_ratio=0.2, start_size=start_size_h, invert=False)
coarse_box.chop(1, count=count_c, c2c_expansion=1, length_ratio=0.6, invert=False)
coarse_box.chop(1, count=count_h, length_ratio=0.2, start_size=start_size_h, invert=True)

coarse_box.chop(2, count=count_s, length_ratio=0.2, start_size=start_size_s, invert=False)
coarse_box.chop(2, count=count_c, c2c_expansion=1, length_ratio=0.6, invert=False)
coarse_box.chop(2, count=count_s, length_ratio=0.2, start_size=start_size_s, invert=True)

#coarse_box.set_patch(['bottom','top','front','back'], 'fluid_surface')
coarse_box.set_patch('top', 'top')
coarse_box.set_patch('front', 'front')
coarse_box.set_patch('back', 'back')
coarse_box.set_patch('bottom', 'bottom')
coarse_box.set_patch('right', 'outlet')
coarse_box.set_patch('left', 'inlet')
"first duct"

mesh.add(coarse_box)
mesh.write(os.path.join (".","system", "fluid", "blockMeshDict"))


# %%
