#%%
import os
import numpy as np
import classy_blocks as cb

mesh = cb.Mesh()

x = 0.15

# angle 1
y1 = -0.017
z1 = -0.017
y_1_2 = -0.015
z_1_2 = -0.015

x_elem = 50 # orig 60



start_size_x = 0.5e-4 
start_size_h = 1e-04 # first cell h - orig 5e-06
start_size_s = 1e-04 # first cell s - orig 5e-06
count_h = 14 # hartmann - orig 32
#count_s = 12# side wall - orig 32
count_c = 16 # core - orig 52

"""
box definition - 4 angles, 4 connections - hex mesh
"""
box_angle1 = cb.Box([-x, -0.017, -0.017], [x, -0.015, -0.015])
box_angle2 = cb.Box([-x, -0.017, 0.015], [x, -0.015, 0.017])
box_angle3 = cb.Box([-x, 0.015, -0.017], [x, 0.017, -0.015])
box_angle4 = cb.Box([-x, 0.015, 0.015], [x, 0.017, 0.017])


box_12 = cb.Box([-x, -0.015, -0.017], [x, 0.015, -0.015 ])
box_13 = cb.Box([-x, -0.017, -0.015], [x, -0.015, 0.015 ])
box_23 = cb.Box([-x, -0.015, 0.015], [x, 0.015, 0.017 ])
box_34 = cb.Box([-x, 0.015, -0.015], [x, 0.017, 0.015 ])

"""
box angles 1,2,3,4 mesh
"""

box_angle1.chop(0, count=x_elem)
box_angle1.chop(1, count=count_h, length_ratio=1, start_size=start_size_h, invert=True)
box_angle1.chop(2, count=count_h, length_ratio=1, start_size=start_size_h, invert=True)
# -------------

box_angle2.chop(0, count=x_elem)
box_angle2.chop(1, count=count_h, length_ratio=1, start_size=start_size_h, invert=True)
box_angle2.chop(2, count=count_h, length_ratio=1, start_size=start_size_h, invert=False)

# --------------
box_angle3.chop(0, count=x_elem)
box_angle3.chop(1, count=count_h, length_ratio=1, start_size=start_size_h, invert=False)
box_angle3.chop(2, count=count_h, length_ratio=1, start_size=start_size_h, invert=True)


#----------------
box_angle4.chop(0, count=x_elem)
box_angle4.chop(1, count=count_h, length_ratio=1, start_size=start_size_h, invert=False)
box_angle4.chop(2, count=count_h, length_ratio=1, start_size=start_size_h, invert=False)
#----------------

"""
box connections 12,13,23,34 mesh
"""
box_12.chop(0, count=x_elem)
box_12.chop(1, count=count_h, length_ratio=0.2, start_size=start_size_h, invert=False)
box_12.chop(1, count=count_c, c2c_expansion=1, length_ratio=0.6, invert=False)
box_12.chop(1, count=count_h, length_ratio=0.2, start_size=start_size_h, invert=True)

box_12.chop(2, count=count_h, length_ratio=1, start_size=start_size_h, invert=True)


# --------------------

box_13.chop(0, count=x_elem)
box_13.chop(1, count=count_h, length_ratio=1, start_size=start_size_h, invert=True)

box_13.chop(2, count=count_h, length_ratio=0.2, start_size=start_size_h, invert=False)
box_13.chop(2, count=count_c, c2c_expansion=1, length_ratio=0.6, invert=False)
box_13.chop(2, count=count_h, length_ratio=0.2, start_size=start_size_h, invert=True)


# --------------------
box_23.chop(0, count=x_elem)
box_23.chop(1, count=count_h, length_ratio=0.2, start_size=start_size_h, invert=False)
box_23.chop(1, count=count_c, c2c_expansion=1, length_ratio=0.6, invert=False)
box_23.chop(1, count=count_h, length_ratio=0.2, start_size=start_size_h, invert=True)
box_23.chop(2, count=count_h, length_ratio=1, start_size=start_size_h, invert=False)

# --------------------
box_34.chop(0, count=x_elem)
box_34.chop(1, count=count_h, length_ratio=1, start_size=start_size_h, invert=False)

box_34.chop(2, count=count_h, length_ratio=0.2, start_size=start_size_h, invert=False)
box_34.chop(2, count=count_c, c2c_expansion=1, length_ratio=0.6, invert=False)
box_34.chop(2, count=count_h, length_ratio=0.2, start_size=start_size_h, invert=True)

# --------------------

# box_angle1.set_patch('back', 'back1')
# box_angle3.set_patch('front', 'front3')
# box_12.set_patch('front', 'mix-back1')
# box_12.set_patch('back', 'mix-front1')

box_angle1.set_patch('right', 'right')
box_angle1.set_patch('left', 'left')
box_angle2.set_patch('right', 'right')
box_angle2.set_patch('left', 'left')
box_angle3.set_patch('right', 'right')
box_angle3.set_patch('left', 'left')
box_angle4.set_patch('right', 'right')
box_angle4.set_patch('left', 'left')

box_12.set_patch('top', 'bottom')
box_13.set_patch('back', 'back')
box_23.set_patch('bottom', 'top')
box_34.set_patch('front', 'front')

box_12.set_patch('right', 'right')
box_12.set_patch('left', 'left')
box_13.set_patch('right', 'right')
box_13.set_patch('left', 'left')
box_23.set_patch('right', 'right')
box_23.set_patch('left', 'left')
box_34.set_patch('right', 'right')
box_34.set_patch('left', 'left')


mesh.add(box_angle1)
mesh.add(box_angle2)
mesh.add(box_angle3)
mesh.add(box_angle4)
mesh.add(box_12)
mesh.add(box_13)
mesh.add(box_23)
mesh.add(box_34)

mesh.set_default_patch('walls', 'wall')
# mesh.merge_patches('back1','mix-back1')
# mesh.merge_patches('front3','mix-front1')

mesh.write(os.path.join (".","system","solid" , "blockMeshDict"))


# %%
