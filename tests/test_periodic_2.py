#!/usr/bin/python3
import random

from fibergen import *

fg = FiberGen()
fg.simbox_xsize = 5.0
fg.simbox_ysize = 5.0
fg.simbox_zsize = 5.0
fg.step_num = [60, 60, 60]
fg.show_viewer = True
fg.num_threads = 1

# capsuleR1 = 0.5
# c1 = Capsule()
# c1.a = [0.0, 0.0, 0.0]
# c1.b = [0.0, 0.0, 3.2]
# c1.r = capsuleR1
# fg.capsules.append(c1)

num_fibers = 5
radius = 0.25
radius_m = 0.05

sp_jitter = 0.42
ep_jitter = 0.65

for ix in range(num_fibers):
    for iy in range(num_fibers):
        xc = random.uniform(0.0, fg.simbox_xsize)
        yc = random.uniform(0.0, fg.simbox_ysize)
        rc = random.uniform(radius - radius_m, radius + radius_m)
        lc = random.random()*0.85*fg.simbox_zsize
        zc = -2+2*random.random()
        capsule = Capsule([xc+random.random()*sp_jitter, yc+random.random()*sp_jitter, zc],
                        [xc+random.random()*ep_jitter, yc+random.random()*ep_jitter, zc +lc], rc)
        fg.capsules.append(capsule)


fg.generate_stl("periodic_2.stl")


