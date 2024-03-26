#!/usr/bin/python3

from fibergen import *

fg = FiberGen()
fg.simbox_xsize = 5.0
fg.simbox_ysize = 5.0
fg.simbox_zsize = 5.0
fg.step_num = [40, 40, 40]

capsuleR1 = 0.5
c1 = Capsule()
c1.a = [0.0, 0.0, 0.0]
c1.b = [0.0, 0.0, 3.2]
c1.r = capsuleR1
fg.capsules.append(c1)

fg.show_viewer = True

fg.generate_stl("periodic_1.stl")


