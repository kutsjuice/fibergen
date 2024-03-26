#!/usr/bin/python3

from fibergen import *

fg = FiberGen()
fg.simbox_xsize = 5.0
fg.simbox_ysize = 5.0
fg.simbox_zsize = 5.0
fg.step_num = [40, 40, 40]

c1 = Capsule()
c1.a = [2.5, 2.5, -2.2]
c1.b = [2.5, 2.5, 2.2]
c1.r = 0.5

c2 = Capsule()
c2.a = [2.5, -2.2, 2.5]
c2.b = [2.5,  2.2, 2.5]
c2.r = 0.5

c3 = Capsule()
c3.a = [-2.2, 2.5, 2.5]
c3.b = [2.2, 2.5, 2.5]
c3.r = 0.5

fg.capsules.append(c1)
fg.capsules.append(c2)
fg.capsules.append(c3)

fg.show_viewer = True

fg.generate_stl("cross1.stl")


