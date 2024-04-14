#!/usr/bin/python3

from fibergen import *

fg = FiberGen()
fg.num_threads = 4
fg.simbox_xsize = 5.0
fg.simbox_ysize = 5.0
fg.simbox_zsize = 5.0
fg.step_num = [80, 80, 80]

c1 = Capsule()
c1.a = [2.5, 2.5, -2.2]
c1.b = [2.5, 2.5, 2.2]
c1.r = 0.35

c2 = Capsule()
c2.a = [2.5, -2.2, 2.5]
c2.b = [2.5,  2.2, 2.5]
c2.r = 0.35

c3 = Capsule()
c3.a = [-2.2, 2.5, 2.5]
c3.b = [2.2, 2.5, 2.5]
c3.r = 0.35

fg.capsules.append(c1)
fg.capsules.append(c2)
fg.capsules.append(c3)
fg.mc_field_value = 1.0
fg.fc_normal_weight = 0.0
fg.fc_mutate_weight = 1.0
fg.show_viewer = True
fg.field_power = 6

fg.generate_stl("cross12.stl")


