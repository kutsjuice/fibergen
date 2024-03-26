import random, json, os
class Capsule:
    def __init__(self, ca=None, cb=None, cr=None):
        self.a = [0.0, 0.0, 0.0] if ca is None else ca
        self.b = [0.0, 0.0, 0.0] if cb is None else cb
        self.r = 0.2 if cr is None else cr

class FiberGen:
    def __init__(self):
        self.output_mesh_name = "out.stl"
        self.fibergen_path = "../bin/fibergen"
        self.simbox_xsize = 4.0
        self.simbox_ysize = 4.0
        self.simbox_zsize = 4.0
        self.mc_field_value = 1.1
        self.num_additional_steps = 5
        self.periodic = True
        self.show_viewer = True
        self.step_num = [50, 50, 50]
        self.show_axes = False
        self.num_threads = 1
        self.capsules = []

    def append_capsule(self, ca, cb, cr):
        self.capsules.append(Capsule(ca, cb, cr))

    def generate_stl(self, stlname = None):
        output_mesh_name = self.output_mesh_name if stlname is None else stlname
        tmpfilename = "tmpfile{}.json".format(random.randint(0, 100000))
        tmpfile = open(tmpfilename, 'w')
        print("tmpfilename is {}".format(tmpfilename))
        jdict = dict()

        jdict["simbox_xsize"] = self.simbox_xsize
        jdict["simbox_ysize"] = self.simbox_ysize
        jdict["simbox_zsize"] = self.simbox_zsize
        jdict["mc_field_value"] = self.mc_field_value
        jdict["num_additional_steps"] = self.num_additional_steps
        jdict["periodic"] = self.periodic
        jdict["output_mesh_name"] = self.output_mesh_name
        jdict["show_viewer"] = self.show_viewer
        jdict["step_num"] = self.step_num
        jdict["show_axes"] = self.show_axes
        jdict["num_threads"] = self.num_threads
        jdict["capsules"] = []

        for capsule in self.capsules:
            capsule_inst = dict()
            capsule_inst["start"] = capsule.a
            capsule_inst["end"] = capsule.b
            capsule_inst["radius"] = capsule.r
            jdict["capsules"].append(capsule_inst)

        json.dump(jdict, tmpfile, indent=4)
        tmpfile.close()
        command = "{0} {1}; rm ./{1}".format(self.fibergen_path, tmpfilename)
        print(command)
        os.system(command)

        pass

if __name__ == "__main__":
    pass