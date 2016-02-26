import os, glob, subprocess

flags = ["--hdf5-fix-global-id", "--hdf5-fix-cartesian-point",
    "--hdf5-ref-attributes", "IfcCartesianPointList3D.CoordList", "IfcDiscreteParameterValueList.Values", "IfcContinuousParameterValueList.Values", "IfcDiscreteGridOffsetList.Offsets", "IfcContiousGridOffsetList.Offsets"]

for i in range(2):
    for fn in glob.glob("generated_files/spf/*.ifc"):
        fn2 = fn.replace("/spf", "/hdf").replace(".ifc", "-compressed.hdf" if i else ".hdf")
        if not os.path.exists(fn2):
            subprocess.call(["bin/ifcconvert/IfcConvert", fn, fn2] + flags)
    flags += ["--hdf5-compress"]
