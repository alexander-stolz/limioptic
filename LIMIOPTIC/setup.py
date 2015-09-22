from py2exe.build_exe import py2exe 
from distutils.core import setup
import matplotlib
import shutil
import os
import sys

setup( 
      author="Alexander Stolz",
      console=["setup_helper.py", "_update.py"],
      py_modules=["___RUN__LIMIOPTIC_2", "ams_spicker", "beamprofile", "importsrc", "limioptic", "syntax", "_update"],
      data_files=matplotlib.get_py2exe_datafiles(),
      options={'py2exe': { 'includes' : ["matplotlib.backends.backend_tkagg", "pyqtgraph", "vtk"], 'excludes': ['_gtkagg', '_tkagg'], "bundle_files" : 3}},
      )

ext_modules = ["liblimioptic-win.so", "liblimioptic-linux.so", "libgcc_s_dw2-1.dll", "libstdc++-6.dll"]
for module in ext_modules:
	try:
		shutil.copy(module, "dist")
	except:
		pass

try:
	sys.winver
	thisos = "windows"
except:
	thisos = "linux"

cwd = os.getcwd()
os.rename(cwd + "/dist/setup_helper.exe", cwd + "/dist/___LIMIOPTIC.exe")
shutil.move(cwd + "/dist", cwd.replace("\\", "/").rsplit("/", 1)[0].rsplit("/", 1)[0] + "/precompiled for {}".format(thisos))
shutil.rmtree(cwd + "/build")
