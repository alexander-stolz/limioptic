from time import localtime

open("version", "w").write("{0:0>4}{1:0>2}{2:0>2}".format(*localtime()))
print ">%s<" % open("version", "r").read()
