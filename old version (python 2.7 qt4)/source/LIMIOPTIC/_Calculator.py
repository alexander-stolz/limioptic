#-------------------------------------------------------------------------------
# Name:        math helper
# Purpose:     easy peasy function calculator and plotting tool
#
# Author:      Alexander Stolz (amstolz@gmail.com)
#
# Created:     17.08.2013
# Copyright:   (c) Alexander Stolz 2013
# Licence:     free to use
#-------------------------------------------------------------------------------

import re
from pylab import *


helptext = """function calculator
----------------------------------
try:
    4 + 5
    5. / 7.
    f(x) = x**2
    plot f, -2 - 2
    g(x, y) = x * y
    g(4, 3)
    cos(pi)
    a = pi
    b = e
    a + b
-----------------------------------\n

USE exit() command to return to LIMIOPTIC.\n"""


def main():
    print helptext
    lastInput = ""
    while lastInput != "close":
        try:
            lastInput = raw_input("> ")
            re_function = re.search("^(\w+)\((.+)\)\s*=\s*(.+)", lastInput)
            re_plot     = re.search("^plot\s(.+),\s*(-*\d+)\s*-\s*(-*\d+)", lastInput)
            if re_function:
                functionCommand = "{0[0]} = lambda {0[1]}: {0[2]}".format(re_function.groups())
                # print re_function.groups(), "-->", functionCommand
                print "--> ", functionCommand
                exec(functionCommand)
            elif re_plot:
                plotCommand = "plot(linspace({0[1]}, {0[2]}, 100), map({0[0]}, linspace({0[1]}, {0[2]}, 100))); grid(); show()".format(re_plot.groups())
                # print re_plot.groups(), "-->", plotCommand
                exec(plotCommand)
            else:
                exec("__a = " + lastInput + "; print __a")
        except Exception, e:
            print e


if __name__ == '__main__':
    main()
