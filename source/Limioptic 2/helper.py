Drift = """Drift(length)

length in [m]

for more options use AddDrift(num, gamma2, length)"""


AddDrift = """AddDrift(num, gamma2, length)

num    = number of drifts with length 'length' [unitless]
gamma2 = relativistic factor
length = longitudinal length [m]"""


Particle = AddParticle = """[Add]Particle(x, a, y, b, dk, dm)

x  = radial deviation [mm]
a  = radial angle [mrad]
y  = axial deviation [mm]
b  = axial angle [mrad]
dk = relative energy deviation [permille]
dm = relative mass deviation [permille]"""


AddBeamX = """AddBeamX(xmax, amax, ymax, bmax, dk, dm, delta)

produces an X-shaped ion beam.
it is better to use BeamX(xmax, amax, ymax, bmax, dk, dm, num).

x     = maximal radial deviation [mm]
a     = maximal radial angle [mrad]
y     = maximal axial deviation [mm]
b     = maximal axial angle [mrad]
dk    = relative energy deviation [permille]
dm    = relative mass deviation [permille]
delta = degree steps of the emittance elipse.
        smaller delta --> more particles [1..360]"""


BeamX = """BeamX(xmax, amax, ymax, bmax, dk, dm, num)

produces an X-shaped ion beam. it is better to use BeamX.

x   = maximal radial deviation [mm]
a   = maximal radial angle [mrad]
y   = maximal axial deviation [mm]
b   = maximal axial angle [mrad]
dk  = relative energy deviation [permille]
dm  = relative mass deviation [permille]
num = number of particles [unitless]"""
