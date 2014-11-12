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

produces an X-shaped ion beam.
it is better to use BeamX.

x   = maximal radial deviation [mm]
a   = maximal radial angle [mrad]
y   = maximal axial deviation [mm]
b   = maximal axial angle [mrad]
dk  = relative energy deviation [permille]
dm  = relative mass deviation [permille]
num = number of particles [unitless]"""


AddBeam3d = """AddBeamX(xmax, amax, ymax, bmax, dk, dm, delta)

produces a nice 3d ion beam.
it is better to use Beam3d(xmax, amax, ymax, bmax, dk, dm, num).

x     = maximal radial deviation [mm]
a     = maximal radial angle [mrad]
y     = maximal axial deviation [mm]
b     = maximal axial angle [mrad]
dk    = relative energy deviation [permille]
dm    = relative mass deviation [permille]
delta = degree steps of the emittance elipse.
        smaller delta --> more particles [1..360]"""


Beam3d = """BeamX(xmax, amax, ymax, bmax, dk, dm, num)

produces a nice 3d ion beam.
it is better to use BeamX.

x   = maximal radial deviation [mm]
a   = maximal radial angle [mrad]
y   = maximal axial deviation [mm]
b   = maximal axial angle [mrad]
dk  = relative energy deviation [permille]
dm  = relative mass deviation [permille]
num = number of particles [unitless]"""


AddGaussBeam = """AddGaussBeam(strag_x, strag_a, strag_y, strag_b, x=0., a=0., y=0., b=0., dk=0., dm=0., strag_k=0., strag_m=0., num=250.)

produces a gaussian distributed ion beam.
it is better to use GaussBeam.

strag_x   = radial straggling [mm]
strag_a   = radial angular straggling [mrad]
strag_y   = axial straggling [mm]
strag_b   = axial angular straggling [mrad]
x, y      = radial/axial deviation [mm]
a, b      = radial/axial angle [mrad]
strag_k  = relative energy straggling [permille]
strag_m  = relative mass straggling [permille]
num = number of particles [unitless]"""