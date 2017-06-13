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

xmax   = maximal radial deviation [mm]
amax   = maximal radial angle [mrad]
ymax   = maximal axial deviation [mm]
bmax   = maximal axial angle [mrad]
dk     = relative energy deviation [permille]
dm     = relative mass deviation [permille]
delta  = degree steps of the emittance elipse.
         smaller delta --> more particles [1..360]"""


BeamX = """BeamX(xmax, amax, ymax, bmax, dk, dm, num)

produces an X-shaped ion beam.
it is better to use BeamX.

xmax = maximal radial deviation [mm]
amax = maximal radial angle [mrad]
ymax = maximal axial deviation [mm]
bmax = maximal axial angle [mrad]
dk   = relative energy deviation [permille]
dm   = relative mass deviation [permille]
num  = number of particles [unitless]"""


AddBeam3d = """AddBeamX(xmax, amax, ymax, bmax, dk, dm, delta)

produces a nice 3d ion beam.
it is better to use Beam3d(xmax, amax, ymax, bmax, dk, dm, num).

xmax   = maximal radial deviation [mm]
amax   = maximal radial angle [mrad]
ymax   = maximal axial deviation [mm]
bmax   = maximal axial angle [mrad]
dk     = relative energy deviation [permille]
dm     = relative mass deviation [permille]
delta  = degree steps of the emittance elipse.
        smaller delta --> more particles [1..360]"""


Beam3d = """BeamX(xmax, amax, ymax, bmax, dk, dm, num)

produces a nice 3d ion beam.
it is better to use BeamX.

xmax = maximal radial deviation [mm]
amax = maximal radial angle [mrad]
ymax = maximal axial deviation [mm]
bmax = maximal axial angle [mrad]
dk   = relative energy deviation [permille]
dm   = relative mass deviation [permille]
num  = number of particles [unitless]"""


AddGaussBeam = """AddGaussBeam(strag_x, strag_a, strag_y, strag_b, x=0., a=0., y=0., b=0., dk=0., dm=0., strag_k=0., strag_m=0., num=250.)

produces a gaussian distributed ion beam.
it is better to use GaussBeam.

strag_x  = radial straggling [mm]
strag_a  = radial angular straggling [mrad]
strag_y  = axial straggling [mm]
strag_b  = axial angular straggling [mrad]
x, y     = radial/axial deviation [mm]
a, b     = radial/axial angle [mrad]
strag_k  = relative energy straggling [permille]
strag_m  = relative mass straggling [permille]
num      = number of particles [unitless]"""


GaussBeam = """GaussBeam(strag_x, strag_a, strag_y, strag_b, x=0., a=0., y=0., b=0., dk=0., dm=0., num=250, strag_k=0., strag_m=0., sigma=1.)

produces a gaussian distributed ion beam.
it is better to use GaussBeam.

strag_x  = radial straggling [mm]
strag_a  = radial angular straggling [mrad]
strag_y  = axial straggling [mm]
strag_b  = axial angular straggling [mrad]
x, y     = radial/axial deviation [mm]
a, b     = radial/axial angle [mrad]
strag_k  = relative energy straggling [permille]
strag_m  = relative mass straggling [permille]
num      = number of particles [unitless]
sigma    = how many standard deviations is the straggling"""


BeamRandomGauss = AddBeamRandomGauss = """[Add]BeamRandomGauss(xmax, amax, ymax, bmax, dk, dm, num, sigma=1.)

like GaussBeam, but the beam is different for every calculation.

xmax   = maximal radial deviation [mm]
amax   = maximal radial angle [mrad]
ymax   = maximal axial deviation [mm]
bmax   = maximal axial angle [mrad]
dk     = relative energy deviation [permille]
dm     = relative mass deviation [permille]
num    = number of particles"""


AddBeam = """AddBeam(xmax, amax, ymax, bmax, dk, dm, delta=10)

produces a simple beam for fast calculations.
better use Beam(xmax, amax, ymax, bmax, dk, dm, num=250, x=0, y=0)

xmax   = maximal radial deviation [mm]
amax   = maximal radial angle [mrad]
ymax   = maximal axial deviation [mm]
bmax   = maximal axial angle [mrad]
dk     = relative energy deviation [permille]
dm     = relative mass deviation [permille]
delta  = degree steps of the emittance elipse.
         smaller delta --> more particles [1..360]"""


Beam = """Beam(xmax, amax, ymax, bmax, dk, dm, num=250, x=0, y=0)

produces a simple beam for fast calculations.
better use Beam(xmax, amax, ymax, bmax, dk, dm, num=250, x=0, y=0)

xmax   = maximal radial deviation [mm]
amax   = maximal radial angle [mrad]
ymax   = maximal axial deviation [mm]
bmax   = maximal axial angle [mrad]
dk     = relative energy deviation [permille]
dm     = relative mass deviation [permille]
num    = number of particles,
x, y   = radial/axial deviation"""


Source = AddSource = """[Add]Source()

use particle-file from srim (TRANSMIT.txt) or in vector-format."""


Matrix = AddMatrix = """AddMatrix(num, mat, length)

add your own matrix num times. length is the length of the corresponding element."""


AMSAcc = AddAMSAcc = """AddAMSAcc(v_qsnout, v_terminal, v_injection, q)

add CologneAMS accelerator.

v_qsnout    = q-snount voltage [V]
v_terminal  = terminal voltage [V]
v_injection = injection voltage [V]
q           = ion charge state [e]"""


VBFN = AddVBFN = """AddVBFN(extraktion, deltaV, laenge=.276, b=1.13, b1=-1., b2=-1., segment=0)

add FN-preacceleration.

extraktion = extraction energy [keV]
deltaV     = preacceleration voltage [v]
laenge     = length of the preacceleration [m]
b          = correction factor for input and output aperture if b1 and b2 are not set
b1, b2     = correction factors for input and output
segment    = just used internally"""


MSA = """MSA(r, alpha, geo=30., solid=True, korrektur=None, B_ist=1., B_soll=1.)

add Homogeneous Magnet (bends in X direction)

r     = bending radius
alpha = bending angle [deg]
geo   = pole shoe distance from symmetry axis [mm]
solid = does it cut the beam? [True, False]

the rest is used internally."""


MSA_Y = """MSA_Y(r, alpha, geo=30., solid=True, korrektur=None, B_ist=1., B_soll=1.)

add Homogeneous Magnet (bends in Y direction)

r     = bending radius
alpha = bending angle [deg]
geo   = pole shoe distance from symmetry axis [mm]
solid = does it cut the beam? [True, False]

the rest is used internally."""
