import sys, numpy, math
from scipy.optimize import leastsq

# Murnaghan equation of state
def eos_murnaghan(params, vol):
    'From Phys. Rev. B 28, 5480 (1983)'
    E0, B0, Bp, V0 = params 
    E = E0 + B0/Bp * vol * ((V0/vol)**Bp/(Bp-1.0)+1.0) - V0*B0/(Bp-1.0)
    return E


# Customized input with default and accepted values
def myinput(prompt, default, accepted):
    while True:
        res = input(prompt + " [default=%s]: " % (default))
        if res == '': res = default
        if res in accepted:
            break
        else:
            print("accepted values:", accepted)
    return res


print("Welcome to Murnaghan-fit.py")
print()
fname = input("Filename containing energy vs volume [volume.dat]: ")
if fname == '': fname = 'volume.dat'

try:
    f = open(fname, 'rt')
except IOError:
    sys.stderr.write("Error opening or reading file %s\n" % (fname))
    sys.exit(1)

# read data from file
print()
print("Data read from file:")
vol = []
ene = []
while True:
    line = f.readline().strip()
    if line == '': break
    if line[0] == '#' or line[0] == '!': continue
    v, e = [float(x) for x in line.split()[:2]]
    vol.append(v)
    ene.append(e)
    print(v, e)
print()
f.close()

# transform to numpy arrays
vol = numpy.array(vol)
ene = numpy.array(ene)

# further questions
is_volume = True
ans = myinput("Volume or lattice spacing (vol/latt)", "vol", ["vol", "latt"])
is_volume = (ans == 'vol')

if is_volume:
    vol_unit = myinput("Volume units (ang3/bohr3)", "ang3", ["ang3", "bohr3"])    
    if vol_unit == "bohr3": vol *= 0.14818471    # convert to angstrom^3
else:
    vol_unit = myinput("Lattice units (ang/bohr)", "ang", ["ang", "bohr"])    
    if vol_unit == "bohr": vol *= 0.52917721     # convert to angstrom
    latt = myinput("Lattice type (sc/fcc/bcc)", "sc", ["sc", "bcc", "fcc"])
    fact = 1.0
    if latt == "fcc": fact = 0.25
    if latt == "bcc": fact = 0.5
    # lattice to volume
    vol = fact * vol**3

ene_unit = myinput("Energy units (eV/Ry/Ha)", "eV", ["eV", "Ry", "Ha"])
if ene_unit == "Ry": ene *= 13.605692            # convert to eV
if ene_unit == "Ha": ene *= 27.211383            # convert to eV
print()

# fit a parabola to the data and get inital guess for equilibirum volume
# and bulk modulus
a, b, c = numpy.polyfit(vol, ene, 2)
V0 = -b/(2*a)
E0 = a*V0**2 + b*V0 + c
B0 = 2*a*V0
Bp = 4.0

# initial guesses in the same order used in the Murnaghan function
x0 = [E0, B0, Bp, V0]

def print_params(label, params):
    E0, B0, Bp, V0 = params
    print(label, ": E0 = %f eV" % (E0))
    print(label, ": B0 = %f GPa" % (B0*160.21765))
    print(label, ": Bp = %f" % (Bp))
    if is_volume:
        print(label, ": V0 = %f angstrom^3" % (V0))
    else:
        print(label, ": V0 = %f angstrom^3, a0 = %f angstrom" % (V0, (V0/fact)**(1./3.)))
    print()

# fit the equations of state
target = lambda params, y, x: y - eos_murnaghan(params, x)
murn, ier = leastsq(target, x0, args=(ene,vol))
print_params("Murnaghan", murn)


try:
    import pylab
except ImportError:
    sys.stderr.write("pylab module non available, skipping plot")
    sys.exit(0)

# plotting
ans = myinput("Do you want to plot the result (yes/no)", "yes", ["yes", "no"])
if ans == "no": sys.exit(0)

import pylab
vfit = numpy.linspace(min(vol),max(vol),100)

pylab.plot(vol, ene-E0, 'ro')
pylab.plot(vfit, eos_murnaghan(murn,vfit)-E0, label='Murnaghan')
pylab.xlabel('Volume ($\AA^3$)')
pylab.ylabel('Energy (eV)')
pylab.legend(loc='best')
pylab.show()
quit()
