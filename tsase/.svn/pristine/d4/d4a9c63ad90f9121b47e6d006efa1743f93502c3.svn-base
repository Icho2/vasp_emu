#!/usr/bin/env python
from numpy.distutils.core import setup, Extension
import subprocess
from sys import exit
import os


packages = []
for dirname, dirnames, filenames in os.walk('tsase'):
    if '__init__.py' in filenames:
        packages.append(dirname.replace('/', '.'))


al_ext = Extension(
    name='tsase.calculators.al.al_',
    sources=[
        'tsase/calculators/al/code/sumembforce.f',
        'tsase/calculators/al/code/potinit.f',
        'tsase/calculators/al/code/gagafeDblexp.f',
        'tsase/calculators/al/code/forces.f',
        'tsase/calculators/al/code/fofrhoDblexp.f',
        'tsase/calculators/al/code/embedenergy.f',
        'tsase/calculators/al/code/dfrhoDblexp.f',
    ],
    include_dirs=['tsase/calculators/al/code/commonblks']
)


lepspho_ext = Extension(
    name='tsase.calculators.lepspho.lepspho_',
    sources=[
        'tsase/calculators/lepspho/code/lepspho.f',
    ],
)


lj_ext = Extension(
    name='tsase.calculators.lj.lj_',
    sources=[
        'tsase/calculators/lj/code/lj.f',
    ],
)


morse_ext = Extension(
    name='tsase.calculators.morse.morse_',
    sources=[
        'tsase/calculators/morse/code/morse.f',
    ],
)


package_dir = {'tsase': 'tsase'}

scripts = [
    "bin/dump2xdat",
    "bin/kmc",
    "bin/lmp2con",
    "bin/lmp2pos",
    "bin/mobfil",
    "bin/neighbors",
    "bin/oldexpectra",
    "bin/pdf-make.py",
    "bin/soc2pos",
    "bin/splitxdat",
    "bin/temgui",
    "bin/tsase",
    "bin/water_solvate",
    "bin/water_solvate_z",
    "bin/xyz",
]

package_data = {'tsase': ['xyz/xyz.glade',
                          'xyz/xyz.help',
                          'xyz/*.png',
                          'calculators/al/al_.so',
                          'calculators/cuo/ffield.comb',
                          'calculators/cuo/ffield.comb3',
                          'calculators/cuo/in.lammps',
                          'calculators/lepspho/lepspho_.so',
                          'calculators/lisi/LiSi.meam',
                          'calculators/lisi/in.lammps',
                          'calculators/lisi/library.meam',
                          'calculators/lj/lj_.so',
                          'calculators/mo/Mo.set',
                          'calculators/mo/in.lammps',
                          'calculators/morse/morse_.so',
                          'calculators/si/Si.meam',
                          'calculators/si/library.meam',
                          'calculators/w/W.set',
                          'calculators/w/in.lammps']}

setup(
    name='tsase',
    version='1.0',
    description='Library based upon ASE for transition state theory calculations.',
    author='Henkelman Research Group',
    author_email='henkelman@utexas.edu',
    url='http://www.henkelmanlab.org',
    #      packages=['tsase'],
    scripts=scripts,
    packages=packages,
    package_data=package_data,
    ext_modules=[al_ext, lepspho_ext, lj_ext, morse_ext]
)
