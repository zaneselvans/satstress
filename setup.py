#!/usr/bin/env python

from distutils.core import setup
import SatStress

# patch distutils if it can't cope with the "classifiers" or
# "download_url" keywords
from sys import version, exit
if version < SatStress.__pythonrequiredversion__:
    exit("""
ERROR: SatStress requires Python %s or greater!
You are currently using Python %s
    """ % (SatStress.__pythonrequiredversion__, version,) )

setup(packages     = ['SatStress',],
      package_dir  = {'SatStress': 'SatStress'},
      scripts      = ['SatStress/Love/JohnWahr/calcLoveWahr4Layer',],
      package_data = {'SatStress': ['input/*.satellite', 'input/*.grid',]},
      name             = SatStress.__name__,
      version          = SatStress.__version__,
      description      = SatStress.__description__,
      long_description = SatStress.__long_description__,
      license          = SatStress.__license__,
      url              = SatStress.__projecturl__,
      download_url     = SatStress.__downloadurl__,
      author           = SatStress.__author__,
      author_email     = SatStress.__contact__,

      requires    = [
          "scipy",
          "numpy",
          "netCDF3",
          ],

      obsoletes   = [
          "%s (<%s)" % (SatStress.__name__, SatStress.__version__),
          ],

      classifiers = [
          'Development Status :: 2 - Pre-Alpha',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License (GPL)',
          'Natural Language :: English',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: POSIX',
          'Operating System :: Unix',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: Physics',
          'Topic :: Scientific/Engineering :: GIS',
          'Topic :: Software Development :: Libraries :: Python Modules',
          ],
     )
