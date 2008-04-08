#!/usr/bin/env python

from distutils.core import setup
import satstress

# patch distutils if it can't cope with the "classifiers" or
# "download_url" keywords
from sys import version, exit
if version < satstress.__pythonrequiredversion__:
    exit("""
ERROR: satstress requires Python %s or greater!
You are currently using Python %s
    """ % (satstress.__pythonrequiredversion__, version,) )

setup(packages     = ['satstress',],
      package_dir  = {'satstress': 'satstress'},
      scripts      = ['satstress/Love/JohnWahr/calcLoveWahr4Layer',],
      package_data = {'satstress': ['input/*.satellite', 'input/*.grid',]},
      name             = satstress.__name__,
      version          = satstress.__version__,
      description      = satstress.__description__,
      long_description = satstress.__long_description__,
      license          = satstress.__license__,
      url              = satstress.__projecturl__,
      download_url     = satstress.__downloadurl__,
      author           = satstress.__author__,
      author_email     = satstress.__contact__,

      requires    = [
          "scipy",
          "numpy",
          "netCDF3",
          ],

      obsoletes   = [
          "%s (<%s)" % (satstress.__name__, satstress.__version__),
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
