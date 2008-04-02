#!/usr/bin/env python

from distutils.core import setup
import SatStress

setup(name         = 'SatStress',
      version      = SatStress.__version__,
      description  = 'Tools for modeling tidal stresses on icy satellites',
      author       = SatStress.__author__,
      author_email = SatStress.__contact__,
      url          = SatStress.__projecturl__,
      packages=['SatStress',],
     )
