#!/usr/bin/env python

from distutils.core import setup
import SatStress

setup(name         = 'SatStress',
      version      = SatStress.__version__,
      description  = 'Tools for modeling tidal stresses on icy satellites',
      author       = SatStress.__author__,
      author_email = SatStress.__contact__,
      url          = SatStress.__projecturl__,
      packages     = ['SatStress',],
      package_dir  = {'SatStress': 'SatStress'},
      scripts      = ['SatStress/Love/JohnWahr/calcLoveWahr4Layer',],
      package_data = {'SatStress': ['input/*.satellite', 'input/*.grid',]},
     )
      #package_data ={'SatStress': ['Love/JohnWahr/love']},
      #data_files   = [('SatStress/Love/JohnWahr',['love',]),]
