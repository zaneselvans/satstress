"""
Tools for analysing the relationship between tidal stresses and tectonics on
icy satellites.

Written by U{Zane Selvans <http://zaneselvans.org>}
(C{U{zane.selvans@colorado.edu <mailto:zane.selvans@colorado.edu>}}) as part of
his Ph.D. dissertation research.

C{SatStress} is released under GNU General Public License (GPL) version 3.  For
the full text of the license, see: U{http://www.gnu.org/}

The project is hosted at Google Code: U{http://code.google.com/p/satstress}

1 Installation
==============
Hopefully getting C{SatStress} to work on your system is a relatively painless
process, however, the software does assume you have basic experience with the
Unix shell and programming within a Unix environment (though it should work on
Windows too).  In particular, this installation information assumes you already
have and are able to use:

  - compilers for both C and Fortran.  Development has been done on Mac OS X
    (10.5) using the GNU compilers C{gcc} and C{g77}, so those should
    definitely work.  On other systems, with other compilers, your mileage may
    vary.

  - the C{make} utility, which manages dependencies between files.

1.1 Other Required and Recommended Software
-------------------------------------------
To get the L{SatStress} package working, you'll need to install some other
(free) software first:

  - B{Python 2.5} or later (U{http://www.python.org}).  If you're running a
    recent install of Linux, or Apple's Leopard operating system (OS X 10.5.x),
    you already have this.  Python is also available for Microsoft Windows, and
    just about any other platform you can think of.

  - B{SciPy} (U{http://www.scipy.org}), a collection of scientific libraries
    that extend the capabilities of the Python language.

In addition, if you want to use L{GridCalc}, you'll need:

  - B{netCDF} (U{http://www.unidata.ucar.edu/software/netcdf/}, a library of
    routines for storing, retrieving, and annotating regularly gridded
    multi-dimensional datasets.  Developed by U{Unidata
    <http://www.unidata.ucar.edu>}

  - B{netcdf4-python} (U{http://code.google.com/p/netcdf4-python/}, a Python
    interface to the netCDF library.

If you want to actually view L{GridCalc} output, you'll need a netCDF file
viewing program.  Many commercial software packages can read netCDF files, such
as ESRI ArcGIS and Matlab (from the Mathworks).  A simple and free reader for
OS X is U{Panoply <http://www.giss.nasa.gov/tools/panoply/>}, from NASA.  If
you want to really be able to interact with the outputs from this model, you
should install:

  - B{Matplotlib/Pylab} (U{http://matplotlib.sourceforge.net/}), a Matlab-like
    interactive plotting and analysis package, which uses Python as its
    "shell".

1.2 Building and Installing SatStress
-------------------------------------
Once you have the required software prerequisites installed, you should be able
to C{cd} into the directory containing the L{SatStress} module, and simply type
C{make all} at the command line.  This will compile the Love number code and
run the small L{SatStress.test} program embedded within L{SatStress}, just to
make sure that everything is in working order (or not).  If you're not using
the GNU Fortran 77 compiler C{g77}, you'll need to edit the C{Makefile} for the
Love number code::

  ./Love/JohnWahr/Makefile

and tell it what Fortran compiler it ought to be using.

2 Design Overview
=================
A few notes on the general architecture of the C{SatStress} package.

2.1 Who is the Audience?
------------------------
  In writing this software and documentation, my hope is that an undergraduate
  research assistant who has been hired for the summer, and who has at least
  some experience with programming (though not necessarily in Python), should
  be able to understand how the system works, and make fruitful use of it.  So
  if it seems like things are sometimes over-explained or over-commented,
  that's why.

2.2 A Toolkit, not a Program
----------------------------
  The C{SatStress} package is not itself a stand-alone program (or not much of
  one anyway).  Instead it is a set of tools with which you can build programs
  that need to know about the stresses on the surface of a satellite, and how
  they compare to tectonic features, so you can do your own hypothesizing and
  testing.

2.3 Object Oriented
-------------------
  The package attempts to make use of U{object oriented programming
  <http://en.wikipedia.org/wiki/Object-oriented_programming>} (OOP) in order to
  maximize the re-usability and extensibility of the code.  Many scientists are
  more familiar with the U{imperative programming style
  <http://en.wikipedia.org/wiki/Imperative_programming>} of languages like
  Fortran and C, but as more data analysis and hypothesis testing takes place
  inside computers, and as many scientists become highly specialized and
  knowledgeable software engineers (even if they don't want to admit it), the
  advantages of OOP become significant.  If the object orientation of this
  module seems odd at first glance, don't despair, it's worth learning.

2.4 Written in Python
---------------------
  U{Python <http://www.python.org>} is a general purpose, high-level scripting
  language.  It is an interpreted language (as opposed to compiled languages
  like Fortran or C) and so Python code is very portable, meaning it is usable
  on a wide variety of computing platforms without any alteration.  It is
  relatively easy to learn and easy to read, and it has a very active
  development community.  It also has a large base of friendly, helpful
  scientific users and an enormous selection of pre-existing libraries designed
  for scientific applications.  For those tasks which are particularly
  computationally intensive, Python allows you to extend the language with code
  written in C and Fortran.  Python is also U{Free Software
  <http://www.gnu.org/philosophy/free-sw.html>}.  If you are a scientist and
  you write code, Python is a great choice.

2.5 Open Source
---------------
  Because science today is intimately intertwined with computation, it is
  important for researchers to share the code that their scientific results are
  based on.  No matter how elegant and accurate your derivation is, if your
  implementation of the model in code is wrong, your results will be flawed.
  As our models and hypotheses become more complex, our code becomes vital
  primary source material, and it needs to be open to peer review.  Opening our
  source:

    - allows bugs to be found and fixed more quickly
    - facilitates collaboration and interoperability
    - reduces duplicated effort
    - enhances institutional memory
    - encourages better software design and documentation

  Of course, it also means that other people can use our code to write their
  own scientific papers, but I{that is the fundamental nature of science}.  We
  are all "standing on the shoulders of giants".  Nobody re-derives quantum
  mechanics when they just want to do a little spectroscopy.  Why should we all
  be re-writing each others code I{ad nauseam}?  Opening scientific source code
  will ultimately increase everyone's productivity.  Additionally, a great deal
  of science is funded by the public, and our code is a major product of that
  funding.  It is unethical to make it proprietary.
"""
__all__ = ["SatStress", "physcon", "GridCalc"]
