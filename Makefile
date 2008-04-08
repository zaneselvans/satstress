# Makefile for satstress
SHELL = /bin/sh
SSDIR = satstress

# Python modules for which we are making documentation:
EPYDOC_MODS = $(SSDIR)/satstress.py\
           $(SSDIR)/gridcalc.py\
           $(SSDIR)/__init__.py

# Source files that are included in the current release:
PUB_SRC = $(EPYDOC_MODS)\
           $(SSDIR)/physcon.py\
           $(SSDIR)/love/JohnWahr/love.f\
           test/test_nsr_diurnal.py\
           test/test_nsr_diurnal.pkl\
           input/Europa.satellite\
           input/NSR_Diurnal_exhaustive.grid

EPYDOC_OPTS = --verbose\
              --css=doc/css/satstress.css\
              --name=satstress\
              --url=http://code.google.com/p/satstress\
              --inheritance=grouped\
              --no-private

# love is built down in this directory, which has its own Makefile:
love : satstress/love/JohnWahr/love.f
	(cd satstress/love/JohnWahr; make love)

# Generate all of the HTML documentation based on the Python __doc__ strings:
htmldoc : $(EPYDOC_MODS) doc/css/satstress.css
	epydoc --html --output=doc/html $(EPYDOC_OPTS) $(EPYDOC_MODS)

# Make PDF documentation based on the Python __doc__ strings:
pdfdoc : $(EPYDOC_MODS)
	epydoc --pdf --output=doc/pdf $(EPYDOC_OPTS) $(EPYDOC_MODS)

# Check the documentation for completeness:
doccheck : $(EPYDOC_MODS)
	epydoc --check $(EPYDOC_OPTS) $(EPYDOC_MODS)

# Get rid of all the documentation:
docclean : 
	rm -f doc/html/* doc/pdf/*

# Make all forms of documentation, and check it:
doc : htmldoc pdfdoc doccheck

# See if satstress is working:
check : love $(PUB_SRC)
	python test/test_nsr_diurnal.py

# An alias for check:
test : check

# Get rid of random junk:
clean :
	rm -rf *~ *.pyc lovetmp-* *.nc
	(cd $(SSDIR); rm -rf *~ *.pyc lovetmp-* *.nc)
	rm -rf dist
	rm -rf build
	rm -f MANIFEST

# Get rid of all non-source files
realclean : clean docclean
	rm -rf $(SSDIR)/love/*/calcLove*

# Restore the package to pristene condition:
distclean : realclean htmldoc

dist : distclean $(PUB_SRC)
	python setup.py sdist

install : love $(PUB_SRC)
	python setup.py install

uninstall :
	rm -rf /usr/local/bin/calcLoveWahr4Layer
	rm -rf /Library/Python/2.5/site-packages/satstress*
