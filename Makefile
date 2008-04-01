# Makefile for SatStress
SHELL = /bin/sh
VERSION = 0.1.1
RELEASE = SatStress-$(VERSION)
SSDIR = SatStress

# Python modules for which we are making documentation:
DOC_MODS = $(SSDIR)/SatStress.py $(SSDIR)/GridCalc.py $(SSDIR)/__init__.py

# Python modules that are included in the current release:
PUB_MODS = $(DOC_MODS) $(SSDIR)/physcon.py

# The modules which are not yet being released to the public:
PRIVATE_MODS = CrossCutGraph.py Lineament.py

# All the files that need to be excised from the release:
PRIVATE_FILES = $(PRIVATE_MODS) todo Love/JamesRoberts

# love is built down in this directory, which has its own Makefile:
love : Love/JohnWahr/love

# Generate all of the HTML documentation based on the Python __doc__ strings:
htmldoc : $(DOC_MODS) doc/css/SatStress.css
	epydoc --html --verbose --output=doc/html --css doc/css/SatStress.css --name SatStress --url=http://code.google.com/p/satstress $(DOC_MODS)

# Make PDF documentation based on the Python __doc__ strings:
pdfdoc : $(DOC_MODS)
	epydoc --pdf --verbose --output=doc/pdf --name SatStress $(DOC_MODS)

# Check the documentation for completeness:
doccheck : $(DOC_MODS)
	epydoc --check --verbose $(DOC_MODS)

# Get rid of all the documentation:
docclean : 
	rm -f doc/html/* doc/pdf/*

# Make all forms of documentation, and check it:
doc : htmldoc pdfdoc doccheck

# Get rid of random junk:
clean :
	rm -rf *~ *.pyc lovetmp-* .*.swp

# Get the package back to its pristene condition:
distclean : clean
	rm -rf Love/*/love

# Get rid of all non-source files:
realclean : distclean docclean

# See if SatStress is working:
check : love $(PUB_MODS)
	python $(SSDIR)/test.py input/Europa.satellite test/SS_test_calc.pkl

# Just another alias for check:
test : check

all : love test

# Create a gzipped tarball containing all the distribution files:
release : realclean htmldoc doccheck
	cp -r ./ ../$(RELEASE)
	rm -rf ../$(RELEASE)/.svn
	rm -rf ../$(RELEASE)/*/.svn
	rm -rf ../$(RELEASE)/*/*/.svn
	(cd ../$(RELEASE); rm -rf $(PRIVATE_FILES)
	(cd ../$(RELEASE); make clean)
	(cd ..; tar --dereference -cvzf $(RELEASE).tar.gz $(RELEASE))
	rm -rf ../$(RELEASE)
