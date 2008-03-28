# Makefile for SatStress

PUB_MODS = SatStress.py GridCalc.py __init__.py

love : Love/JohnWahr/love

# Generate all of the HTML documentation:
htmldoc : $(PUB_MODS) doc/css/SatStress.css
	epydoc --html --no-frames --verbose --output=doc/html --css doc/css/SatStress.css --name SatStress --url=http://code.google.com/p/satstress $(PUB_MODS)

pdfdoc : $(PUB_MODS)
	epydoc --pdf --verbose --output=doc/pdf --name SatStress $(PUB_MODS)

checkdoc : $(PUB_MODS)
	epydoc --check --verbose $(PUB_MODS)

cleandoc : 
	rm -f doc/html/* doc/pdf/*

doc : htmldoc pdfdoc checkdoc

clean :
	rm -rf *~ *.pyc lovetmp-*

realclean : clean cleandoc
	rm -rf Love/*/love

test : love $(PUB_MODS)
	python SatStress.py input/Europa.satellite test/SS_test_calc.pkl

all : love test

