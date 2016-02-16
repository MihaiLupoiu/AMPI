## ##############################################################################################
##
## Copyright 2012 CNRS, INPT
##  
## This file is part of qr_mumps.
##  
## qr_mumps is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as 
## published by the Free Software Foundation, either version 3 of 
## the License, or (at your option) any later version.
##  
## qr_mumps is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##  
## You can find a copy of the GNU Lesser General Public License
## in the qr_mumps/doc directory.
##
## ##############################################################################################


## -*- Makefile -*-
##
## $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
## $Author: abuttari $
## $Version: 1.1$
## $Revision: 1980 $
##

# -*- Makefile -*-
include Make.inc


dprec: 
	(cd src; $(MAKE) dprec)
	(cd include; $(MAKE) dprec)

sprec: 
	(cd src; $(MAKE) sprec)
	(cd include; $(MAKE) sprec)

zprec: 
	(cd src; $(MAKE) zprec)
	(cd include; $(MAKE) zprec)

cprec: 
	(cd src; $(MAKE) cprec)
	(cd include; $(MAKE) cprec)

dtest: 
	(cd test; $(MAKE) dprec)

stest: 
	(cd test; $(MAKE) sprec)

ztest: 
	(cd test; $(MAKE) zprec)

ctest: 
	(cd test; $(MAKE) cprec)



.PHONY: clean distclean

clean:
	(cd src; $(MAKE) clean)

distclean:
	(cd src; $(MAKE) distclean)
	(cd lib; rm -f *.a)
	(cd include; $(MAKE) distclean)


