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

include ../Make.inc

OBJS =  

INCLUDES=-I../include -I.

sprec: 
	(cd C; $(MAKE) sprec)
	(cd modules; $(MAKE) sprec)
	(cd utils; $(MAKE) sprec)
	(cd analysis; $(MAKE) sprec)
	(cd factorization; $(MAKE) sprec)
	(cd solve; $(MAKE) sprec)
	(cd methods; $(MAKE) sprec)

dprec: 
	(cd C; $(MAKE) dprec)
	(cd modules; $(MAKE) dprec)
	(cd utils; $(MAKE) dprec)
	(cd analysis; $(MAKE) dprec)
	(cd factorization; $(MAKE) dprec)
	(cd solve; $(MAKE) dprec)
	(cd methods; $(MAKE) dprec)

cprec: 
	(cd C; $(MAKE) cprec)
	(cd modules; $(MAKE) cprec)
	(cd utils; $(MAKE) cprec)
	(cd analysis; $(MAKE) cprec)
	(cd factorization; $(MAKE) cprec)
	(cd solve; $(MAKE) cprec)
	(cd methods; $(MAKE) cprec)

zprec: 
	(cd C; $(MAKE) zprec)
	(cd modules; $(MAKE) zprec)
	(cd utils; $(MAKE) zprec)
	(cd analysis; $(MAKE) zprec)
	(cd factorization; $(MAKE) zprec)
	(cd solve; $(MAKE) zprec)
	(cd methods; $(MAKE) zprec)


.PHONY: clean distclean

distclean:
	(cd modules; $(MAKE) distclean)
	(cd utils; $(MAKE) distclean)
	(cd analysis; $(MAKE) distclean)
	(cd factorization; $(MAKE) distclean)
	(cd solve; $(MAKE) distclean)
	(cd C; $(MAKE) distclean)
	(cd methods; $(MAKE) distclean)
	(rm -f *.o )

clean: 
	(cd modules; $(MAKE) clean)
	(cd utils; $(MAKE) clean)
	(cd analysis; $(MAKE) clean)
	(cd factorization; $(MAKE) clean)
	(cd solve; $(MAKE) clean)
	(cd C; $(MAKE) clean)
	(cd methods; $(MAKE) clean)
	(rm -f sprec dprec)

