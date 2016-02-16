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

include ../../Make.inc
LIB=../../lib/lib${PREF}qrm.a
CLIB=../../lib/libqrm_common.a

OBJS  = 
COBJS = qrm_colamd_wrap.o qrm_get_wtime.o qrm_hwloc_utils.o \
	qrm_metis_wrap.o

INCLUDES=-I../../include -I.

sprec: 
	($(MAKE) objs PREC=-Dsprec)
	$(AR) $(CLIB) $(COBJS) 

dprec: 
	($(MAKE) objs PREC=-Ddprec)
	$(AR) $(CLIB) $(COBJS) 

cprec: 
	($(MAKE) objs PREC=-Dcprec)
	$(AR) $(CLIB) $(COBJS) 

zprec: 
	($(MAKE) objs PREC=-Dzprec)
	$(AR) $(CLIB) $(COBJS) 

objs:	$(OBJS) $(COBJS)


.PHONY: clean

distclean:
	(rm -f *.o )

clean: distclean

