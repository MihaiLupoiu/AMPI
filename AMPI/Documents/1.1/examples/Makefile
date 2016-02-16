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

QRM_LIB=../lib/lib${ARITH}qrm.a

TPSRC = qrm_test.F90 
CPSRC = dqrm_test_loop_c.c 
SPSRC = qrm_test_small.F90


TFSRC   = $(addprefix $(ARITH), $(TPSRC))
CFSRC   = $(addprefix $(ARITH), $(CPSRC))
SFSRC   = $(addprefix $(ARITH), $(SPSRC))

TFOBJS  = $(TFSRC:.F90=.o)
CFOBJS  = $(CFSRC:.c=.o)
SFOBJS  = $(SFSRC:.F90=.o)
BFOBJS  = qrm_test_ds.o


INCLUDES=-I../include -I.


dprec: 
	($(MAKE) all PREC=-Ddprec ARITH=d)

sprec: 
	($(MAKE) all PREC=-Dsprec ARITH=s)

zprec: 
	($(MAKE) all PREC=-Dzprec ARITH=z)

cprec: 
	($(MAKE) all PREC=-Dcprec ARITH=c)


all: qrm_test qrm_test_small qrm_test_c


qrm_test:   $(TFOBJS)
	$(LINK) $(LDFLAGS) -o $(ARITH)$@ $(TFOBJS) \
	$(QRM_LIB) \
	$(QRM_COMMONLIB) \
	$(LLAPACK) \
	$(LBLAS) \
	$(LCOLAMD) \
	$(LMETIS) \
	$(LSCOTCH) 

qrm_test_small: $(SFOBJS)
	$(LINK) $(LDFLAGS) -o $(ARITH)$@ $(SFOBJS) \
	$(QRM_LIB) \
	$(QRM_COMMONLIB) \
	$(LLAPACK) \
	$(LBLAS) \
	$(LCOLAMD) \
	$(LMETIS) \
	$(LSCOTCH) 

qrm_test_c: $(CFOBJS)
	$(LINK)  $(LDFLAGS) -o $(ARITH)$@ $(CFOBJS) \
	$(QRM_LIB) \
	$(QRM_COMMONLIB) \
	$(LLAPACK) \
	$(LBLAS) \
	$(LCOLAMD) \
	$(LMETIS) \
	$(LSCOTCH) 

qrm_test_ds: $(BFOBJS) 
	$(LINK) $(LDFLAGS) -o $@ $(BFOBJS) \
	$(QRM_SLIB) \
	$(QRM_DLIB) \
	$(QRM_COMMONLIB) \
	$(LLAPACK) \
	$(LBLAS) \
	$(LCOLAMD) \
	$(LMETIS) \
	$(LSCOTCH) 


$(TFSRC): $(ARITH)%.F90 : %.F90
	perl ../include/$(ARITH).pl $< > $@

$(CFSRC): $(ARITH)%.c : %.c
	perl ../include/$(ARITH).pl $< > $@

$(SFSRC): $(ARITH)%.F90 : %.F90
	perl ../include/$(ARITH).pl $< > $@


lib_dprec:
	(cd ..; $(MAKE) dprec)

lib_sprec:
	(cd ..; $(MAKE) sprec)

lib_zprec:
	(cd ..; $(MAKE) zprec)

lib_cprec:
	(cd ..; $(MAKE) cprec)


.PHONY: clean distclean

distclean: clean
	(rm -f d*.F90 s*.F90 z*.F90 c*.F90  d*.c s*.c z*.c c*.c)

clean:
	(rm -f *.o *.mod dqrm_test sqrm_test cqrm_test zqrm_test \
	dqrm_test_c sqrm_test_c zqrm_test_c cqrm_test_c dqrm_test_small \
	sqrm_test_small zqrm_test_small cqrm_test_small qrm_test_ds)

