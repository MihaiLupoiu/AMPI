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

PSRC  = qrm_coverage_mod.F90 qrm_coverage.F90 qrm_test_err.F90 \
	qrm_test_ord.F90 qrm_test_c.F90 qrm_test_sing.F90 \
	qrm_test_facto.F90 qrm_test_slv.F90

CPSRC = qrm_test_solve_c.c

COBJS = 

FSRC   = $(addprefix $(ARITH), $(PSRC))
CFSRC  = $(addprefix $(ARITH), $(CPSRC))

FOBJS  = $(FSRC:.F90=.o)
CFOBJS = $(CFSRC:.c=.o)

INCLUDES=-I../include -I.


dprec: 
	($(MAKE) qrm_coverage PREC=-Ddprec ARITH=d)

sprec: 
	($(MAKE) qrm_coverage PREC=-Dsprec ARITH=s)

zprec: 
	($(MAKE) qrm_coverage PREC=-Dzprec ARITH=z)

cprec: 
	($(MAKE) qrm_coverage PREC=-Dcprec ARITH=c)



qrm_coverage:   $(FOBJS) $(CFOBJS) $(COBJS)
	$(LINK) $(LDFLAGS) -o $(ARITH)$@ $(FOBJS) $(CFOBJS) $(COBJS) \
	$(QRM_LIB) \
	$(QRM_COMMONLIB) \
	$(LLAPACK) \
	$(LBLAS) \
	$(LCOLAMD) \
	$(LMETIS) \
	$(LSCOTCH) 



$(FSRC): $(ARITH)%.F90 : %.F90
	perl ../include/$(ARITH).pl $< > $@
$(CFSRC): $(ARITH)%.c : %.c
	perl ../include/$(ARITH).pl $< > $@


lib_dprec:
	(cd ..; $(MAKE) dprec)

lib_sprec:
	(cd ..; $(MAKE) sprec)

lib_zprec:
	(cd ..; $(MAKE) zprec)

lib_cprec:
	(cd ..; $(MAKE) cprec)


$(addprefix $(ARITH),qrm_coverage.o) : $(addprefix $(ARITH),qrm_coverage_mod.o qrm_test_err.o \
	qrm_test_ord.o qrm_test_c.o)


$(addprefix $(ARITH),qrm_test_err.o) : $(addprefix $(ARITH),qrm_coverage_mod.o)

$(addprefix $(ARITH),qrm_test_ord.o) : $(addprefix $(ARITH),qrm_coverage_mod.o)

$(addprefix $(ARITH),qrm_test_c.o) : $(addprefix $(ARITH),qrm_coverage_mod.o)

# qrm_test_ds.o: qrm_test_utils.o

.PHONY: clean distclean

distclean: clean
	(rm -f d*.F90 s*.F90 z*.F90 c*.F90 d*.c s*.c z*.c c*.c)

clean:
	(rm -f *.o *.mod dqrm_coverage sqrm_coverage cqrm_coverage zqrm_coverage)

