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
LIB  =../../lib/lib${ARITH}qrm.a
CLIB =../../lib/libqrm_common.a

MDIR  = ../modules/

MODS  = $(addprefix $(MDIR),$(addprefix $(ARITH), qrm_spmat_mod.o \
	qrm_analysis_mod.o \
	qrm_mod.o \
	qrm_fdata_mod.o \
	qrm_sdata_mod.o) \
	qrm_error_mod.o \
	qrm_mem_mod.o \
	qrm_common_mod.o \
	qrm_trace_mod.o \
	qrm_adata_mod.o \
	qrm_task_mod.o \
	qrm_sort_mod.o)

PSRC  = qrm_readmat.F90 qrm_matmul.F90 qrm_matnrm.F90 qrm_vecnrm.F90 \
	qrm_remap_pnt.F90

FSRC  = $(addprefix $(ARITH), $(PSRC))

FOBJS = $(FSRC:.F90=.o)

COBJS = qrm_prnt_array.o 


INCLUDES=-I../../include -I.

sprec: 
	($(MAKE) objs PREC=-Dsprec  ARITH=s)

dprec: 
	($(MAKE) objs PREC=-Ddprec  ARITH=d)

cprec: 
	($(MAKE) objs PREC=-Dcprec  ARITH=c)

zprec: 
	($(MAKE) objs PREC=-Dzprec  ARITH=z)


objs:	$(FOBJS) $(COBJS)
	$(AR) $(LIB) $(FOBJS) 
	$(AR) $(CLIB) $(COBJS) 


$(FSRC): $(ARITH)%.F90 : %.F90
	perl ../../include/$(ARITH).pl $< > $@



$(FOBJS): $(MODS)


.PHONY: clean

distclean: clean
	(rm -f s*.F90 d*.F90 c*.F90 z*.F90)

clean: 
	(rm -f *.o *.mod)

