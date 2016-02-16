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
LIB=../../lib/lib${ARITH}qrm.a
CLIB=../../lib/libqrm_common.a

PSRC =  qrm_spmat_mod.F90 qrm_mod.F90 qrm_utils_mod.F90 qrm_analysis_mod.F90 \
	qrm_factorization_mod.F90 qrm_fdata_mod.F90 qrm_sdata_mod.F90 \
	qrm_solve_mod.F90 qrm_rfpf_mod.F90 \
	qrm_c_interface.F90 qrm_methods_mod.F90

FSRC  = $(addprefix $(ARITH),$(PSRC))

FOBJS  = $(FSRC:.F90=.o)

COBJS = qrm_error_mod.o qrm_mem_mod.o \
	qrm_adata_mod.o qrm_common_mod.o qrm_trace_mod.o \
	qrm_queue_mod.o qrm_c_comm_interface.o qrm_task_mod.o \
	qrm_string_mod.o qrm_sort_mod.o qrm_const_mod.o

INCLUDES=-I../../include -I.


sprec: 
	($(MAKE) objs PREC=-Dsprec  ARITH=s)

dprec: 
	($(MAKE) objs PREC=-Ddprec  ARITH=d)

cprec: 
	($(MAKE) objs PREC=-Dcprec  ARITH=c)

zprec: 
	($(MAKE) objs PREC=-Dzprec  ARITH=z)



objs:	$(COBJS) $(FOBJS)
	(cp *.mod ../../include)
	$(AR) $(LIB) $(FOBJS) 
	$(AR) $(CLIB) $(COBJS) 




$(FSRC): $(ARITH)%.F90 : %.F90
	perl ../../include/$(ARITH).pl $< > $@


$(addprefix $(ARITH),qrm_spmat_mod.o): $(addprefix $(ARITH),qrm_fdata_mod.o qrm_sdata_mod.o) \
	qrm_adata_mod.o qrm_mem_mod.o \
	qrm_error_mod.o 


$(addprefix $(ARITH),qrm_mod.o): $(addprefix $(ARITH),qrm_spmat_mod.o \
	qrm_utils_mod.o qrm_analysis_mod.o qrm_factorization_mod.o \
	qrm_solve_mod.o qrm_methods_mod.o) \
	qrm_mem_mod.o \
	qrm_common_mod.o \
	qrm_error_mod.o \


$(addprefix $(ARITH),qrm_methods_mod.o): $(addprefix $(ARITH),qrm_spmat_mod.o \
	qrm_utils_mod.o qrm_analysis_mod.o qrm_factorization_mod.o qrm_solve_mod.o) \
	qrm_mem_mod.o \
	qrm_common_mod.o \
	qrm_error_mod.o \


qrm_queue_mod.o: qrm_mem_mod.o qrm_const_mod.o

qrm_error_mod.o: qrm_const_mod.o

qrm_common_mod.o: qrm_adata_mod.o qrm_mem_mod.o qrm_string_mod.o qrm_error_mod.o qrm_const_mod.o


qrm_adata_mod.o: qrm_mem_mod.o qrm_const_mod.o

$(addprefix $(ARITH),qrm_fdata_mod.o): qrm_mem_mod.o qrm_const_mod.o

$(addprefix $(ARITH),qrm_analysis_mod.o): $(addprefix $(ARITH),qrm_spmat_mod.o) qrm_adata_mod.o  qrm_const_mod.o

$(addprefix $(ARITH),qrm_factorization_mod.o): $(addprefix $(ARITH),qrm_fdata_mod.o qrm_spmat_mod.o)  qrm_const_mod.o

$(addprefix $(ARITH),qrm_solve_mod.o): $(addprefix $(ARITH),qrm_fdata_mod.o qrm_spmat_mod.o)  qrm_const_mod.o


qrm_mem_mod.o: qrm_error_mod.o  qrm_const_mod.o

qrm_utils_mod.o: $(addprefix $(ARITH),qrm_spmat_mod.o) qrm_mem_mod.o  qrm_const_mod.o

$(addprefix $(ARITH),qrm_spmat_mod.o): qrm_error_mod.o \
	qrm_mem_mod.o qrm_string_mod.o  qrm_const_mod.o

.PHONY: clean distclean

distclean: clean
	(rm -f s*.F90 d*.F90 c*.F90 z*.F90)

clean: 
	(rm -f *.o *.mod)

