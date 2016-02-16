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

MDIR  = ../modules/

MODS  = $(addprefix $(MDIR),$(addprefix $(ARITH), qrm_spmat_mod.o \
	qrm_analysis_mod.o \
	qrm_mod.o)\
	qrm_common_mod.o \
	qrm_mem_mod.o \
	qrm_error_mod.o \
	qrm_string_mod.o \
	qrm_sort_mod.o \
	qrm_adata_mod.o)


PSRC  = qrm_analyse.F90 qrm_compute_graph.F90 \
	qrm_detect_singletons.F90 \
	qrm_do_ordering.F90 qrm_do_colamd.F90 qrm_elim_tree.F90 \
	qrm_rowcount.F90 qrm_rowperm.F90 qrm_do_metis.F90 qrm_ata_graph.F90 \
	qrm_do_scotch.F90 qrm_symbolic.F90 qrm_attach_singletons.F90

FSRC   = $(addprefix $(ARITH), $(PSRC))

FOBJS  = $(FSRC:.F90=.o)

COBJS  = qrm_print_tree.o qrm_postorder.o qrm_amalg_tree.o \
	qrm_compress_data.o \
	qrm_reorder_tree.o qrm_check_cperm.o

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


$(addprefix $(ARITH),qrm_analyse.o): $(addprefix $(ARITH),qrm_compute_graph.o \
	qrm_detect_singletons.o qrm_do_ordering.o qrm_elim_tree.o qrm_rowcount.o \
	qrm_rowperm.o qrm_symbolic.o qrm_attach_singletons.o) \
	$(COBJS) \
	$(MODS) 

$(addprefix $(ARITH),qrm_compute_graph.o): $(MODS)


$(addprefix $(ARITH),qrm_do_ordering.o): $(addprefix $(ARITH),qrm_do_colamd.o qrm_do_metis.o qrm_do_scotch.o)


$(addprefix $(ARITH),qrm_do_metis.o qrm_do_scotch.o): $(addprefix $(ARITH),qrm_ata_graph.o)


$(FOBJS): $(MODS)


.PHONY: clean distclean

distclean: clean
	(rm -f s*.F90 d*.F90 c*.F90 z*.F90)

clean: 
	(rm -f *.o *.mod)

