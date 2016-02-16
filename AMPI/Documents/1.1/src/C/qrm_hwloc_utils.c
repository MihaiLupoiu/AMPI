/* ##############################################################################################
**
** Copyright 2012 CNRS, INPT
**  
** This file is part of qr_mumps.
**  
** qr_mumps is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as 
** published by the Free Software Foundation, either version 3 of 
** the License, or (at your option) any later version.
**  
** qr_mumps is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**  
** You can find a copy of the GNU Lesser General Public License
** in the qr_mumps/doc directory.
**
** ##############################################################################################*/


/*##############################################################################################*/
/** @file qrm_hwloc_utils.c
 * FIXME: add comments
 *
 * $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
 * $Author: abuttari $
 * $Version: 1.1$
 * $Revision: 1980 $
 *
 **/
/*##############################################################################################*/


#if defined(hwloc)
#include <hwloc.h>

/* Wrapper routines for hwloc */
void qrm_hwloc_bind(int id)
{
  int depth, ret;
  unsigned i, n;
  int topodepth;
  hwloc_topology_t topology;
  hwloc_cpuset_t cpuset;
  hwloc_obj_t obj;
  
  hwloc_topology_init(&topology);
  
  hwloc_topology_load(topology);
  
  obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_CORE, id); 
  
  ret = hwloc_set_cpubind(topology, obj->cpuset, HWLOC_CPUBIND_THREAD);
  if (ret) {
    printf("Couldn't bind to core %d\n", id); 
    assert(0);
  } else {
    printf("Bound to core %d\n", id); 
  }    
  
  hwloc_topology_destroy(topology);
  
  return;
}


void qrm_hwloc_info(int *ncores, int *nnodes, int *cnode)
{
  int depth, ret;
  unsigned i, n, j;
  int topodepth, numa;
  hwloc_topology_t topology;
  hwloc_cpuset_t cpuset;
  hwloc_obj_t obj, cobj;
  hwloc_obj_type_t otype;
  
  hwloc_topology_init(&topology);
  
  hwloc_topology_load(topology);

  /* get the number os cores and NUMA nodes */
  *ncores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE);
  /* printf("ncores: %d\n",*ncores); */

  *nnodes = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_NODE);
  if(*nnodes == 0){
    otype = HWLOC_OBJ_SOCKET;
    /* printf("grouping with sockets\n"); */
    *nnodes = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_SOCKET);
  } else {
    otype = HWLOC_OBJ_NODE;
    /* printf("grouping with NUMA nodes\n"); */
  }

  /* get the handle for the first NUMA node */
  obj = hwloc_get_obj_by_type(topology, otype, 0); 
  
  /* get the number of cores in one NUMA node (supposedly the same for all nodes) */
  *cnode = hwloc_get_nbobjs_inside_cpuset_by_type(topology, obj->cpuset, HWLOC_OBJ_CORE);
  
  hwloc_topology_destroy(topology);
  return;
}


void qrm_hwloc_topo(int *nodes, int *topo)
{
  int depth, ret;
  unsigned i, n, j, ncores, nnodes, cnode;
  int topodepth, numa;
  hwloc_topology_t topology;
  hwloc_cpuset_t cpuset;
  hwloc_obj_t obj, cobj;
  hwloc_obj_type_t otype;
  
  hwloc_topology_init(&topology);
  
  hwloc_topology_load(topology);

  /* get the number os cores and NUMA nodes */
  ncores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE);
  printf("ncores: %d\n",ncores);

  nnodes = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_NODE);
  if(nnodes == 0){
    otype = HWLOC_OBJ_SOCKET;
    printf("grouping with sockets\n");
    nnodes = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_SOCKET);
  } else {
    otype = HWLOC_OBJ_NODE;
    printf("grouping with NUMA nodes\n");
  }

  /* get the handle for the first NUMA node */
  obj = hwloc_get_obj_by_type(topology, otype, 0); 
  
  /* get the number of cores in one NUMA node (supposedly the same for all nodes) */
  cnode = hwloc_get_nbobjs_inside_cpuset_by_type(topology, obj->cpuset, HWLOC_OBJ_CORE);
  
  for(i=0; i<nnodes; i++){
    /* get the handle for the first i-th node */
    obj = hwloc_get_obj_by_type(topology, otype, i);
    /* get the number of cores in i-th NUMA node (supposedly the same for all nodes) */
    cnode = hwloc_get_nbobjs_inside_cpuset_by_type(topology, obj->cpuset, HWLOC_OBJ_CORE);

    /* get the first core in this node */
    cobj = hwloc_get_next_obj_inside_cpuset_by_type(topology, obj->cpuset, HWLOC_OBJ_CORE, NULL);
    topo[(i*cnode)] = cobj->logical_index;
    /* printf("%2d -- group:  %2d",i,cobj->logical_index); */
    for(j=1; j<cnode; j++){
      cobj = hwloc_get_next_obj_inside_cpuset_by_type(topology, obj->cpuset, HWLOC_OBJ_CORE, cobj);
      topo[(i*cnode)+j] = cobj->logical_index;
      /* printf(" %2d",cobj->logical_index); */
    }
    /* printf("\n"); */
  }
  
  
  hwloc_topology_destroy(topology);
  return;
}



























/* void qrm_hwloc_topo(int *topo) */
/* { */
/*   int depth, ret; */
/*   unsigned i, n, cnode, j, ncores, nnodes; */
/*   int topodepth, dryrun; */
/*   hwloc_topology_t topology; */
/*   hwloc_cpuset_t cpuset; */
/*   hwloc_obj_t obj, cobj; */
  
/*   hwloc_topology_init(&topology); */
/*   hwloc_topology_load(topology); */
  
/*   /\* get the number os cores and NUMA nodes *\/ */
/*   ncores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE); */
/*   nnodes = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_NODE); */
  
/*   for(i=0; i<nnodes; i++){ */
/*     /\* get the handle for the first i-th node *\/ */
/*     obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_NODE, i); */
/*     /\* get the number of cores in i-th NUMA node (supposedly the same for all nodes) *\/ */
/*     cnode = hwloc_get_nbobjs_inside_cpuset_by_type(topology, obj->cpuset, HWLOC_OBJ_CORE); */

/*     /\* get the first core in this node *\/ */
/*     cobj = hwloc_get_next_obj_inside_cpuset_by_type(topology, obj->cpuset, HWLOC_OBJ_CORE, NULL); */
/*     topo[(i*cnode)] = cobj->logical_index; */
/*     for(j=1; j<cnode; j++){ */
/*       cobj = hwloc_get_next_obj_inside_cpuset_by_type(topology, obj->cpuset, HWLOC_OBJ_CORE, cobj); */
/*       topo[(i*cnode)+j] = cobj->logical_index; */
/*     } */
/*   } */

/*   hwloc_topology_destroy(topology); */
/*   return; */
/* } */


/* void qrm_hwloc_info(int *ncores, int *nnodes, int *cnode) */
/* { */
/*   int depth, ret; */
/*   unsigned i, n, icnode, j; */
/*   int topodepth, dryrun; */
/*   hwloc_topology_t topology; */
/*   hwloc_cpuset_t cpuset; */
/*   hwloc_obj_t obj, cobj; */
  
/*   hwloc_topology_init(&topology); */
  
/*   hwloc_topology_load(topology); */
  
/*   /\* get the number os cores and NUMA nodes *\/ */
/*   *ncores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE); */
/*   *nnodes = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_NODE); */

/*   if(*nnodes == 0){ */
/*     *nnodes = 1; */
/*     *cnode  = *ncores; */
/*     hwloc_topology_destroy(topology); */
/*     return; */
/*   } */
/*   /\* get the handle for the first NUMA node *\/ */
/*   obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_NODE, 0);  */
  
/*   /\* get the number of cores in one NUMA node (supposedly the same for all nodes) *\/ */
/*   *cnode = hwloc_get_nbobjs_inside_cpuset_by_type(topology, obj->cpuset, HWLOC_OBJ_CORE); */

/*   hwloc_topology_destroy(topology); */
/*   return; */
/* } */


#endif
