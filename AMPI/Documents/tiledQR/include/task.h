/*
 *  task.h
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/1/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */

#ifndef HEADER_TASK
#define HEADER_TASK

#include "common.h"
#include "main_scheduler.h"

#include <msg/msg.h>
#include <xbt/set.h>

typedef enum {F,H,Z,V,ZS,VS} task_type_t;



typedef struct s_task_t{
  //XBT_SET_HEADERS;
  char name[20];
  

  
  void (*p_fptr)(task_t p_task);
  task_type_t e_type;
  m_host_t p_proc_map;
  unsigned char uc_done;
  unsigned int ui_tile_size;
  double f_unit_cost;
  union task_u {
    struct V_s {
      /*line holding the final triangle*/
      unsigned int ui_i;
      /*line zeroed out*/
      unsigned int ui_ii;
      /*step, or column onto which Z was applied*/
      unsigned int ui_j;
      /* comlumn onto which the update needs to be applied*/
      unsigned int ui_k;
    } V;
    
    struct Z_s {
      /*line holding the final triangle*/
      unsigned int ui_i;
      /*line getting zeroed out*/
      unsigned int ui_ii;
      /*step, or column onto which Z is to be applied*/
      unsigned int ui_j;
    } Z;
    
    struct H_s {
      /*line holding the triangle*/
      unsigned int ui_i;
      /*step, or column onto which F was applied*/
      unsigned int ui_j;
      /* comlumn onto which the update needs to be applied*/
      unsigned int ui_k;
    } H;
    
    struct F_s {
      /*line holding the triangle*/
      unsigned int ui_i;
      /*step, or column onto which F is to be applied*/
      unsigned int ui_j;
    } F;
  } task;
  
  xbt_dynar_t a_ancestors;
} s_task_t;


task_t new_task(task_type_t e_type, msched_args_t ps_args, task_t p_ancestor, unsigned int ui_a1, unsigned int ui_a2, unsigned int ui_a3, unsigned int ui_a4);

void delete_task(task_t p_task);

#endif
