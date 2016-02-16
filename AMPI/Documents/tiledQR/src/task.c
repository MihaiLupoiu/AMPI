/*
 *  task.c
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/1/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */

#include "task.h"

#include "xbt/sysdep.h"
#include "xbt/dynar.h"

extern int operation_cnt;


task_t new_task(task_type_t e_type, msched_args_t ps_args, task_t p_ancestor, unsigned int ui_a1, unsigned int ui_a2, unsigned int ui_a3, unsigned int ui_a4)
{
  task_t p_task = (task_t)xbt_new(s_task_t,1);
  p_task->uc_done = 0;
  p_task->p_proc_map = NULL;
  p_task->e_type = e_type;
  p_task->ui_tile_size = ps_args->ui_tile_size;
  p_task->a_ancestors = xbt_dynar_new(sizeof(task_t),NULL);
  
  switch (e_type) {
    case F:
      p_task->task.F.ui_i = ui_a1;
      p_task->task.F.ui_j = ui_a2;
      if (ps_args->uc_coarse) {
        p_task->f_unit_cost = 0;  
      }
      else {
        p_task->f_unit_cost = 4;        
      }

      sprintf(p_task->name,"F(%u,%u)",ui_a1,ui_a2);
      xbt_fifo_push(ps_args->a_tile_updates[p_task->task.F.ui_i*ps_args->ui_q+p_task->task.F.ui_j], p_task);
      break;
    case H:
      p_task->task.H.ui_i = ui_a1;
      p_task->task.H.ui_j = ui_a2;
      p_task->task.H.ui_k = ui_a3;
      if (ps_args->uc_coarse) {
      p_task->f_unit_cost = 0;
      }else{
      p_task->f_unit_cost = 6;
      }    
      sprintf(p_task->name,"H(%u,%u,%u)",ui_a1,ui_a2,ui_a3);
      /*place the update in the tile fifo*/
      xbt_fifo_push(ps_args->a_tile_updates[p_task->task.H.ui_i*ps_args->ui_q+p_task->task.H.ui_k], p_task);      
      xbt_dynar_push_as(p_task->a_ancestors,task_t,p_ancestor);
      break;
    case Z:
      p_task->task.Z.ui_i = ui_a1;
      p_task->task.Z.ui_ii = ui_a2;
      p_task->task.Z.ui_j = ui_a3;
      if (ps_args->uc_coarse) {
      p_task->f_unit_cost = 1;
      }else{
      p_task->f_unit_cost = 2;
      }
      xbt_fifo_push(ps_args->a_tile_updates[p_task->task.Z.ui_i*ps_args->ui_q+p_task->task.Z.ui_j], p_task);
      xbt_fifo_push(ps_args->a_tile_updates[p_task->task.Z.ui_ii*ps_args->ui_q+p_task->task.Z.ui_j], p_task);
      sprintf(p_task->name,"Z(%u,%u,%u)",ui_a1,ui_a2,ui_a3);
      break;
    case ZS:
      p_task->task.Z.ui_i = ui_a1;
      p_task->task.Z.ui_ii = ui_a2;
      p_task->task.Z.ui_j = ui_a3;
      if (ps_args->uc_coarse) {
      p_task->f_unit_cost = 1;
      }else{
      p_task->f_unit_cost = 6;
      }
      /*place the update in the tile fifo*/
      xbt_fifo_push(ps_args->a_tile_updates[p_task->task.Z.ui_i*ps_args->ui_q+p_task->task.Z.ui_j], p_task);
      xbt_fifo_push(ps_args->a_tile_updates[p_task->task.Z.ui_ii*ps_args->ui_q+p_task->task.Z.ui_j], p_task);
      sprintf(p_task->name,"ZS(%u,%u,%u)",ui_a1,ui_a2,ui_a3);      
      break;      
    case V:
      p_task->task.V.ui_i = ui_a1;
      p_task->task.V.ui_ii = ui_a2;
      p_task->task.V.ui_j = ui_a3;
      p_task->task.V.ui_k = ui_a4;      
      if (ps_args->uc_coarse) {
        if(ps_args->fptr_handle_finish!=NULL){
          p_task->f_unit_cost = 0;
        }
        else{
          p_task->f_unit_cost = 1;
        }
      }else{
      p_task->f_unit_cost = 6;
      }
      xbt_dynar_push_as(p_task->a_ancestors,task_t,p_ancestor);
      /*place the update in the tile fifo*/
      xbt_fifo_push(ps_args->a_tile_updates[p_task->task.V.ui_i*ps_args->ui_q+p_task->task.V.ui_k], p_task);
      xbt_fifo_push(ps_args->a_tile_updates[p_task->task.V.ui_ii*ps_args->ui_q+p_task->task.V.ui_k], p_task);
      sprintf(p_task->name,"V(%u,%u,%u,%u)",ui_a1,ui_a2,ui_a3,ui_a4);
      break;
    case VS:
      p_task->task.V.ui_i = ui_a1;
      p_task->task.V.ui_ii = ui_a2;
      p_task->task.V.ui_j = ui_a3;
      p_task->task.V.ui_k = ui_a4;     
      if (ps_args->uc_coarse) {
        if(ps_args->fptr_handle_finish!=NULL){
          p_task->f_unit_cost = 0;
        }
        else{
          p_task->f_unit_cost = 1;
        }
      }else{
      p_task->f_unit_cost = 12;
      }
      xbt_dynar_push_as(p_task->a_ancestors,task_t,p_ancestor);
      /*place the update in the tile fifo*/
      xbt_fifo_push(ps_args->a_tile_updates[p_task->task.V.ui_i*ps_args->ui_q+p_task->task.V.ui_k], p_task);
      xbt_fifo_push(ps_args->a_tile_updates[p_task->task.V.ui_ii*ps_args->ui_q+p_task->task.V.ui_k], p_task);
      sprintf(p_task->name,"VS(%u,%u,%u,%u)",ui_a1,ui_a2,ui_a3,ui_a4);
      break;      
    default:
      exit(-1);
      break;
  }
  
  operation_cnt+=ceil(p_task->f_unit_cost);
  xbt_dynar_push_as(ps_args->a_ET, task_t,p_task);
  xbt_fifo_push(ps_args->a_general_workqueue, p_task);
  

  return p_task;
}

void delete_task(task_t p_task){
  if(p_task == NULL){
    return;
  }
  xbt_dynar_free(&p_task->a_ancestors);
  /*free the task itself and cleanup the pointer*/
  xbt_free_f(p_task);
}
