/*
 *  worker.c
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/4/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */

#include "worker.h"

#include "task.h"
#include "common.h"




void worker(){
  worker_t p_me = MSG_host_get_data(MSG_host_self());

  while(1){
    m_task_t p_task = MSG_TASK_UNINITIALIZED;
    
    
    msg_comm_t comm =	MSG_task_irecv (&p_task, p_me->channel);
    MSG_comm_wait(comm, -1);
    MSG_comm_destroy(comm);
    
    task_t p_todo_task = MSG_task_get_data(p_task);
    if(p_todo_task!=FINALIZE){
    
      
      MSG_task_destroy(p_task);
      double delay = 1.0 + p_me->lf_delay*(double)rand()/(double)(RAND_MAX);
      if(delay<0){
        delay=0;
      }
      
      if (p_me->ps_ms_args->uc_coarse) {
        p_task = MSG_task_create(p_todo_task->name, delay*p_todo_task->f_unit_cost*MSG_get_host_speed(MSG_host_self()), 0, NULL);
      }
      else{
        //p_task = MSG_task_create(p_todo_task->name, delay*p_todo_task->f_unit_cost*pow(p_todo_task->ui_tile_size,3)/3.0, 0, NULL);
        p_task = MSG_task_create(p_todo_task->name, delay*p_todo_task->f_unit_cost*MSG_get_host_speed(MSG_host_self()), 0, NULL);
      }
      MSG_task_execute(p_task);

      xbt_fifo_shift(p_me->a_workqueue);
      p_me->p_last_task = p_todo_task;
      
    p_todo_task->uc_done = 1; 
    
      
    if(p_todo_task->e_type == Z || p_todo_task->e_type == ZS ){
        /*put the topmost line again in the triangles*/
        xbt_dynar_push_as(p_me->ps_ms_args->a_NT[p_todo_task->task.Z.ui_j],unsigned int,p_todo_task->task.Z.ui_ii);
    }
    else if (p_todo_task->e_type == F) {
      if(p_todo_task->task.F.ui_i>=p_todo_task->task.F.ui_j){
        /*push the line in NT*/
        xbt_dynar_push_as(p_me->ps_ms_args->a_NT[p_todo_task->task.F.ui_j],unsigned int,p_todo_task->task.F.ui_i);
      }
    }

      switch (p_todo_task->e_type) {
        case F:
          /*remove F from updates*/
          xbt_fifo_shift(p_me->ps_ms_args->a_tile_updates[p_todo_task->task.F.ui_i*p_me->ps_ms_args->ui_q+p_todo_task->task.F.ui_j]);
          break;
        case H:
          /*remove H from updates*/
          xbt_fifo_shift(p_me->ps_ms_args->a_tile_updates[p_todo_task->task.H.ui_i*p_me->ps_ms_args->ui_q+p_todo_task->task.H.ui_k]);
          break;
        case Z:
          /*remove Z from updates*/        
          xbt_fifo_shift(p_me->ps_ms_args->a_tile_updates[p_todo_task->task.Z.ui_i*p_me->ps_ms_args->ui_q+p_todo_task->task.Z.ui_j]);        
          xbt_fifo_shift(p_me->ps_ms_args->a_tile_updates[p_todo_task->task.Z.ui_ii*p_me->ps_ms_args->ui_q+p_todo_task->task.Z.ui_j]);        
          break;
        case V:
          /*remove V from updates*/
          xbt_fifo_shift(p_me->ps_ms_args->a_tile_updates[p_todo_task->task.V.ui_i*p_me->ps_ms_args->ui_q+p_todo_task->task.V.ui_k]);
          xbt_fifo_shift(p_me->ps_ms_args->a_tile_updates[p_todo_task->task.V.ui_ii*p_me->ps_ms_args->ui_q+p_todo_task->task.V.ui_k]);
          break;
        case ZS:
          /*remove Z from updates*/        
          xbt_fifo_shift(p_me->ps_ms_args->a_tile_updates[p_todo_task->task.Z.ui_i*p_me->ps_ms_args->ui_q+p_todo_task->task.Z.ui_j]);        
          xbt_fifo_shift(p_me->ps_ms_args->a_tile_updates[p_todo_task->task.Z.ui_ii*p_me->ps_ms_args->ui_q+p_todo_task->task.Z.ui_j]);        
          break;
        case VS:
          /*remove V from updates*/
          xbt_fifo_shift(p_me->ps_ms_args->a_tile_updates[p_todo_task->task.V.ui_i*p_me->ps_ms_args->ui_q+p_todo_task->task.V.ui_k]);
          xbt_fifo_shift(p_me->ps_ms_args->a_tile_updates[p_todo_task->task.V.ui_ii*p_me->ps_ms_args->ui_q+p_todo_task->task.V.ui_k]);
          break;        
      }  

    MSG_task_destroy(p_task);
    p_task = MSG_task_create("handshake", 0, 0, p_todo_task);
    MSG_task_put(p_task, MSG_get_host_by_name("master"), FINISHED_TASK_CHANNEL);
    
    }
    else {
      MSG_task_destroy(p_task);
      xbt_fifo_free(p_me->a_workqueue);
      xbt_free_ref(&p_me);
      break;
    }

  }
  
}
