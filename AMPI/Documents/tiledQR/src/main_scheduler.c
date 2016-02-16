/*
 *  main_scheduler.c
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/3/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */

#include "main_scheduler.h"
#include "worker.h"


unsigned int * a_ZTops;


void Scheduler_activate_tasks()
{
  msched_args_t ps_args = (msched_args_t)MSG_host_get_data(MSG_host_self());
  xbt_fifo_item_t p_cur_item;
  task_t p_cur_task;
  
  xbt_fifo_t a_tmp_fifo = xbt_fifo_copy(ps_args->a_general_workqueue);
  xbt_fifo_foreach(a_tmp_fifo,p_cur_item,p_cur_task,task_t){
    switch (p_cur_task->e_type) {
      case F:
        /*lookup in the tile updates fifo whether the task is in first position or not (if it is, then it is ready)*/
        if(xbt_fifo_get_item_content(xbt_fifo_get_first_item(ps_args->a_tile_updates[p_cur_task->task.F.ui_i*ps_args->ui_q+p_cur_task->task.F.ui_j]))==p_cur_task){
          xbt_dynar_remove_at(ps_args->a_ET, xbt_dynar_search(ps_args->a_ET, &p_cur_task), NULL);
          xbt_dynar_push_as(ps_args->a_RT, task_t,p_cur_task);
          xbt_fifo_remove(ps_args->a_general_workqueue, p_cur_task);            
        }
        break;
      case H:
        /*lookup in the tile updates fifo whether the task is in first position or not (if it is, then it is ready)*/
        if(xbt_fifo_get_item_content(xbt_fifo_get_first_item(ps_args->a_tile_updates[p_cur_task->task.H.ui_i*ps_args->ui_q+p_cur_task->task.H.ui_k]))==p_cur_task){
          /*check whether F is finished*/
            unsigned int ui_is_f = 0;
            task_t p_ancestor = xbt_dynar_get_as(p_cur_task->a_ancestors,0,task_t);
            if(xbt_fifo_is_in(ps_args->a_tile_updates[p_cur_task->task.H.ui_i*ps_args->ui_q+p_cur_task->task.H.ui_j], p_ancestor) ){
               if(p_ancestor->e_type == F){
                  if(p_ancestor->task.F.ui_j == p_cur_task->task.H.ui_j && p_ancestor->task.F.ui_i == p_cur_task->task.H.ui_i){
                    ui_is_f = 1;
                  }
               }
            }

          if(!ui_is_f){
            xbt_dynar_remove_at(ps_args->a_ET, xbt_dynar_search(ps_args->a_ET, &p_cur_task), NULL);
            xbt_dynar_push_as(ps_args->a_RT, task_t,p_cur_task);
            xbt_fifo_remove(ps_args->a_general_workqueue, p_cur_task);
          }
        }
        break;
      case Z:        
        /*lookup in both tile updates fifos whether the task is in first position or not (if it is, then it is ready)*/
        if(xbt_fifo_get_item_content(xbt_fifo_get_first_item(ps_args->a_tile_updates[p_cur_task->task.Z.ui_i*ps_args->ui_q+p_cur_task->task.Z.ui_j]))==p_cur_task){
          if(xbt_fifo_get_item_content(xbt_fifo_get_first_item(ps_args->a_tile_updates[p_cur_task->task.Z.ui_ii*ps_args->ui_q+p_cur_task->task.Z.ui_j]))==p_cur_task){
            a_ZTops[p_cur_task->task.Z.ui_i*ps_args->ui_q + p_cur_task->task.Z.ui_j] = floor(MSG_get_clock()) + p_cur_task->f_unit_cost;
            a_ZTops[p_cur_task->task.Z.ui_ii*ps_args->ui_q + p_cur_task->task.Z.ui_j] = floor(MSG_get_clock()) + p_cur_task->f_unit_cost;
            
            xbt_dynar_remove_at(ps_args->a_ET, xbt_dynar_search(ps_args->a_ET, &p_cur_task), NULL);
            xbt_dynar_push_as(ps_args->a_RT, task_t,p_cur_task);
            xbt_fifo_remove(ps_args->a_general_workqueue, p_cur_task);
          }
        }        
        break;
      case V:
        /*lookup in both tile updates fifos whether the task is in first position or not (if it is, then it is ready)*/
        if(xbt_fifo_get_item_content(xbt_fifo_get_first_item(ps_args->a_tile_updates[p_cur_task->task.V.ui_i*ps_args->ui_q+p_cur_task->task.V.ui_k]))==p_cur_task){
          if(xbt_fifo_get_item_content(xbt_fifo_get_first_item(ps_args->a_tile_updates[p_cur_task->task.V.ui_ii*ps_args->ui_q+p_cur_task->task.V.ui_k]))==p_cur_task){
            /*check whether Z is finished*/
            unsigned int ui_is_z = 0;
            task_t p_ancestor = xbt_dynar_get_as(p_cur_task->a_ancestors,0,task_t);
            if(xbt_fifo_is_in(ps_args->a_tile_updates[p_cur_task->task.V.ui_i*ps_args->ui_q+p_cur_task->task.V.ui_j], p_ancestor) ){
               if(p_ancestor->e_type == Z){
                  if(p_ancestor->task.Z.ui_j == p_cur_task->task.V.ui_j && p_ancestor->task.Z.ui_i == p_cur_task->task.V.ui_ii && p_ancestor->task.Z.ui_j == p_cur_task->task.V.ui_ii){
                    ui_is_z = 1;
                  }
               }
            }

            if(xbt_fifo_is_in(ps_args->a_tile_updates[p_cur_task->task.V.ui_ii*ps_args->ui_q+p_cur_task->task.V.ui_j], p_ancestor) ){
               if(p_ancestor->e_type == Z){
                  if(p_ancestor->task.Z.ui_j == p_cur_task->task.V.ui_j && p_ancestor->task.Z.ui_i == p_cur_task->task.V.ui_ii && p_ancestor->task.Z.ui_j == p_cur_task->task.V.ui_ii){
                    ui_is_z = 1;
                  }
               }
            }

            if(!ui_is_z){
              xbt_dynar_remove_at(ps_args->a_ET, xbt_dynar_search(ps_args->a_ET, &p_cur_task), NULL);
              xbt_dynar_push_as(ps_args->a_RT, task_t,p_cur_task);
              xbt_fifo_remove(ps_args->a_general_workqueue, p_cur_task);
            }
          }
        }         
        break;
      case ZS:
        
        /*lookup in both tile updates fifos whether the task is in first position or not (if it is, then it is ready)*/
        if(xbt_fifo_get_item_content(xbt_fifo_get_first_item(ps_args->a_tile_updates[p_cur_task->task.Z.ui_i*ps_args->ui_q+p_cur_task->task.Z.ui_j]))==p_cur_task){
          if(xbt_fifo_get_item_content(xbt_fifo_get_first_item(ps_args->a_tile_updates[p_cur_task->task.Z.ui_ii*ps_args->ui_q+p_cur_task->task.Z.ui_j]))==p_cur_task){ 
            a_ZTops[p_cur_task->task.Z.ui_i*ps_args->ui_q + p_cur_task->task.Z.ui_j] = floor(MSG_get_clock()) + p_cur_task->f_unit_cost;
            a_ZTops[p_cur_task->task.Z.ui_ii*ps_args->ui_q + p_cur_task->task.Z.ui_j] = floor(MSG_get_clock()) + p_cur_task->f_unit_cost;
            
            xbt_dynar_remove_at(ps_args->a_ET, xbt_dynar_search(ps_args->a_ET, &p_cur_task), NULL);
            xbt_dynar_push_as(ps_args->a_RT, task_t,p_cur_task);
            xbt_fifo_remove(ps_args->a_general_workqueue, p_cur_task);
          }
        }        
        break;
      case VS:
        /*lookup in both tile updates fifos whether the task is in first position or not (if it is, then it is ready)*/
        if(xbt_fifo_get_item_content(xbt_fifo_get_first_item(ps_args->a_tile_updates[p_cur_task->task.V.ui_i*ps_args->ui_q+p_cur_task->task.V.ui_k]))==p_cur_task){
          if(xbt_fifo_get_item_content(xbt_fifo_get_first_item(ps_args->a_tile_updates[p_cur_task->task.V.ui_ii*ps_args->ui_q+p_cur_task->task.V.ui_k]))==p_cur_task){          
            /*check whether ZS is finished*/
            if(!xbt_fifo_is_in(ps_args->a_tile_updates[p_cur_task->task.V.ui_i*ps_args->ui_q+p_cur_task->task.V.ui_j], xbt_dynar_get_as(p_cur_task->a_ancestors,0,task_t)) && !xbt_fifo_is_in(ps_args->a_tile_updates[p_cur_task->task.V.ui_ii*ps_args->ui_q+p_cur_task->task.V.ui_j], xbt_dynar_get_as(p_cur_task->a_ancestors,0,task_t))){
              xbt_dynar_remove_at(ps_args->a_ET, xbt_dynar_search(ps_args->a_ET, &p_cur_task), NULL);
              xbt_dynar_push_as(ps_args->a_RT, task_t,p_cur_task);
              xbt_fifo_remove(ps_args->a_general_workqueue, p_cur_task);
            }
          }
        }         
        break;
    }
    
  }
  
  xbt_fifo_free(a_tmp_fifo);
  
}


void main_scheduler(){  

  msched_args_t ps_args = (msched_args_t)MSG_host_get_data(MSG_host_self());
  
  a_ZTops = xbt_new0( unsigned int , ps_args->ui_p * ps_args->ui_q );
  
  
  
  /*call the user init function*/
  xbt_assert0(ps_args->fptr_init!=NULL,"User must provide a initialization function");
  ps_args->fptr_init();
  
  while(1){
    Scheduler_activate_tasks();
    FCFS_map();
    
    if(xbt_dynar_length(ps_args->a_ET)==0 && xbt_dynar_length(ps_args->a_RT)==0 && xbt_dynar_length(ps_args->a_running_tasks)==0){   
      /*send finalize*/
      m_host_t p_cur_proc;
      unsigned int i;
      xbt_dynar_foreach(ps_args->a_processors,i,p_cur_proc){
        m_task_t p_task = MSG_task_create("FINALIZE", 0, 0, FINALIZE);
        worker_t p_worker_data = MSG_host_get_data(p_cur_proc);
        MSG_task_isend (p_task, p_worker_data->channel);
      }
      
      /*call the user function*/
      if(ps_args->fptr_handle_finish!=NULL){
        ps_args->fptr_handle_finish(FINALIZE);
      }

      if(ps_args->uc_show_elim){
        for (unsigned int i=0;i<ps_args->ui_p;i++) {
          printf(" ");
          for (unsigned int j=0; j<ps_args->ui_q; j++) {
            
            if(i==j){
              printf(" \\s &");
            }
            else {
              printf(" %3d &",a_ZTops[i*ps_args->ui_q+j]);
            }
          }
          printf("\\\\\n");
        }        
      }
      
      xbt_free_ref(&a_ZTops);
      
      break;
    }
    
    
    m_task_t p_task = MSG_TASK_UNINITIALIZED;
    MSG_task_get(&p_task, FINISHED_TASK_CHANNEL);
    task_t p_finished_task = MSG_task_get_data(p_task);
    
    
    
    /*call the user function*/
    if(ps_args->fptr_handle_finish!=NULL){
      ps_args->fptr_handle_finish(p_finished_task);
    }
    
    
    
    /*delete the old task*/
    MSG_task_destroy(p_task);
    xbt_dynar_remove_at(ps_args->a_running_tasks, xbt_dynar_search(ps_args->a_running_tasks, &p_finished_task), NULL);
    delete_task(p_finished_task);
    
  }  
  
}
