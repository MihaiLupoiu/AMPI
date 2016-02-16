/*
 *  mapping_policy.c
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/6/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */

#include "mapping_policy.h"

#include <msg/msg.h>

#include "worker.h"
#include "main_scheduler.h"


int load_asc(const void * p_hostA, const void * p_hostB){
  
  float lf_loadA,lf_loadB;
	worker_t p_procA;
	worker_t p_procB;
  
	p_procA = MSG_host_get_data(*(m_host_t*)p_hostA);
 	p_procB = MSG_host_get_data(*(m_host_t*)p_hostB);
  
  lf_loadA = 0;
  xbt_fifo_item_t cur_item = NULL;
  task_t cur_task;
  xbt_fifo_foreach(p_procA->a_workqueue,cur_item,cur_task,task_t){
    lf_loadA+=cur_task->f_unit_cost;
  }
  
  lf_loadB = 0;
  cur_item = NULL;
  cur_task;
  xbt_fifo_foreach(p_procB->a_workqueue,cur_item,cur_task,task_t){
    lf_loadB+=cur_task->f_unit_cost;
  }  
  
	if(lf_loadA > lf_loadB){
		return 1;
	}
	else if (lf_loadA < lf_loadB){
		return -1;
	}
	else {
		return 0;
	}
}



void FCFS_map()
{
  msched_args_t ps_args = (msched_args_t)MSG_host_get_data(MSG_host_self());
  
  /*map tasks from the ready task list onto idle processors*/
  /*map ready tasks*/
  unsigned int j = 0;
  for (unsigned int i = 0; i<xbt_dynar_length(ps_args->a_RT); i++,j++) {
    task_t p_cur_task;

    /*sort the processors in non descending order of load*/
    xbt_dynar_sort(ps_args->a_processors, load_asc);
    m_host_t p_proc =xbt_dynar_get_as(ps_args->a_processors,0,m_host_t);
    worker_t p_worker_data = MSG_host_get_data(p_proc);
    if(xbt_fifo_size(p_worker_data->a_workqueue)){
      i++;
      break;
    }
    
    
    xbt_dynar_remove_at(ps_args->a_RT, i, &p_cur_task);
    i--;    
    
    
    xbt_dynar_push_as(ps_args->a_running_tasks,task_t,p_cur_task);
    
    m_task_t task = MSG_task_create(p_cur_task->name, 0, 0, p_cur_task);
    
    
    p_cur_task->p_proc_map = p_proc;

    xbt_fifo_push(p_worker_data->a_workqueue, p_cur_task);
    
    
    MSG_task_isend (task, p_worker_data->channel);
    
  }
}

