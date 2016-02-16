/*
 *  worker.h
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/4/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */

#ifndef HEADER_WORKER
#define HEADER_WORKER

#include "task.h"
#include <xbt/fifo.h>
#include <xbt/dict.h>

#include "main_scheduler.h"

typedef struct s_worker_t {
  char channel[20];
  task_t p_last_task;
  xbt_fifo_t a_workqueue;
  double lf_delay;  
  msched_args_t ps_ms_args;
} s_worker_t, * worker_t;

void worker();

#endif
