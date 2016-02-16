/*
 *  main_scheduler.h
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/3/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */

#ifndef HEADER_MSCHED
#define HEADER_MSCHED

#include "task.h"
#include "common.h"
#include <xbt/dynar.h>
#include <xbt/fifo.h>

typedef struct s_msched_args_t{
    /*list of processors*/
    xbt_dynar_t a_processors;

    unsigned char uc_show_elim;
    unsigned char uc_coarse;

    void (*fptr_handle_finish)(task_t p_task);
    void (*fptr_init)(void);

    xbt_fifo_t a_general_workqueue;
    xbt_dynar_t a_running_tasks;
    xbt_fifo_t * a_tile_updates;

    /*array of ui_q list of line indexes : a_NT(x)
      contains lines having x zeros under the diagonal*/ 
    xbt_dynar_t * a_NT;

    /*list of created tasks*/
    xbt_dynar_t a_ET;
    /*list of ready tasks*/
    xbt_dynar_t a_RT;

    /*matrix leading dimensions*/
    unsigned int ui_p,ui_q,ui_tile_size;
    unsigned int ui_BS;
    double lf_delay;
} s_msched_args_t;

void Scheduler_activate_tasks();
void main_scheduler();

#endif
