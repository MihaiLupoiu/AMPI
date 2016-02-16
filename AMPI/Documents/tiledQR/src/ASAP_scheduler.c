/*
 *  ASAP_scheduler.c
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 2/25/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */

#include "ASAP_scheduler.h"

#include <xbt/log.h>

#include "common.h"
#include "main_scheduler.h"
#include "worker.h"
#include "mapping_policy.h"

/* Create a log channel to have nice outputs. */


int ui_desc(const void * p_ui_A, const void * p_ui_B){

    if(*((unsigned int *)p_ui_A) < *((unsigned int *)p_ui_B)){
        return 1;
    }
    else if (*((unsigned int *)p_ui_A) > *((unsigned int *)p_ui_B)){
        return -1;
    }
    else {
        return 0;
    }
}


int ui_asc(const void * p_ui_A, const void * p_ui_B){

    if(*((unsigned int *)p_ui_A) > *((unsigned int *)p_ui_B)){
        return 1;
    }
    else if (*((unsigned int *)p_ui_A) < *((unsigned int *)p_ui_B)){
        return -1;
    }
    else {
        return 0;
    }
}

void ASAP_init(){
    msched_args_t ps_args = (msched_args_t)MSG_host_get_data(MSG_host_self());  
    /*create the first tasks*/
    for(unsigned int i = 0;i<ps_args->ui_p;i++){
        task_t newtask = new_task(F, ps_args,NULL,i, 0,0,0);
    }

}

void ASAP_handle_finish(task_t p_finished_task){
    msched_args_t ps_args = (msched_args_t)MSG_host_get_data(MSG_host_self());


    if(p_finished_task==FINALIZE){
        return;
    }

    /*create new tasks and put it the either the ET or the RT*/
    switch (p_finished_task->e_type) {
        case F:
            /*create H tasks*/
            for(unsigned int k=p_finished_task->task.F.ui_j+1;k<ps_args->ui_q;k++){
                task_t newtask = new_task(H,ps_args,p_finished_task, p_finished_task->task.F.ui_i,p_finished_task->task.F.ui_j, k,0);         
            }
            break;
        case H:
            break;
        case Z:
            /*create V tasks*/
            for(unsigned int k=p_finished_task->task.Z.ui_j+1;k<ps_args->ui_q;k++){
                task_t newtask = new_task(V,ps_args,p_finished_task, p_finished_task->task.Z.ui_i,p_finished_task->task.Z.ui_ii, p_finished_task->task.Z.ui_j,k);         
            }

            break;
        case V:

            /*create F tasks if k = j+1*/
            if(p_finished_task->task.V.ui_k == p_finished_task->task.V.ui_j +1){

                xbt_assert0(xbt_fifo_size(ps_args->a_tile_updates[p_finished_task->task.V.ui_i*ps_args->ui_q+p_finished_task->task.V.ui_k])==0,"there should not be any updates at this point");
                task_t newtask = new_task(F,ps_args,NULL, p_finished_task->task.V.ui_i, p_finished_task->task.V.ui_k,0,0);   
            }
            break;
    }


    /*create Z tasks*/
    if(p_finished_task->e_type == F || p_finished_task->e_type == Z){
        unsigned int k = min( ps_args->ui_q, ps_args->ui_p);  
        for(unsigned int ui_nt = 0;ui_nt< k;ui_nt++){
            unsigned int ui_half_list;


            if(xbt_dynar_length(ps_args->a_NT[ui_nt])>=2){
                /*sort the line indexes having a triangle in the ui_nt th column*/
                xbt_dynar_sort(ps_args->a_NT[ui_nt], ui_asc);
            }


            xbt_dynar_t a_new_tasks = xbt_dynar_new(sizeof(task_t), NULL);
            while(xbt_dynar_length(ps_args->a_NT[ui_nt])>=2){
                /*take two lines and eliminate them*/
                unsigned int ui_n;


                ui_n = xbt_dynar_length(ps_args->a_NT[ui_nt]);
                ui_half_list = floor((double)ui_n/2.0);
                unsigned int ui_tl,ui_bl;
                xbt_dynar_remove_at(ps_args->a_NT[ui_nt], xbt_dynar_length(ps_args->a_NT[ui_nt])-1, &ui_bl);
                xbt_dynar_remove_at(ps_args->a_NT[ui_nt], ui_half_list + ui_n - 2*ui_half_list -1, &ui_tl);
                task_t newtask = new_task(Z,ps_args,NULL, ui_bl, ui_tl, ui_nt,0);
            }
            unsigned int length = xbt_dynar_length(a_new_tasks);

            xbt_dynar_free(&a_new_tasks);

        }

    }
}



