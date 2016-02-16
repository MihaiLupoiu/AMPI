/*
 *  TS-FT_scheduler.c
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/5/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */

#include "TS-FT_scheduler.h"

#include <xbt/log.h>

#include "common.h"
#include "main_scheduler.h"
#include "worker.h"
#include "mapping_policy.h"

/* Create a log channel to have nice outputs. */
XBT_LOG_EXTERNAL_DEFAULT_CATEGORY(scheduler_log);

void TSFT_init(){
  msched_args_t ps_args = (msched_args_t)MSG_host_get_data(MSG_host_self()); 
  task_t facto,update,elim;

  for (unsigned int l = 0; l < ps_args->ui_q ; l++ ){
    
    facto = new_task(F, ps_args,NULL,l, l,0,0);
    
    
    for (unsigned int j = l+1; j < ps_args->ui_q ; j++ ){
      update = new_task(H,ps_args,facto, l,l, j,0);
      
    }
    
    for (unsigned int i = l+1; i < ps_args->ui_p ; i++ ){
      
      elim = new_task(ZS,ps_args,NULL, i, l, l,0);
      
      for (unsigned int j = l+1; j < ps_args->ui_q ; j++ ){
        update = new_task(VS,ps_args,elim, i,l, l,j);
      }
      
    }
    
  }
}




