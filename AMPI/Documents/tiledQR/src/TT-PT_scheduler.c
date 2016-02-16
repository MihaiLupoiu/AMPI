/*
 *  TT-PT_scheduler.c
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/10/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */


#include "TT-PT_scheduler.h"
#include "task.h"
#include "common.h"
#include "main_scheduler.h"

void TTPT_init(){
  msched_args_t ps_args = (msched_args_t)MSG_host_get_data(MSG_host_self()); 
  task_t facto,update,elim;
  
  unsigned int BS = ps_args->ui_BS;
  unsigned int p = ps_args->ui_p;
  unsigned int q = ps_args->ui_q;
  int K, M, RD; 
  int k,mmm, nnn;
  
  K = min( p, q); 
  
  for (k = 0; k < K; k++) {
    for (M = k;
         //M < p-1 || M == k;  // No bottom single-row subdomain
         M < p || M == k; // The line above is PLASMA code, the line here is our variants, note: we need square tiles
         M += BS) {
      
      facto = new_task(F, ps_args,NULL,M, k,0,0);
      
      for (nnn = k+1; nnn < q; nnn++) {
        update = new_task(H,ps_args,facto, M,k, nnn,0);
      }
      
      for (mmm = M+1;
           //(mmm < M+BS && mmm < p) || mmm == p-1; // Suck in bottom single-row domain
           (mmm < M+BS && mmm < p); // The line above is PLASMA code, the line here is our variants, note: we need square tiles
           mmm++) {
        
        facto = new_task(F, ps_args,NULL,mmm, k,0,0);
        elim = new_task(Z,ps_args,NULL, mmm, M, k,0);
    
        for (nnn = k+1; nnn < q; nnn++) {
          update = new_task(H,ps_args,facto, mmm,k, nnn,0);
          update = new_task(V,ps_args,elim, mmm,M, k,nnn);
        }
      }
    }
    for (RD = BS; RD < p-k; RD *= 2) {
      for (M = k;
           //M+RD < p-1; // No reduction with bottom single-row subdomain
           M+RD < p; // The line above is PLASMA code, the line here is our variants, note: we need square tiles
           M += 2*RD) {
        
        elim = new_task(Z,ps_args,NULL, M+RD, M, k,0);

        for (nnn = k+1; nnn < q; nnn++) {
          update = new_task(V,ps_args,elim, M+RD,M, k,nnn);
        }
        
      }
      
    }    
  }
}
