/*
 *  GREEDY_scheduler.c
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/10/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */


#include "GREEDY_scheduler.h"
#include "task.h"
#include "common.h"
#include "main_scheduler.h"

void GREEDY_init(){
  
  msched_args_t ps_args = (msched_args_t)MSG_host_get_data(MSG_host_self()); 
  task_t facto,update,elim;
  unsigned int p,q;
  unsigned int j, ii, jj, kk, ll, k;
  unsigned int olddone, *nT, *nZ, flag = 1;
  
  p = ps_args->ui_p;
  q = ps_args->ui_q;
  nT = (unsigned int*)malloc(q*sizeof(unsigned int)); for(j=0;j<q;j++) nT[j] = 0;
  nZ = (unsigned int*)malloc(q*sizeof(unsigned int)); for(j=0;j<q;j++) nZ[j] = 0;
  

  //  for ( time=0; time<2000; time++ ) {
  while (flag){
    for ( j=q-1; j>0; j-- ){
      
      if (( j == q-1 ) && ( nZ[j] == p-q ) && ( p != q )){ 
        flag = 0;
        break;
      }
      
      if (( p == q ) && (nT[p-1] == 1)){
        flag = 0;
        break;
      }
      
    
        olddone = nZ[j];
        nZ[j] = nZ[j] + ( nT[j] - nZ[j] ) / 2 ; 
        for ( k = olddone; k < nZ[j]; k++ ){
          ll = p-k-nZ[j]+olddone-1;
          ii = p-k-1;
          
          elim = new_task(Z,ps_args,NULL, ii, ll, j,0);
          for ( jj = j+1; jj < q ; jj++ ){
            update = new_task(V,ps_args,elim, ii,ll, j,jj); 
            
          }
          
        }
        
        olddone = nT[j];
        nT[j] = nZ[j-1];
        for ( k = olddone; k < nT[j]; k++){
          kk = p-k-1;
          facto = new_task(F, ps_args,NULL,kk, j,0,0);
          
          for ( jj = j+1; jj < q ; jj++ ){
            update = new_task(H,ps_args,facto, kk,j, jj,0);         
          }
        }

      
    }
    
    j = 0;
    
    if (( j == q-1 ) && ( nZ[j] == p-q ) && ( p != q )){
      flag = 0;
      break;
    }
    olddone = nZ[j];
    nZ[j] = nZ[j] + ( nT[j] - nZ[j] ) / 2 ;
    for ( k = olddone; k < nZ[j] ; k++ ){
      ll = p-k-nZ[j]+olddone-1;
      ii = p-k-1;
      elim = new_task(Z,ps_args,NULL, ii, ll, 0,0);
      
      for ( jj = 1; jj < q ; jj++ ){
        update = new_task(V,ps_args,elim, ii,ll, 0,jj);  
      }
      
    }
    
    if ( p-nT[j] > 0 ){
      for ( k = p; k > 0 ; k-- ){
        facto = new_task(F, ps_args,NULL,k-1, 0,0,0);
        
        for ( jj = 1; jj < q ; jj++ ){
          update = new_task(H,ps_args,facto, k-1,0, jj,0);
        }
      }
      nT[j] = nT[j] + ( p - nT[j] );
    }
    
  }
  
  free(nT);
  free(nZ);
}

