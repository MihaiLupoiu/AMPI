/*
 *  main.c
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/10/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */


#include <stdio.h>
#include <getopt.h>

#include "common.h"

#include "GREEDY_scheduler.h"
#include "TT-FT_scheduler.h"
#include "TS-FT_scheduler.h"
#include "ASAP_scheduler.h"
#include "TS-PT_scheduler.h"
#include "TT-PT_scheduler.h"
#include "MC_scheduler.h"
#include "main_scheduler.h"
#include "worker.h"


int operation_cnt;


int read_arguments(int argc, char *argv[],msched_args_t ps_args) {
  int c;

  /*set default options*/
  ps_args->ui_BS = 1;
  ps_args->fptr_init = GREEDY_init;
  ps_args->fptr_handle_finish = NULL;
  ps_args->lf_delay = 0.0;
  ps_args->ui_tile_size = 1;
  ps_args->uc_show_elim = 0;
  ps_args->uc_coarse = 0;
  
  while (1) {
    static struct option long_options[] = {
      /* These options set a flag.   */
      { "help", no_argument, 0, 'h' },
      { "show-elim", no_argument, 0, 'e' },
      { "coarse", no_argument, 0, 'c' },
      { "rowcount", required_argument, 0, 'p' },
      { "colcount", required_argument, 0, 'q' },      
      { "blocksize", required_argument, 0, 'b' },
      { "algorithm", required_argument, 0, 'a' },
      { "domainsize", required_argument, 0, 't' },
      { 0, 0, 0, 0 }
    };
    /* getopt_long stores the option index here.   */
    int option_index = 0;
    c = getopt_long(argc, argv, "hp:q:b:d:a:t:", long_options, &option_index);
    
    
    /* Detect the end of the options.   */
    if (c == -1)
      break;
    
    
    switch (c) {
      case 'h':
        printf("Usage: %s  platform_file.xml  --rowcount P --colcount Q [ --algorithm ALGO [--domainsize BS] --blocksize BLK --coarse --show-elim]\n\
                    platform_file.xml is a file generated using the generate_platform.pl script\n\
                    -a or --algorithm can be either GREEDY, ASAP, MC, FT-TT, FT-TS, PT-TT or PT-TS, default is GREEDY\n\
                    -t or --domainsize let the user specify the domain size for PT algorithms, default value is 1\n\
                    -c or --coarse toggle the coarse-grain model\n\
                    -e or --show-elim display the elimination time steps for each tile\n\
                    -h or --help print this message\n", argv[0]);
        exit(0);
        break;
      case 'b':
        ps_args->ui_tile_size = atoi(optarg);
        break;
      case 'e':
        ps_args->uc_show_elim = 1;
        break;
      case 'c':
        ps_args->uc_coarse = 1;
        break;        
      case 'p':
        ps_args->ui_p = atoi(optarg);
        break;
      case 'q':
        ps_args->ui_q = atoi(optarg);
        break;
      case 'a':
      {
        
        if(strcmp(optarg,"FT-TT")==0){
          ps_args->fptr_init = TTFT_init;
          ps_args->fptr_handle_finish = NULL;
        }
        else if(strcmp(optarg,"FT-TS")==0){
          ps_args->fptr_init = TSFT_init;
          ps_args->fptr_handle_finish = NULL;
        }
        else if(strcmp(optarg,"ASAP")==0){
          ps_args->fptr_init = ASAP_init;
          ps_args->fptr_handle_finish = ASAP_handle_finish;
        }
        else if(strcmp(optarg,"GREEDY")==0){
          ps_args->fptr_init = GREEDY_init;
          ps_args->fptr_handle_finish = NULL;
        }
        else if(strcmp(optarg,"MC")==0){
          ps_args->fptr_init = MC_init;
          ps_args->fptr_handle_finish = NULL;
        }
        else if(strcmp(optarg,"PT-TS")==0){
          ps_args->fptr_init = TSPT_init;
          ps_args->fptr_handle_finish = NULL;
        }
        else if(strcmp(optarg,"PT-TT")==0){
          ps_args->fptr_init = TTPT_init;
          ps_args->fptr_handle_finish = NULL;
        } 
      }
        break;
      case 't':
        ps_args->ui_BS = atoi(optarg);
        break;
      default:
        exit(0);
    }
    
  }
  return 0;
}










int main (int argc, const char * argv[]) {
  
  MSG_error_t res = MSG_OK;
  
  srand(time(0));
  
  
  MSG_global_init(&argc, argv);
  if (argc < 6) {
        printf("Usage: %s  platform_file.xml  --rowcount P --colcount Q [ --algorithm ALGO [--domainsize BS] --blocksize BLK --coarse --show-elim]\n\
                    platform_file.xml is a file generated using the generate_platform.pl script\n\
                    -a or --algorithm can be either GREEDY, ASAP, MC, FT-TT, FT-TS, PT-TT or PT-TS, default is GREEDY\n\
                    -t or --domainsize let the user specify the domain size for PT algorithms, default value is 1\n\
                    -c or --coarse toggle the coarse-grain model\n\
                    -e or --show-elim display the elimination time steps for each tile\n\
                    -h or --help print this message\n", argv[0]);
    exit(1);
  }
  
  operation_cnt = 0;
  
  MSG_set_channel_number(MAX_CHANNEL);
  MSG_create_environment(argv[1]);


  /*create main scheduler args*/
  msched_args_t ps_args = xbt_new0(s_msched_args_t,1);
  ps_args->a_processors = xbt_dynar_new(sizeof(m_host_t), NULL);
  ps_args->a_general_workqueue = xbt_fifo_new();
  ps_args->a_running_tasks = xbt_dynar_new(sizeof(task_t), NULL);
  
  
  m_host_t * a_hosts = MSG_get_host_table();
  for(int i = 0;i<MSG_get_host_number();i++){
    if(strcmp(MSG_host_get_name(a_hosts[i]), "master")!=0){
      xbt_dynar_push_as(ps_args->a_processors,m_host_t,a_hosts[i]);
    }
  }

  /*parse the options*/
  read_arguments(argc, argv, ps_args);
  
  ps_args->ui_p = ps_args->ui_p/ps_args->ui_tile_size;
  ps_args->ui_q = ps_args->ui_q/ps_args->ui_tile_size;
  ps_args->a_ET = xbt_dynar_new(sizeof(task_t),NULL);
  ps_args->a_RT = xbt_dynar_new(sizeof(task_t),NULL);   
  ps_args->a_NT = xbt_new0(xbt_dynar_t,ps_args->ui_p);
  int cnt = ps_args->ui_p*ps_args->ui_q;
  ps_args->a_tile_updates = xbt_new0(xbt_fifo_t,cnt);
  
  for(unsigned int i = 0;i<ps_args->ui_p;i++){
    /*create a dynar per line*/
    ps_args->a_NT[i] = xbt_dynar_new(sizeof(unsigned int), NULL);

    /*create a fifo per tile (in order to store updates H & V )*/
    for(unsigned int j = 0;j<ps_args->ui_q;j++){
      ps_args->a_tile_updates[i*ps_args->ui_q + j] = xbt_fifo_new();
    }
  }
  
  

  MSG_host_set_data(MSG_get_host_by_name("master"), ps_args);
  MSG_process_create("main_scheduler", main_scheduler, NULL, MSG_get_host_by_name("master"));
  
  unsigned int i;
  m_host_t p_cur_host;
  xbt_dynar_foreach(ps_args->a_processors,i,p_cur_host){
    worker_t p_worker_data = xbt_new0(s_worker_t,1);
    sprintf(p_worker_data->channel,"%s",MSG_host_get_name(p_cur_host));
    p_worker_data->a_workqueue = xbt_fifo_new();
    p_worker_data->lf_delay = ps_args->lf_delay;
    p_worker_data->ps_ms_args = ps_args;
    MSG_host_set_data(p_cur_host, p_worker_data);
    MSG_process_create("worker", worker, NULL, p_cur_host);
  }
  
  res = MSG_main();
  

  double f_time = MSG_get_clock();
  printf("Simulation time : %f\n",f_time);
  printf("Weighted task count : %d\n",operation_cnt);

  
  
  res = MSG_clean();


  /*free main scheduler args*/
  
  for(unsigned int i = 0;i<ps_args->ui_p;i++){
    xbt_dynar_free(&ps_args->a_NT[i]);
    for(unsigned int j = 0;j<ps_args->ui_q;j++){
      xbt_fifo_free(ps_args->a_tile_updates[i*ps_args->ui_q + j]);
    }
  }
  xbt_fifo_free(ps_args->a_general_workqueue);
  xbt_free_ref(&ps_args->a_tile_updates);
  xbt_dynar_free(&ps_args->a_processors);
  xbt_dynar_free(&ps_args->a_running_tasks);
  xbt_dynar_free(&ps_args->a_ET);
  xbt_dynar_free(&ps_args->a_RT); 
  xbt_free_ref(&ps_args->a_NT);
  
  
  xbt_free_ref(&ps_args);

  
  if (res == MSG_OK)
    return 0;
  else
    return 1;
}
