/*
 *  common.h
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/4/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */

#ifndef HEADER_COMMON
#define HEADER_COMMON

#include "msg/msg.h"            /* Yeah! If you want to use msg, you need to include msg/msg.h */
#include "xbt/dict.h"
#include "xbt/fifo.h"
#include "xbt/sysdep.h"         /* calloc, printf */
#include "xbt/log.h"
#include "xbt/asserts.h"


#define FINALIZE ((void*)221297)        /* a magic number to tell people to stop working */



typedef enum {
  GREEDY_CHANNEL = 0,
  FINISHED_TASK_CHANNEL,
  MAX_CHANNEL
} channel_t;

struct s_task_t;
typedef struct s_task_t * task_t;
struct s_msched_args_t;
typedef struct s_msched_args_t * msched_args_t;




#endif
