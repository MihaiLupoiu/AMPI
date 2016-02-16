/*
 *  ASAP_scheduler.h
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 2/25/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */

#ifndef HEADER_ASAP
#define HEADER_ASAP

#include "task.h"




void ASAP_handle_finish(task_t p_finished_task);
void ASAP_init();
int ui_desc(const void * p_ui_A, const void * p_ui_B);
int ui_asc(const void * p_ui_A, const void * p_ui_B);
#endif
