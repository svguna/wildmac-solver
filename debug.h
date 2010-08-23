/*
 * wildmac-solver - returns the proper configuration of the wildmac protocol,
 * given a desired detection latency and probability.
 * Copyright (C) 2010  Stefan Guna
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see 
 * http://www.gnu.org/licenses/gpl-3.0-standalone.html.
 */
#ifndef __DEBUG_H
#define __DEBUG_H

#define TRACE

#ifdef TRACE
#define at() printf(" AT: %s:%d: %s()\n", __FILE__, __LINE__, __FUNCTION__); fflush(NULL)
#define intro() printf(" IN: %s:%d: %s() <-----\n", __FILE__, __LINE__, __FUNCTION__); fflush(NULL)
#define outro() printf("OUT: %s:%d: %s() ----->\n", __FILE__, __LINE__, __FUNCTION__); fflush(NULL)
#endif


#endif
