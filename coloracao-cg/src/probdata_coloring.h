/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probdata_coloring.h
 * @brief  Problem data for coloring problem
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file handles the main problem data used in that project. For more details see \ref PROBLEMDATA page.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_COLORING__
#define __SCIP_PROBDATA_COLORING__

#include "scip/scip.h"

#define EPSILON 0.000001
typedef struct{
  int i;
  int j;
} edgeT;

typedef struct column
{
   int nconsids; /**< total of constraints*/
   int* consids; /**< index of PMR constraints */
  struct column* next;
} columnType;

/** sets up the problem data */
extern
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,         /**< SCIP data structure */
   const char*           probname,     /**< problem name */
   int                   n,            /**< total of vertices */
   int                   m,            /**< total of edges */
   int**                 A,            /**< matrix of adjacency */
   edgeT*                E             /**< array of edges */
);


/** returns total of itens */
extern
int SCIPprobdataGetN(
  SCIP_PROBDATA*        probdata            /**< problem data */
		     );
extern
int SCIPprobdataGetM(
  SCIP_PROBDATA*        probdata            /**< problem data */
		     );

extern
int** SCIPprobdataGetA(
  SCIP_PROBDATA*        probdata            /**< problem data */
		     );

/** returns array of edges */
extern
edgeT* SCIPprobdataGetE(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns number of variables */
extern
int SCIPprobdataGetNVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of constrains */
extern
SCIP_CONS** SCIPprobdataGetConss(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of vars */
extern
SCIP_VAR** SCIPprobdataGetVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );


#endif
