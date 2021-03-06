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

/**@file   pricer_coloring.h
 * @ingroup PRICERS
 * @brief  coloring variable pricer
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRICER_COLORING_H__
#define __SCIP_PRICER_COLORING_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef DEBUG
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...) 
#endif

#define IS_EDGE(i,j) ( A[i][j] != -1  && A[i][j] < m)
#define ROW_EDGE(i,j) ( n + EDGE(i,j) ) 
#define W(i,j) ( edges[A[i][j]].iTreeArc )
#define X(i,j) ( A[i][j] )
#define X0(i,j) (edges[A[i][j]].iArc0)
#define VAR_Y(i) (i)
#define VAR_Z(i,k) (n + n*(i)+k) /* i em U1+w e k em V*/
#define VAR_U(i)   (n + n*(nU1+nW) + 1 + i ) /* +1) ??*/
#define VAR_X(e)   (2*n +n*(nU1+nW) + 1 + e)
#define VAR_W(e)   (2*n +n*(nU1+nW) + 1 + nArcs + e) /* 2015-09-30 EDNA 2*m + e) */
#define VAR_F(i)   (2*n +n*(nU1+nW) + 1 + nArcs + nTreeArcs + i) /* 2015-09-30 EDNA 2*m + nTreeArcs + i)  */
#define ROW_OUT_TREE(i) (i)
#define ROW_IN_TREE(i) (n+i)
#define ROW_OUT_RING(i) (n+nU1+nW+i)
#define ROW_IN_RING(i) (2*n+nU1+nW+i)
#define ROW_DISJ(i) (3*n+nU1+nW+i)
#define ROW_DISJ_DEPOT (3*n+2*nU1+2*nW)
#define ROW_CAPACITY (3*n+2*nU1+2*nW+1)
#define ROW_U_TREE(e) (3*n+2*nU1+2*nW+2+e)
#define ROW_U_RING(e) (3*n+2*nU1+2*nW+2+nTreeArcs0 + e)
#define ROW_RING_LENGHT (5*n-2*nU2+nTreeArcs0 + nArcs0)
#define ROW_TREE_CAPACITY(k) (5*n-2*nU2+nTreeArcs0 + nArcs0 + k +1)
#define ROW_DEPOT_NOT_FACILITY (5*n-2*nU2+nTreeArcs0 + nArcs0 + n +1)
#define ROW_FACILITY(k) (5*n-2*nU2+nTreeArcs0 + nArcs0 + n +2 + k)
#define ROW_SAMEROOT(e) (5*n-2*nU2+nTreeArcs0 + nArcs0 + n +2 + n-1 + e)
#define ROW_SAMEROOT2(e) (5*n-2*nU2+nTreeArcs0 + nArcs0 + n +2 + n-1 + nTreeArcs + e)
#define ROW_FORCEF(i) (5*n-2*nU2+nTreeArcs0 + nArcs0 + n +2 + n-1 + 2*nTreeArcs + i)

/*
 * Data structures
 */

/** @brief Variable pricer data used in the \ref pricer_coloring.c "pricer" */
struct SCIP_PricerData
{
   SCIP* subscip;         /**< SCIP problem related to the pricing problem */
   SCIP_CONSHDLR*        conshdlr;           /**< comstraint handler for "same" and "diff" constraints */
  int iter; /**< total of column generation iterations that was run */
  int colsgeradas; /**< total of columns generated to current node */
  SCIP_NODE* currentNode; /**< current node */
};

extern int heur;
extern int iter;
extern double tempo;
extern double LB;
extern double UB;
extern double tempoHeur;

/** creates the coloring variable pricer and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePricerColoring(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
SCIP_RETCODE SCIPpricerColoringActivate(
   SCIP*                 scip               /**< SCIP data structure */
				    );
extern
SCIP_RETCODE pricing(
    SCIP*                 scip,               /**< SCIP data structure */
    SCIP_PRICER * pricer,
    SCIP_Bool                  isfarkas,            /**< whether we perform Farkas pricing */
    SCIP_RESULT *result
		     );

int temRestricaoBranchAtiva(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler for branching data */
			    );
  int pricingExact(
    SCIP*                 scip,               /**< SCIP data structure */
    SCIP_PRICERDATA* pricerdata,
    SCIP_Bool                  isfarkas,            /**< whether we perform Farkas pricing */
    columnType** poolOfColumns,
    int* nPoolOfColumns,
    int n,
    int m,
    int** A,
    edgeT* edges
		   );
  int heuristicOfPricing(
    SCIP*                 scip,               /**< SCIP data structure */
    SCIP_PRICERDATA* pricerdata,
    SCIP_Bool             isfarkas,            /**< whether we perform Farkas pricing */
    columnType** poolOfColumns,
    int* nPoolOfColumns,
    int n,
    int m,
    int** A,
    edgeT* edges
			 );



#ifdef __cplusplus
}
#endif

#endif
