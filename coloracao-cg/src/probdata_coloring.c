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

/**@file   probdata_coloring.c
 * @brief  Problem data for coloring problem
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file handles the main problem data used in that project. For more details see \ref PROBLEMDATA page.
 *
 * @page PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the coloring problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
 * \code
 *  ** @brief Problem data which is accessible in all places
 *  *
 *  *   This problem data is used to store the input of the coloring instance, all variables which are created, and all
 *  *   constraints.
 *  *
 * \endcode
 *
 * The function SCIPprobdataCreate(), which is called in the \ref reader_coloring.c "reader plugin" after the input file was
 * parsed, initializes the problem data structure and creates the problem in the SCIP environment. For this, it creates
 * for each item of the binpacking problem one set covering constraint and creates an initial set of variables for the
 * packings. Note that the set covering constraints have to have the <code>modifiable</code>-flag set to TRUE. This is
 * necessary to tell the solver that these constraints are not completed yet. This means, during the search new
 * variables/packings might be added.  The solver needs this information because certain reductions are not allowed.
 * See the body of the function SCIPprobdataCreate() for more details.
 *
 * A list of all interface methods can be found in probdata_coloring.h.
 **/

/*

Modificações:

*/

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* 
#define SCIP_DEBUG 
*/
#include <string.h>

#include "probdata_coloring.h"
#include "vardata_coloring.h"
#include "pricer_coloring.h"

#include "scip/cons_setppc.h"
#include "scip/scip.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the coloring, all variables which are created, and all
 * constrsaints.
 */
struct SCIP_ProbData
{
   SCIP_VAR**            vars;         /**< all exiting variables in the problem */
   SCIP_CONS**           conss;        /**< all constraints */
   int                   nvars;        /**< number of generated variables */
   int                   varssize;     /**< size of the variable array */
   int                   ncons;        /**< number of constraints */
   int                   n;            /**< total of vertices */
   int                   m;            /**< total of edges */
   int**                 A;            /**< matrix of adjacency */
   edgeT*                E;             /**< array of edges */
   const char*           probname;     /**< filename of the instance */
};

/**@name Event handler properties
 *
 * @{
 */

#define EVENTHDLR_NAME         "addedvar"
#define EVENTHDLR_DESC         "event handler for catching added variables"

/**@} */

/**@name Callback methods of event handler
 *
 * @{
 */

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecAddedVar)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARADDED);

   /*SCIPdebugMessage("exec method of event handler for added variable to probdata\n");*/

   /* add new variable to probdata */
   SCIP_CALL( SCIPprobdataAddVar(scip, SCIPgetProbData(scip), SCIPeventGetVar(event)) );

   return SCIP_OKAY;
}

/**@} */


/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   SCIP_VAR**            vars,               /**< all exist variables */
   SCIP_CONS**           conss,              /**< all pmr constraints */
   int                   nvars,              /**< number of variables */
   int                   ncons,              /**< number of constraints */
   int                   n,            /**< total of vertices */
   int                   m,            /**< total of edges */
   int**                 A,            /**< matrix of adjacency */
   edgeT*                E,             /**< array of edges */
   const char*           probname
   )
{
   assert(scip != NULL);
   assert(probdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, probdata) );

   if( nvars > 0 )
   {
      /* copy variable array */
     (*probdata)->vars=vars;
     SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*probdata)->vars, vars, nvars) ); /* NEEDED for transformed problem*/
   }
   else
      (*probdata)->vars = NULL;

   /* duplicate arrays */
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*probdata)->conss, conss, ncons) ); /* NEEDED for transformed problem */

   (*probdata)->nvars = nvars;
   (*probdata)->varssize = nvars;
   (*probdata)->ncons = ncons;
   (*probdata)->n = n;
   (*probdata)->m = m;
   (*probdata)->A = A;
   (*probdata)->E = E;
   (*probdata)->probname = probname;

   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
   int i;

   assert(scip != NULL);
   assert(probdata != NULL);

   /* release all variables */
   for( i = 0; i < (*probdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
   }

   /* release all constraints */
   for( i = 0; i < (*probdata)->ncons; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->conss[i]) );
   }

   /* free memory of arrays */
   SCIPfreeMemoryArray(scip, &(*probdata)->vars);
   SCIPfreeMemoryArray(scip, &(*probdata)->conss);
   SCIPfreeMemoryArray(scip, &(*probdata)->A);
   SCIPfreeMemoryArray(scip, &(*probdata)->E);

   /* free probdata */
   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


/** create initial columns (populate initial basis) */
static
SCIP_RETCODE createInitialColumns(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,            /**< problem data */
   SCIP_VAR** vars
   )
{
   SCIP_CONS** conss;
   SCIP_VARDATA* vardata;
   SCIP_VAR* var;
   char name[SCIP_MAXSTRLEN];

   int j, k;
   int initialColumns;
   columnType* columns, *pColumn;

   int n;
   int m;
   int** A;
   edgeT* E;

   assert(scip != NULL);
   assert(probdata != NULL);

   /* get probdata */   
   conss = SCIPprobdataGetConss(probdata);
   A = SCIPprobdataGetA(probdata);
   n = SCIPprobdataGetN(probdata);
   m = SCIPprobdataGetM(probdata);
   E = SCIPprobdataGetE(probdata);
   
   //   conss = probdata->conss;

   columns=NULL;
   initialColumns=0;
#ifdef PRIMALHEUR   
   //   initialColumns=initialHeuristic(scip, probdata, probdata->probname, n, m, A, E, &columns);
#endif

   if( initialColumns == 0 )
   {
     return SCIP_OKAY;
   }

   /* create start solution with each column */
   pColumn=columns;
   while(pColumn!=NULL){
         
     SCIP_CALL( SCIPvardataCreateColoring(scip, &vardata, pColumn->consids, pColumn->nconsids) ); 
     
     SCIP_CALL( SCIPcreateVarColoring(scip, &var, name, 1.0, FALSE, TRUE, vardata) );
     //SCIP_CALL( SCIPcreateVarColoring(scip, &var, name, 1.0, TRUE, TRUE, NULL) );

     /* add variable to correponding constraints of the PMR */
     for( j=0; j < pColumn->nconsids; j++)
       {
	 k = pColumn->consids[j];
	 SCIPdebugMessage("Nova coluna tem vertice: %d\n", k);
	 SCIP_CALL( SCIPaddCoefLinear(scip, conss[k], var, 1.0) );
       }   
	  
     /* add variable to the problem PMR */
     SCIP_CALL( SCIPaddVar(scip, var) );
     
     SCIP_CALL( SCIPprobdataAddVar(scip, probdata, var) );
     
     SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );
     pColumn = pColumn->next;
   }

   return SCIP_OKAY;
}

/**@} */


/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigColoring)
{
  /*SCIPdebugMessage("free original problem data\n");*/

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) @@@ ?*/
static
SCIP_DECL_PROBTRANS(probtransColoring)
{
   /* create transform probdata */
  SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->vars, sourcedata->conss, sourcedata->nvars, sourcedata->ncons, sourcedata->n, sourcedata->m, sourcedata->A, sourcedata->E, sourcedata->probname) );

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->ncons, (*targetdata)->conss, (*targetdata)->conss) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nvars, (*targetdata)->vars, (*targetdata)->vars) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransColoring)
{
  /*SCIPdebugMessage("free transformed problem data\n");*/

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** solving process initialization method of transformed data (called before the branch and bound process begins) */
static
SCIP_DECL_PROBINITSOL(probinitsolColoring)
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(probdata != NULL);

   /* catch variable added event */
   eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
static
SCIP_DECL_PROBEXITSOL(probexitsolColoring)
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(probdata != NULL);

   /* drop variable added event */
   eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, -1) );


   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(   
   SCIP*                 scip,         /**< SCIP data structure */
   const char*           probname,     /**< problem name */
   int                   n,            /**< total of vertices */
   int                   m,            /**< total of edges */
   int**                 A,            /**< matrix of adjacency */
   edgeT*                E             /**< array of edges */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS** conss;
   SCIP_VAR** vars, *var;
   
   char name[SCIP_MAXSTRLEN];
   int i;
   int ncons;
   int nvars;

   int j, nconsids, consids[2], *pconsids;
   SCIP_VARDATA* vardata;

   assert(scip != NULL);

   /* create event handler if it does not exist yet */
   if( SCIPfindEventhdlr(scip, EVENTHDLR_NAME) == NULL )
   {
      SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecAddedVar, NULL) );
   }

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigColoring) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransColoring) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransColoring) );
   SCIP_CALL( SCIPsetProbInitsol(scip, probinitsolColoring) );
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolColoring) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /* tell SCIP that the objective will be always integral */
   SCIP_CALL( SCIPsetObjIntegral(scip) );

   SCIP_CALL( SCIPallocBufferArray(scip, &conss, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, n) );

   for (i=0;i<n;i++){
     /* create cover constraint  */
     (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cover%d", i);
     SCIPdebugMessage("create constraint %s\n", name);
     SCIP_CALL( SCIPcreateConsBasicLinear (scip, &conss[i], name, 0, NULL, NULL, 1.0, SCIPinfinity(scip)) );
      /* declare constraint modifiable for adding variables during pricing */
      SCIP_CALL( SCIPsetConsModifiable(scip, conss[i], TRUE) );
     SCIP_CALL( SCIPaddCons(scip, conss[i]) );
   }
   
   nvars=0;
   ncons=n; 

   for(i=0; i<n; i++){
     (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lambda_%d", i);
     SCIPdebugMessage("create variable %s\n", name);

     nconsids = 0;
     consids[nconsids++]=i;
     
     SCIP_CALL( SCIPvardataCreateColoring(scip, &vardata, consids, nconsids) ); 
     
     SCIP_CALL( SCIPcreateVarColoring(scip, &var, name, 1.0, FALSE, TRUE, vardata) );

     assert(var != NULL);
     
     vars[nvars++]=var;
     
     /* add variable to corresponding constraint */
     SCIP_CALL( SCIPaddCoefLinear(scip, conss[i], var, 1.0) );     

     /* add variable to the problem */
     SCIP_CALL( SCIPaddVar(scip, var) );
         
     SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );

   }

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, vars, conss, nvars, ncons, n, m, A, E, probname) );

   SCIP_CALL( createInitialColumns(scip, probdata, vars) );

   SCIP_CALL( SCIPwriteOrigProblem(scip, "coloring.lp", "lp", FALSE) ); /* grava na saida padrao ou em file */

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   /* set user pricing data and activate the coloring pricer */
   SCIP_CALL( SCIPpricerColoringActivate(scip) ); 

   /* free local buffer arrays */
   SCIPfreeBufferArray(scip, &conss);
   SCIPfreeBufferArray(scip, &vars);  

   return SCIP_OKAY;
}

/** returns total of item */
int SCIPprobdataGetN(
  SCIP_PROBDATA*        probdata            /**< problem data */
  )
{
  return probdata->n;
}

/** returns matrix of adjacency */
int** SCIPprobdataGetA(
  SCIP_PROBDATA*        probdata            /**< problem data */
  )
{
  return probdata->A;
}

/** returns total of edges */
int SCIPprobdataGetM(
  SCIP_PROBDATA*        probdata            /**< problem data */
  )
{
  return probdata->m;
}

/** returns array of edges */
edgeT* SCIPprobdataGetE(
  SCIP_PROBDATA*        probdata            /**< problem data */
  )
{
  return probdata->E;
}

/** returns array of all variables itemed in the way they got generated */
SCIP_VAR** SCIPprobdataGetVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->vars;
}

/** returns number of variables */
int SCIPprobdataGetNVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->nvars;
}

/** returns array of set partitioning constrains */
SCIP_CONS** SCIPprobdataGetConss(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->conss;
}

/** returns array of set partitioning constrains */
int SCIPprobdataGetNcons(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->ncons;
}

/** returns Probname of the instance */
char* SCIPprobdataGetProbname(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->probname;
}

/** adds given variable to the problem data (used by pricer) */
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_VAR*             var                 /**< variables to add */
   )
{
   /* check if enough memory is left */
   if( probdata->vars == NULL )
   {
      probdata->varssize = 100;
      SCIP_CALL( SCIPallocBufferArray(scip, &probdata->vars, probdata->varssize) );
   }
   else  if( probdata->varssize == probdata->nvars )
   {
      probdata->varssize = MAX(100, probdata->varssize * 2);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &probdata->vars, probdata->varssize) );
   }

   /* capture variables */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   probdata->vars[probdata->nvars] = var;
   probdata->nvars++;

   /*SCIPdebugMessage("added variable to probdata; nvars = %d\n", probdata->nvars);*/

   return SCIP_OKAY;
}

/**@} */
