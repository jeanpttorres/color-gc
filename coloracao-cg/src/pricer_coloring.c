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

/**@file   pricer_coloring.c
 * @brief  Coloring variable pricer
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file implements the variable pricer which check if variables exist with negative reduced cost. See
 * @ref PRICER for more details.
 *
 * @page PRICER Pricing new variables
 *
 * The task of the pricer is to search for new variables with negative reduced costs. For this, the following integer
 * program is solved:
 *
 *  \f[
 *  \begin{array}[t]{rll}
 *       \max & \displaystyle \sum_{i=1}^n (\lambda_S)_i y^\star_i\\
 *        & \\
 *        subject \ to & \displaystyle \sum_{i=0}^n (\lambda_S)_i s_i \leq \kappa \\
 *        & \\
 *        & (\lambda_S)_i \in \{0,1\} & \quad \forall i \in \{ 1, \dots , n \} \\
 *  \end{array}
 * \f]
 *
 * where \f$ (\lambda_S)_i \f$ for \f$i\in\{1,\dots,n\}\f$ are binary variables and \f$y^\star_i\f$ given by the dual
 * solution of the restricted master problem. See the \ref PROBLEM "problem description" for more details.
 *
 * To solve the above integer program, we create a new SCIP instance within SCIP and use the usual functions to create
 * variables and constraints. Besides, we need the current dual solutions to all set covering constraints (each stands
 * for one item) which are the objective coefficients of the binary variables. Therefore, we use the function
 * SCIPgetDualsolSetppc() which returns the dual solutions for the given set covering constraint.
 *
 * Since we also want to generate new variables during search, we have to care that we do not generate variables over
 * and over again. For example, if we branched or fixed a certain packing to zero, we have to make sure that we do not
 * generate the corresponding variables at that node again. For this, we have to add constraints forbidding to generate
 * variables which are locally fixed to zero. See the function addFixedVarsConss() for more details. While using the
 *
 **/ 

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define SCIP_DEBUG*/

#include <assert.h>
#include <string.h>
#include <time.h>

#include "scip/cons_knapsack.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/scipdefplugins.h"

#include "probdata_coloring.h"
#include "pricer_coloring.h"
#include "vardata_coloring.h"

/**@name Pricer properties
 *
 * @{
 */

#define PRICER_NAME            "coloring"
#define PRICER_DESC            "pricer for coloring"
#define PRICER_PRIORITY        0
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */

/**@} */



/**@name Local methods
 *
 * @{
 */

/** initializes the objective function of the pricing problem */
static
SCIP_RETCODE initPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,         /**< pricer data */
   SCIP*                 subscip,            /**< pricing SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array for the pricing */
   SCIP_Bool isfarkas
   )
{
   SCIP_PROBDATA*        probdata;
   SCIP_CONS** conss;
   SCIP_CONS* cons;
   SCIP_VAR* var;
   SCIP_Real dual;

   int                   n;            /**< total of vertices */

   int e;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   assert(pricerdata != NULL);

   subscip = pricerdata->subscip;
   assert( SCIPgetStage(subscip) == SCIP_STAGE_PROBLEM );

   conss = SCIPprobdataGetConss(probdata);
   n = SCIPprobdataGetN(probdata);

   for( e = 0; e < n; e++ )
   {

      cons = conss[e];
      if(!isfarkas )
      {
         dual = SCIPgetDualsolLinear(scip, cons);
      }
      else{
         dual = SCIPgetDualfarkasLinear(scip, cons);
     }
#ifdef DEBUG2
      if(dual>EPSILON || dual < -EPSILON)
      {
         printf("\ndual[%d]=%lf", e, dual);
      }
#endif
   

      /* set objective coefficient for each variable of pricing problem */
      var = vars[e];
      /* set correctly the objective value coeficient for the variable */
      SCIP_CALL( SCIPchgVarObj(subscip, var, dual) );
     
   }

   return SCIP_OKAY;
}

/**@} */

/**name Callback methods
 *
 * @{
 */

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRICERFREE(pricerFreeColoring)
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   if( pricerdata != NULL)
   {
     /* free pricer MIP @@@ here or in PRICEREXITSOL */
     SCIP_CALL( SCIPfree(&pricerdata->subscip) );

     /* free memory */
     SCIPfreeMemory(scip, &pricerdata);
   }

   return SCIP_OKAY;
}


/** initialization method of variable pricer (called after problem was transformed) */
static
SCIP_DECL_PRICERINIT(pricerInitColoring)
{  /*lint --e{715}*/
  /*   SCIP_PRICERDATA* pricerdata;*/

   assert(scip != NULL);
   assert(pricer != NULL);

  
   return SCIP_OKAY;
}


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolColoring)
{
   SCIP_PRICERDATA* pricerdata;
   /*   SCIP* subscip;*/

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return SCIP_OKAY;
}


/** reduced cost pricing method of variable pricer for feasible LPs @@@ o metodo chamado durante a geracao de colunas */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostColoring)
{  /*lint --e{715}*/
  SCIPdebugMessage("call scip_redcost ...\n");

  /* set result pointer, see above */
  (*result) = SCIP_SUCCESS;

  /* call pricing routine */
  SCIP_CALL( pricing(scip, pricer, FALSE, result) );

  return SCIP_OKAY;
}


int pricingExact(
    SCIP*                 scip,               /**< SCIP data structure */
    SCIP_PRICERDATA* pricerdata,
    SCIP_Bool                  isfarkas,            /**< whether we perform Farkas pricing */
    columnType** poolOfColumns,
    int* nPoolOfColumns,
    int n,
    int m,
    int** A,
    edgeT* E
)
{
   SCIP* subscip;
   SCIP_CONS** pricConss;
   SCIP_SOL** sols;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   int nvars;
   int ncons;
   int nsols;
   int s, o;
   SCIP_VARDATA* vardata;

   SCIP_Real timelimit;
   SCIP_Real memorylimit;

   int i, u, v;
   int j;
   int k;
   char name[SCIP_MAXSTRLEN];

   int* consids;
   int nconss;

   SCIP_CONS** conss; /* PM's constraints */
   SCIP_PROBDATA* PMprobdata;
   int status;

   assert(scip != NULL);

   PMprobdata = SCIPgetProbData(scip);
   assert(PMprobdata != NULL);
   conss = SCIPprobdataGetConss(PMprobdata);

   /* create subscip for pricing problem */
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* create problem in sub SCIP */
   SCIP_CALL( SCIPcreateProbBasic(subscip, "pricing") );
   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(subscip, SCIP_OBJSENSE_MAXIMIZE) );

   nvars = n;
   ncons = m;
   
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   SCIP_CALL( SCIPallocBufferArray(subscip, &pricConss, ncons) );
   SCIP_CALL( SCIPallocBufferArray(subscip, &vars, nvars) );

      /* tell SCIP that the objective will be always integral */
   SCIP_CALL( SCIPsetObjIntegral(subscip) );


   for (i=0;i<m;i++){
     /* create edge constraint  */
     (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "edge%d", i);
     SCIPdebugMessage("create constraint %s\n", name);
     SCIP_CALL( SCIPcreateConsBasicLinear (subscip, &pricConss[i], name, 0, NULL, NULL, -SCIPinfinity(scip), 1.0) );
     SCIP_CALL( SCIPaddCons(subscip, pricConss[i]) );
   }
   
   for(i=0; i<n; i++){
     (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
     SCIPdebugMessage("create variable %s\n", name);

     /* create a basic variable object */
     SCIP_CALL( SCIPcreateVarBasic(subscip, &var, name, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
     assert(var != NULL);
     
     vars[i]=var;
     
     /* add variable to the problem */
     SCIP_CALL( SCIPaddVar(subscip, var) );    
   }

   for(i=0;i<m;i++){
     u=E[i].i;
     v=E[i].j;
     /* add variables x_u and x_v to corresponding edge (u,v) constraint */
     SCIP_CALL( SCIPaddCoefLinear(subscip, pricConss[i], vars[u], 1.0) );     
     SCIP_CALL( SCIPaddCoefLinear(subscip, pricConss[i], vars[v], 1.0) );     
   }
 
   pricerdata->subscip = subscip;

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* get the remaining time and memory limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
     timelimit -= SCIPgetTotalTime(scip);//SCIPgetSolvingTime(scip);

   if(timelimit < EPSILON){
     printf("\n################################\n###Pricing finished but not executed (time limit exceeded\n### total time = %lf\n", SCIPgetTotalTime(scip));
     return 0;
   }
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;

   /* set time and memory limit */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   /* initialization local pricing problem --> @@@ change objective coef with dual */

   SCIP_CALL( initPricing(scip, pricerdata, subscip, vars, isfarkas) );

   SCIP_CALL( SCIPwriteOrigProblem(subscip, "pricing_X.lp", "lp", FALSE) ); 
   /*   printf("\n################################\n###Solving pricing\n### iter=%d node atual=%lld total of nodos=%lld (left=%d) \n", pricerdata->iter, SCIPnodeGetNumber(pricerdata->currentNode), SCIPgetNNodes(scip), SCIPgetNNodesLeft(scip));*/
   
#ifdef DEBUG
   
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricing_%d.lp", pricerdata->iter);

   /* grava lp do pricing */
   SCIP_CALL( SCIPwriteOrigProblem(subscip, name, "lp", FALSE) );
   printf("\n### gravado %s\n", name);
   //   printf("\n################################\n###Solving pricing\n### iter=%d node atual=%lld total of nodos=%lld (left=%d) \n", pricerdata->iter, SCIPnodeGetNumber(pricerdata->currentNode), SCIPgetNNodes(scip), SCIPgetNNodesLeft(scip));
#endif 

   SCIPdebugMessage("solve pricer problem\n");

   /* solve sub SCIP */
   SCIP_CALL( SCIPsolve(subscip) );

   sols = SCIPgetSols(subscip);
   nsols = SCIPgetNSols(subscip);

   /* loop over all solutions and create the corresponding column to master if the reduced cost are negative for master,
    * that is the objective value i greater than 1.0 
    */

#ifdef DEBUG
   printf("### Exact Pricing finished\n### nsols=%d status=%d", nsols, SCIPgetStatus(subscip));
   //   printf("\n### current LP value=%lf iter:%d total de vars do pmr=%d\n################################\n", SCIPgetLPObjval(scip), pricerdata->iter, SCIPgetNVars(scip));
#endif
   
   /* allocate memory for new column */
   SCIP_CALL( SCIPallocBufferArray(subscip, &consids, n) );


   for( s = 0; s < nsols; ++s )
   {
      SCIP_Bool feasible;
      SCIP_SOL* sol;

      /* the soultion should be sorted w.r.t. the objective function value */
      assert(s == 0 || SCIPisFeasGE(subscip, SCIPgetSolOrigObj(subscip, sols[s-1]), SCIPgetSolOrigObj(subscip, sols[s])));

      sol = sols[s];
      assert(sol != NULL);

      /* check if solution is feasible in original sub SCIP */
      /*      SCIP_CALL( SCIPcheckSolOrig(subscip, sol, &feasible, FALSE, FALSE ) );*/
      feasible = TRUE;

      if( !feasible )
      {
         SCIPwarningMessage(scip, "solution in pricing problem is infeasible\n");
         continue;
      }
#ifdef DEBUG
        printf ("\n\n#############\n Solucao %d de cr=%lf e custo na f.o do PMR=%lf", s, 1.0 - SCIPgetSolOrigObj(subscip, sol), 1.0);
#endif

        /* check if 1 - the solution value has a value less than 0.0 */
      if( SCIPisFeasLT(subscip, 1.0 - SCIPgetSolOrigObj(subscip, sol), 0.0) ) 
      {

#ifdef DEBUG
        printf("\n New variable will be create in RMP. ");
        SCIPdebug( SCIP_CALL( SCIPprintSol(subscip, sol, NULL, FALSE) ) );
#endif

         nconss = 0;

         /* check if vertex o belongs to the independent set */
	 PRINTF("\n Vertices: { ");
         for( o = 0; o < n; ++o )
         {
            if( SCIPgetSolVal(subscip, sol, vars[o]) > 0.5 ) 
            {
               consids[nconss] = o;
               nconss++;
	       PRINTF("%d ", o);	       
            }
           else
           {
              assert( SCIPisFeasEQ(subscip, SCIPgetSolVal(subscip, sol, vars[0]), 0.0) );
           }
         }
	 PRINTF("}");

	 SCIP_CALL( SCIPvardataCreateColoring(scip, &vardata, consids, nconss) );
	     
	 (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lambda_it%d_s%d", pricerdata->iter, s);
	     
	 /* create variable for a new column with objective function coefficient 0.0 */
	 SCIP_CALL( SCIPcreateVarColoring(scip, &var, name, 1.0, FALSE, TRUE, vardata) );

	 /* add the variable data to the variable */
	 //	 SCIPvarSetData(var, vardata);  /* já setado com createVarColoring() */

	 /* add the new variable to the pricer store */
	 SCIP_CALL( SCIPaddPricedVar(scip, var, 1.0) );
	 pricerdata->colsgeradas++;
	     
	 SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );
	     
	 /* add variable to correponding constraints of the PMR */
	 for( j=0; j < nconss; j++)
	   {
	     k = consids[j];
	     SCIP_CALL( SCIPaddCoefLinear(scip, conss[k], var, 1.0) );
	   }
#ifdef DEBUG	 
	 SCIPdebug(SCIPprintVar(scip, var, NULL) );
#endif
         SCIP_CALL( SCIPreleaseVar(scip, &var) );
         
      }
   }
   SCIPfreeBufferArray(subscip, &consids);
   status = SCIPgetStatus(subscip);
   SCIP_CALL( SCIPfree(&subscip) );
   
   return ( status == SCIP_STATUS_OPTIMAL || status == SCIP_STATUS_INFEASIBLE );
}

int heuristicOfPricing(
    SCIP*                 scip,               /**< SCIP data structure */
    SCIP_PRICERDATA* pricerdata,
    SCIP_Bool             isfarkas,            /**< whether we perform Farkas pricing */
    columnType** poolOfColumns,
    int* nPoolOfColumns,
    int n,
    int m,
    int** A,
    edgeT* E
)
{
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   int* consids;
   int nconsids;

   /* code here ! */

   return 0;
} 

 
SCIP_RETCODE pricing(
    SCIP*                 scip,               /**< SCIP data structure */
    SCIP_PRICER * pricer,
    SCIP_Bool                  isfarkas,            /**< whether we perform Farkas pricing Quando tem que gerar coluna e não tem dual*/
    SCIP_RESULT *result
    )
  {
   SCIP_PROBDATA*        probdata;
   SCIP_PRICERDATA* pricerdata;
   SCIP_Bool addvar;
   columnType** poolOfColumns;
   int nPoolOfColumns;

   int n;
   int m;
   int** A;
   edgeT* E;

   
   int failedPricingExact;
   clock_t antes, agora;
   int failedPricingHeur;
   
   assert(scip != NULL);
   assert(pricer != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get the pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* get probdata */   
   A = SCIPprobdataGetA(probdata);
   n = SCIPprobdataGetN(probdata);
   m = SCIPprobdataGetM(probdata);
   E = SCIPprobdataGetE(probdata);



   if( pricerdata->currentNode == NULL || pricerdata->currentNode != SCIPgetCurrentNode(scip) )
     {
       pricerdata->iter++;
       pricerdata->colsgeradas = 0;
       pricerdata->currentNode = SCIPgetCurrentNode(scip);
     }
   else
     {
       pricerdata->iter++;
     }

   failedPricingExact=1;
   addvar = FALSE;

#ifdef DEBUG
   printf("\n################################\n### Solving pricing\n### iter=%d node atual=%lld total of nodos=%lld (left=%d) heur=%d", pricerdata->iter, SCIPnodeGetNumber(pricerdata->currentNode), SCIPgetNNodes(scip), SCIPgetNNodesLeft(scip), heur);
   printf("\n### current LP value=%lf iter=%d total de vars do pmr=%d", SCIPgetLPObjval(scip), pricerdata->iter, SCIPgetNVars(scip));
#endif
#ifdef HEUR_PRIC
   antes=clock();
   addvar = heuristicOfPricing(scip, pricerdata, isfarkas, poolOfColumns, &nPoolOfColumns, n, m, A, E);
   agora=clock();
   tempoHeur+=((double) agora-antes);///CLOCKS_PER_SEC;
#endif
   failedPricingHeur=!addvar;
   if(!addvar)
   {
     failedPricingExact = !pricingExact(scip, pricerdata, isfarkas, poolOfColumns, &nPoolOfColumns, n, m, A, E);
   }
   else{
     heur++;
   }



   if( addvar || !failedPricingExact) /* SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )*/
   {
     if( !addvar )
       {
         //	 pricerdata->colsgeradas+=nPoolOfColumns;
#ifdef DEBUG
	 printf("\n\n### Pricing finished! Current LP value=%lf iter=%d total de colunas geradas=%d", SCIPgetLPObjval(scip), pricerdata->iter, pricerdata->colsgeradas);
#endif
       }
      (*result) = SCIP_SUCCESS;
   }
#ifdef DEBUG
   SCIP_CALL( SCIPwriteTransProblem(scip, "coloring_X.lp", "lp", FALSE) ); /* grava na saida padrao ou em file */ 
   printf("\n### coloring_X.lp gravado! resultados do pricing=%d addvar=%d failed=%d\n", *result, addvar, failedPricingExact);
   getchar();
#endif
   iter=pricerdata->iter;
   
   return SCIP_OKAY;
}

/** farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasColoring)
{  /*lint --e{715}*/

  SCIPdebugMessage("call scip_farkas ...\n");

  /* call pricing routine */
  SCIP_CALL( pricing(scip, pricer, TRUE, result) );

  return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */


/** creates the coloring variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerColoring(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata;
   SCIP_PRICER* pricer;

   /* create coloring variable pricer data */
   SCIP_CALL( SCIPallocMemory(scip, &pricerdata) );

   pricerdata->conshdlr = SCIPfindConshdlr(scip, "samediff");

   /* TODO: (optional) create variable pricer specific data here */

   pricer = NULL;
   
   /* include variable pricer */
#if 0
   /* use SCIPincludePricer() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludePricer(scip, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerCopyColoring, pricerFreeColoring, pricerInitColoring, pricerExitColoring, 
         pricerInitsolColoring, pricerExitsolColoring, pricerRedcostColoring, pricerFarkasColoring,
         pricerdata) );
#else
   /* use SCIPincludePricerBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludePricerBasic(scip, &pricer, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerRedcostColoring, pricerFarkasColoring, pricerdata) );
   assert(pricer != NULL);

   /* set non fundamental callbacks via setter functions */
   /*   SCIP_CALL( SCIPsetPricerCopy(scip, pricer, pricerCopyColoring) );*/
   SCIP_CALL( SCIPsetPricerFree(scip, pricer, pricerFreeColoring) );
   SCIP_CALL( SCIPsetPricerInit(scip, pricer, pricerInitColoring) );
   /*   SCIP_CALL( SCIPsetPricerExit(scip, pricer, pricerExitColoring) );*/
   /*   SCIP_CALL( SCIPsetPricerInitsol(scip, pricer, pricerInitsolColoring) );*/
   SCIP_CALL( SCIPsetPricerExitsol(scip, pricer, pricerExitsolColoring) );
#endif

   /* add coloring variable pricer parameters */
   /* TODO: (optional) add variable pricer specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** added problem specific data to pricer and activates pricer @@@ chamado depois que o problema eh lido e pmr eh criado. Aqui, oportunidade para criar o pricing problem e setar os dados do pricing */
SCIP_RETCODE SCIPpricerColoringActivate(
   SCIP*                 scip               /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   pricerdata->iter = 0;
   pricerdata->colsgeradas = 0;
   pricerdata->currentNode = NULL;

   /* store subpric */
   pricerdata->subscip=NULL;

   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, pricer) );

   return SCIP_OKAY;

}

/**@} */
