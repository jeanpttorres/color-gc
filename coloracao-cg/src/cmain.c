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

/**@file   cmain.c
 * @brief  Main file for coloring pricing example
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 *  This the file contains the \ref main() main function of the projects. This includes all the default plugins of
 *  \SCIP and the once which belong to that projects. After that is starts the interactive shell of \SCIP or processes
 *  the shell arguments if given.
 */

/** \mainpage B&P for the vertex coloring problem Index Page
 *
 * \section Introduction
 *
 * This project refers to a branch-and-price implementation to solve the vertex coloring problem.
 *
 */

/* 2015-09-30 EDNA: incluido codigo para gerar um layout de grafo e setar execucao somente no no raiz quando ONLYROOT ativado*/

#include <stdio.h>
#include <time.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"

#include "vardata_coloring.h"
#include "probdata_coloring.h"
#include "pricer_coloring.h"
#include "reader_coloring.h"

#include "parameters_coloring.h"

int heur;
int iter;
clock_t antes, agora;
double tempo, LB, UB, tempoHeur;

void printSol(SCIP* scip, char* filename);


/** print best primal solution in file with extension ".sol"
 *  
 */
void printSol(SCIP* scip, char *filename)
{
  FILE* file, *graphFile; /* 2015-09-30 EDNA */
   SCIP_SOL* bestSolution;
   int* consids, nconsids, i, v;
   SCIP_PROBDATA*        probdata;
   SCIP_VAR** vars, *var;
   SCIP_Real solval;
   SCIP_VARDATA* vardata;
   int n;
   int m;
   edgeT* E;
   
   int nvars;
   int cont;
   int totalCores;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);

   bestSolution = SCIPgetBestSol(scip);
   if( bestSolution == NULL )
     return;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   nvars = SCIPprobdataGetNVars(probdata);
   vars = SCIPprobdataGetVars(probdata);
   n = SCIPprobdataGetN(probdata);
   m = SCIPprobdataGetM(probdata);
   E = SCIPprobdataGetE(probdata);

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s.sol", filename); 

   file = fopen(name, "w");

   if(!file)
     {
       printf("\nProblem to create solution file: %s", name);
       return;
     }
   name[0]='\0';
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s.gr", filename); 
   graphFile = fopen(name, "w");

   if(!graphFile)
     {
       printf("\nProblem to create solution file: %s", name);
       return;
     }

   totalCores = SCIPsolGetOrigObj(bestSolution);
   fprintf(graphFile, "graph {\n");
   fprintf(file, "\nValue: %d", totalCores);
  
   for( v=0, cont=0; v< nvars; v++ )
     {
       var = vars[v];
       solval = SCIPgetSolVal(scip, bestSolution, var);
       if( solval > EPSILON )
	 {
	   cont++;
	   vardata = SCIPvarGetData(var);
	   assert(vardata!=NULL);
	   nconsids = SCIPvardataGetNConsids(vardata);
	   
	   consids = SCIPvardataGetConsids(vardata);
	   assert(consids!=NULL);
	   for( i=0; i<nconsids && consids[i] < n; i++ ) 
	     {
	       fprintf(file, "\n%d (color %d) ", consids[i], cont);
	       fprintf(graphFile, "%d [style=filled, fillcolor=\"%lf %lf %lf\"];\n", consids[i], ((double) cont+1.0)/(totalCores+1), 0.5,0.7); 	       
	     }
	 }
     }
   for( v=0; v< m; v++){
     fprintf(graphFile, "%d--%d [splines=ortho];\n", E[v].i, E[v].j); 	       
   }
   fprintf(graphFile, "}");

   fclose(file);
   fclose(graphFile);
}

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runShell(
   int                        argc,               /**< number of shell parameters */
   char**                     argv,               /**< array with shell parameters */
   const char*                defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;
   char filename[300];
   FILE* fout;

   antes = clock();

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   
   /* include coloring reader */
   SCIP_CALL( SCIPincludeReaderColoring(scip) );
  
  /* include coloring pricer  */
   SCIP_CALL( SCIPincludePricerColoring(scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* for column generation instances, disable restarts */
   SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts",0) );

   /* turn off all separation algorithms */
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable heuristics */
   /*   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) ); */

   /* disable presolving */
   /* SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) ); */  

   SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxrounds", 0) );  
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) ); /* only at root node */
   SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority", 1000000) );

   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxrounds", 0) );


   /* set time limit */
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", TIME_LIMIT) );

#ifdef ONLYROOT
   /* run only at root node */
   SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1) );
#endif
   
   /**********************************
    * Process command line arguments *
    **********************************/
   tempoHeur=0;

   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   agora = clock();

   tempo = (agora-antes)/((double) CLOCKS_PER_SEC);

   filename[0]='\0';
#ifdef HEUR_PRIC
   sprintf(filename, "%s.1", argv[2]);
#else
   sprintf(filename, "%s.0", argv[2]);
#endif
   
#ifdef DEBUG
   strcat(filename, ".deb");
#endif

#ifdef ONLYROOT
   strcat(filename, ".root");
#endif
   
   
   printSol(scip, filename);

   strcat(filename, ".out");
   
#ifdef DEBUG
   printf("DEBUG ATIVADO\n");
#endif


   fout = fopen(filename, "w");
   if(!fout){
     printf("\nProblema na criacao do arquivo: %s", filename);
   }
   else{
     UB = SCIPgetPrimalbound(scip);
     LB = SCIPgetDualbound(scip);
#ifdef HEUR_PRIC
     fprintf(fout,"%s;1;%d;%d;%lli;%lf;%lf;%lf;%lf;%lf;%lf;%lli;%d;%lf;%lf;%lli;%d;%d\n", argv[2],heur, iter, SCIPgetNRootLPIterations(scip), tempo, LB, UB, SCIPgetGap(scip), SCIPgetDualboundRoot(scip),((double)tempoHeur)/CLOCKS_PER_SEC, SCIPgetNTotalNodes(scip), SCIPgetNNodesLeft(scip), SCIPgetSolvingTime(scip),SCIPgetTotalTime(scip),SCIPgetMemUsed(scip),SCIPgetNLPCols(scip),SCIPgetStatus(scip));
#else
     fprintf(fout,"%s;0;%d;%d;%lli;%lf;%lf;%lf;%lf;%lf;%lf;;%lli;%d;%lf;%lf;%lli;%d;%d\n", argv[2],heur, iter, SCIPgetNRootLPIterations(scip), tempo, LB, UB, SCIPgetGap(scip), SCIPgetDualboundRoot(scip), ((double)tempoHeur)/CLOCKS_PER_SEC,SCIPgetNTotalNodes(scip), SCIPgetNNodesLeft(scip), SCIPgetSolvingTime(scip),SCIPgetTotalTime(scip),SCIPgetMemUsed(scip),SCIPgetNLPCols(scip),SCIPgetStatus(scip));
#endif
     fclose(fout);
   }
   /********************
    * Deinitialization *
    ********************/

   /*   SCIP_CALL( SCIPfree(&scip) ); */

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   heur=0;
   retcode = runShell(argc, argv, "scip.set");
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}

