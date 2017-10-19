/*
 *  CellsSystem-F.cpp
 *  Sim3D
 *
 *  Created by Edoardo Milotti on 22/04/10.
 *  Copyright 2010 I.N.F.N.-Sezione di Trieste. All rights reserved.
 *
 *  This file prints flows.
 */


#include "CellsSystem.h"
using namespace vbl;
//  ******************** PrintFlows ********************

// ***************************************************************
// Calculation and printing of extracellular flows (binary only)
// ***************************************************************
//
void CellsSystem::PrintFlows()
{

	// apertura del file di output
	char str[100];	
	sprintf(str,"%lu.bin",nconfiguration);
	std::string flow_record = flow_b_filename+str;
	flow_b_file.open( flow_record.c_str(), std::ios::binary );
	

	// concentrazioni ambientali (comuni a tutte le cellule)
	double concO2env = Env.GetEnvironmentO2()/Env.GetEnvironmentvolume();
	double concGenv = Env.GetEnvironmentG()/Env.GetEnvironmentvolume();
	double concAenv = Env.GetEnvironmentA()/Env.GetEnvironmentvolume();
	double concAcLenv = Env.GetEnvironmentAcL()/Env.GetEnvironmentvolume();
		
		
	for(unsigned long n=0; n<ncells; n++)	// loop sulle cellule
		{
		
		// inizializzazione delle componenti del flusso
		double fO2[] = {0., 0. ,0. };
		double fG[] = {0., 0. ,0. };
		double fA[] = {0., 0. ,0. };
		double fAcL[] = {0., 0. ,0. };
				
		// inizializzazione delle componenti del versore ambientale
		double versor_env[] = {0., 0. ,0. };

		// vettore posizione della cellula n-esima
		double rpos[3];
		rpos[0] = x[n];
		rpos[1] = y[n];
		rpos[2] = z[n];
		
		// concentrazioni locali extracellulari
		double concO2 = O2[n]/volume[n];
		double concG = G_extra[n]/volume_extra[n];
		double concA = A_extra[n]/volume_extra[n];
		double concAcL = AcL_extra[n]/volume_extra[n];
		
		
		for( int k=0; k< neigh[n]; k++)						// loop sulle cellule adiacenti alla cellula n-esima
			{
			// nome della k-esima cellula adiacente
			int neighbor = vneigh[n][k];	

			// vettore posizione della k-esima cellula adiacente
			double rposk[3];
			rposk[0] = x[neighbor];
			rposk[1] = y[neighbor];
			rposk[2] = z[neighbor];
			
			// versore da cellula n-esima a cellula k-esima
			double versornk[] = {0., 0. ,0. };
		
			// costruzione (parziale) del versore, la normalizzazione viene fatta piu' avanti
			double vsum = 0.;
			for(int kk=0; kk<3; kk++)	// loop sulle componenti
				{
				versornk[kk]=rposk[kk]-rpos[kk];	// componente kk-esima
				vsum += versornk[kk]*versornk[kk];	// calcolo del modulo del vettore
				}
			vsum = sqrt(vsum);						// modulo del vettore

			// concentrazioni locali extracellulari della cellula adiacente
			double concO2n = O2[neighbor]/volume[neighbor];
			double concGn = G_extra[neighbor]/volume_extra[neighbor];
			double concAn = A_extra[neighbor]/volume_extra[neighbor];
			double concAcLn = AcL_extra[neighbor]/volume_extra[neighbor];
			
			for(int kk=0; kk<3; kk++)	// loop sulle componenti
				{
				// normalizzazione del versore 
				versornk[kk] /= vsum;										
				
				// calcolo del versore ambientale
				versor_env[kk] += -vcsurf[n][k]*versornk[kk];					
							
				fO2[kk] += (double)(Diff_ES_O2*(concO2-concO2n)*gnk[n][k]*versornk[kk]);	// calcolo dei flussi
				fG[kk] += (double)(Diff_ES_G*(concG-concGn)*gnk[n][k]*versornk[kk]);
				fA[kk] += (double)(Diff_ES_A*(concA-concAn)*gnk[n][k]*versornk[kk]);
				fAcL[kk] += (double)(Diff_ES_AcL*(concAcL-concAcLn)*gnk[n][k]*versornk[kk]);
				}
			}
		
		// modulo del versore ambientale
		double vsum2 = sqrt( versor_env[0]*versor_env[0] + versor_env[1]*versor_env[1] + versor_env[2]*versor_env[2] );

		// ultimo loop per l'inserimento del termine ambientale nel flusso
		for(int kk=0; kk<3; kk++)
			{
			versor_env[kk] /= vsum2;
			fO2[kk] += (double)(Diff_W_O2*(concO2-concO2env)*g_env[n]*versor_env[kk]);	// calcolo dei flussi
			fG[kk] += (double)(Diff_W_G*(concG-concGenv)*g_env[n]*versor_env[kk]);
			fA[kk] += (double)(Diff_W_A*(concA-concAenv)*g_env[n]*versor_env[kk]);
			fAcL[kk] += (double)(Diff_W_AcL*(concAcL-concAcLenv)*g_env[n]*versor_env[kk]);
			}
		
		// printout
		flow_b_file.write( (char*)fO2, 3*sizeof(double) );
		flow_b_file.write( (char*)fG, 3*sizeof(double) );
		flow_b_file.write( (char*)fA, 3*sizeof(double) );
		flow_b_file.write( (char*)fAcL, 3*sizeof(double) );

		}
	
	// chiusura del file di output
	flow_b_file.close();
				
}
