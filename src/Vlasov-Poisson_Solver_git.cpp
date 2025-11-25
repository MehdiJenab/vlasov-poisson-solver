#include <iostream>
#include <cstdlib> 
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <math.h>
#include <cmath> // for abs(double)
#include <fstream> //writing and reading from files
#include <string>
#include <sstream> 
#include <mpi.h> //parallel programming
#include <vector>
#include <ctime>
#include <sys/time.h>
#include <sys/types.h> //to create folder in Contour
#include <sys/stat.h> //to create folder in Contour
#include <jansson.h> // to handle JSON file
#include "StorageClass.h" 
#include "config.h" //to handle preprocessing flags
#include <hdf5.h>
#include <map>
#include <algorithm> //for using max_element from standard library
#include <iomanip> //for std::setprecision
#include "ParametersClass.h" 
#include "MpiClass.h"
#include "PhasePointClass.h"
#include "SpeciesClass.h"
#include "GridClass.h"
#include "naming.h"

//command line: mpic++ -Wall vlasov.cpp -I ../eigen/ -I /usr/include/hdf5/openmpi/ -o vlasov -lz -ljansson
//mpic++ -Wall -I ../eigen/ -I /usr/include/hdf5/openmpi/ vlasov.cpp -lz -ljansson -lhdf5_openmpi -o vlasov
//mpic++ -Wall -I /usr/include/hdf5/openmpi/ vlasov_160402.cpp -lz -ljansson -lhdf5_openmpi -o vlasov

/*
 using std::cerr;
  using std::cout;
  using std::endl;
  using std::fixed;
  using std::flush;
  using std::ofstream;
  using std::scientific;
  using std::stringstream;
  using std::vector;
  using std::setprecision;
  using std:: max_element;
  using std:: min_element;*/
//using namespace std;


class TimingClass {
public:
	struct timeval timer1,timer2;
	ParametersClass * P;
	
	TimingClass(ParametersClass * P_in){P = P_in;}
	void get_initial_time() { gettimeofday(&timer1, 0);}
	void get_final_time  () { gettimeofday(&timer2, 0);}
	
	void print_timing    (const char* subroutine_name) {
		P->file_Info << " <"<<subroutine_name<<"= "<< fixed << setprecision(2) <<(timer2.tv_sec - timer1.tv_sec)+(timer2.tv_usec - timer1.tv_usec)/1e6<<">";
	}
  };

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
class SimulationClass {

public:  
	//(0) Number of species
	int N_species;
	
	//(1) objects which exists in simulation class
	//---- MPI variables ------
	MpiClass * MPI;
	int size, rank,NODE;
	int * Left_cpu;
	int * Right_cpu;
	GridClass * Grid; // one object
	SpeciesClass ** Species; // as many as N_species
	ParametersClass ** P; // as many as N_species 
	ofstream file_Time_Ene_Tot;
	
	//(2) the connection between phase space and spatial space
	vector<double> Den_g_All;
	double Ex;   
	
	//(3) species control flag
	int i_species;
	
	//(4) printing parameters
	vector<ofstream*>  file_Density;
	
	//(5) total stuff for each species  
	vector<int> Num_PhasePoints_Tot;
	double Ene_Tot, Ene_Tot_ini;
	
	// to control the number of HDF5 outputs
	int i_HDF5_contour, i_HDF5_spatial;

	// Timing 
	int rank_timing=0;

	// how to creat a 3D vector with [5,3,2]vector with number 4 in each node
	//<vector<vector<int> > > vec (5,vector<vector<int> >(3,vector <int>(2,4)));
	//vector <SpeciesClass> Species;
	
	// data structure of "leapfrog trapezoidal"
	int x , vx, vx0, vx1;
	double d_Time_half;

	//===============================================================================================
	SimulationClass(){

		MPI = new MpiClass();
		MPI->Initialize_MPI();
		size = MPI->size;
		rank = MPI->rank;
		NODE = MPI->NODE;
		Left_cpu = new int [size];
		Right_cpu = new int [size];
		for (int i=0; i<size; ++i){
			Left_cpu[i]     =  MPI->Left_cpu[i];
			Right_cpu[i]    =  MPI->Right_cpu[i];	
		}
		
		//(0) reading number of species
		Read_Num_Species ();
		
		//(1) creating N_species instances of ParametersClass
		P = new  ParametersClass * [N_species];
		for (i_species=0; i_species<=N_species; i_species++) { 
		P[i_species] = new ParametersClass(i_species, MPI);
		x = P[i_species]->x;
		vx=P[i_species]->vx;
		vx0=P[i_species]->vx0;
		vx1=P[i_species]->vx1;
		d_Time_half= P[i_species]->d_Time/2.0;
		}
		
		//(2) creating N_species instances of SpeciesClass      
		Species= new SpeciesClass * [N_species];
		for (i_species=0; i_species<=N_species; i_species++) { 
		Species[i_species] = new SpeciesClass(P[i_species], MPI);
		}
		
		//(3) creating one instance of GridClass      
		Grid = new GridClass(P[0], MPI);
		
		//(4) initializing of connection between phase space and spatial space
		Den_g_All.resize(P[0]->nX_cpu);
	};
	//===============================================================================================

	//***********************************************************************************************
	void Initialization_SimulationClass(){
		#ifdef debug_initial
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- Initial Step (0): initialization of simulation class");
		#endif  
		
		if (rank==0) Print_Time_Date();
		
		#ifdef Techplot_Usage
		Create_Folders_Techplot ();
		#endif
			
			//(5) creating files for output
		stringstream pathname;
		
		#ifdef print_Temporal_Techplot
		pathname.str("");
		pathname << "./temporal/Ene_Tot.dat";
		file_Time_Ene_Tot.open(pathname.str().c_str());
		#endif
			

		Grid->Initialization_GridClass();
		
		for (i_species=0; i_species<=N_species; i_species++) { 
			Species[i_species]->Initialization_SpeciesClass();
			}
		
		
		
		
		#ifdef debug_initial
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- Initial Step (0): initialization of variables");
		#endif
			Num_PhasePoints_Tot.resize(N_species+1); //Jenab, are you sure?

		//---------------------------------------------------
		for (i_species=0; i_species<=N_species; i_species++){
		// 	Species[i_species]->DF_INITIALIZATION();
			Species[i_species]->DF_INITIALIZATION_2();
			#ifdef debug_initial
			P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- Initial Step (1): DF initialization");
			#endif
			Species[i_species]->RELOCATION_CPU();
			#ifdef debug_initial
			P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- Initial Step (2): DF initialization");
			#endif
		
		}

	//---------------------------------------------------
	
	

		


	}
	//***********************************************************************************************
 
 
	//***********************************************************************************************
	void Convergance_Step ()  {
		vector<double> RO_diff_array(size);
		double RO_diff,RO_diff_Tot;
		
		int i_converage;
		
		
		
		Grid->RO_Zero();
		i_converage = 0;
		RO_diff_Tot = 100.0;
		while (RO_diff_Tot>0.00005) { 
			Grid-> PHI_Static(i_converage);
			i_converage += 1;
			//--- (1) Xpush ----
			Xpush_Simulation();      
			
			//--- (2) calculation of electric field ----
			Grid->RO_Old_Store();
			Grid->RO_Zero();
			//-----------------------------------------------------------------------
			for (i_species=0; i_species<=N_species; i_species++){
				fill(Den_g_All.begin(),Den_g_All.end(),0.0);
				Species[i_species]->Interpolate_Integrate (Den_g_All);
				if (P[i_species]->Boltzmann_Relation == 0) Grid->RO_Calculator(Den_g_All);
			}
			//-----------------------------------------------------------------------
			Grid->Ex_Zero();

			if (N_species == 0) Grid->Quasineutrality(P[0]->Charge);
			Grid->Compare_RO(RO_diff);
			
			
			MPI_Gather(&RO_diff, 1, MPI_DOUBLE, &RO_diff_array.at(0), 1, MPI_DOUBLE,0, MPI_COMM_WORLD);
			if (rank==0) {
				RO_diff_Tot= 0;
				for (int i=0; (i<=NODE); ++i){
					if (RO_diff_Tot < abs(RO_diff_array.at(i))) RO_diff_Tot = abs(RO_diff_array.at(i)) ;
					//cout <<Flag_Transfer_Tot<<endl;
				}
			}
			MPI_Bcast( &RO_diff_Tot, 1,  MPI_DOUBLE, 0, MPI_COMM_WORLD);   

			//--- (3) Vpush-----
			Grid->E_x_Calculator();
			Vpush_Simulation(P[0]->d_Time);
			if (rank == 0) cout<<"convergance iterations = "<<i_converage<<"  "<<RO_diff_Tot<<endl;
			if (i_converage <1000) RO_diff_Tot = 100.0;
			if (i_converage >1250) RO_diff_Tot = 0.0;
		}
	}
	//***********************************************************************************************
  
  
	//***********************************************************************************************
	void Initial_Step ()  {
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (0): Start of initial step");
		#endif
		//      #ifdef print_Contour_HDF5
		//          Write_Contour_HDF5();
		//      #endif
		//--- (1) calculation of electric field -------------------------------------------
		Set_Vx_interplotion(vx); Poisson_Solver_Simulation();
		//---------------------------------------------------------------------------------
		#ifdef print_Contour_HDF5
		Write_Contour_HDF5();
		#endif
		
		//--- (2) vpush half --------------------------------------------------------------
		//      Vpush_Simulation(d_Time_half);
		
		#ifdef print_Spatial_HDF5 
		Write_Spatial_HDF5();
		#endif
			
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (0): End of initial step");
		#endif
	}
	//***********************************************************************************************
	
	//***********************************************************************************************
	void Initial_Step_Leapfrog_Trapezoidal ()  {
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (0): Start of initial step");
		#endif
		//      #ifdef print_Contour_HDF5
		//          Write_Contour_HDF5();
		//      #endif
		//--- (1) calculation of electric field -------------------------------------------
		Set_Vx_interplotion(vx); Poisson_Solver_Simulation();
		Grid->Save_E_in_E0();
		//---------------------------------------------------------------------------------
		#ifdef print_Contour_HDF5
		Write_Contour_HDF5();
		#endif
		
		//--- (2) vpush half --------------------------------------------------------------
		Xpush_Simulation_new(x,x,vx,d_Time_half);
		
		#ifdef print_Spatial_HDF5 
		Write_Spatial_HDF5();
		#endif
			
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (0): End of initial step");
		#endif
	}
	//***********************************************************************************************
	//*********************************************************************
	void Do_Step () {		
		Xpush_Simulation();
		#ifdef Co_moving_frame
		for (i_species=0; i_species<=N_species; i_species++){
			if ((P[0]->i_Time)%(P[0]->nTime_Co_Moving) == 0){
				Species[i_species]->Moving_box_simulation(P[i_species], P[0]->d_Time);
			}
		}
		#endif
		Set_Vx_interplotion(vx); Poisson_Solver_Simulation();
		Vpush_Simulation(P[0]->d_Time);	
	}
	//*********************************************************************

	//*********************************************************************
	void Do_Step_Leapfrog_trapezoidal () {
		
		Set_Vx_interplotion(vx); Poisson_Solver_Simulation();
		
		Vpush_Simulation_new(vx0,vx,P[0]->d_Time);
		
		Grid->Average_E_E0_store_in_E();
		Vpush_Simulation_new(vx,vx,d_Time_half);
		
		Simulation_average_V_V0(vx1, vx, vx0);
		Xpush_Simulation_new(x, x, vx1,d_Time_half);
		
		Set_Vx_interplotion(vx0); Poisson_Solver_Simulation();
		
		Vpush_Simulation_new(vx0,vx1,P[0]->d_Time);
		
		Grid->Average_E_E0_store_in_E();
		Vpush_Simulation_new(vx,vx,d_Time_half);
		
		Simulation_average_V_V0(vx1, vx, vx0);
		Xpush_Simulation_new(x, x, vx1,d_Time_half);
		
		#ifdef Co_moving_frame
		for (i_species=0; i_species<=N_species; i_species++){
			if ((P[0]->i_Time)%(P[0]->nTime_Co_Moving) == 0){
				Species[i_species]->Moving_box_simulation(P[i_species], P[0]->d_Time);
			}
		}
		#endif
	
	}
	//*********************************************************************
	
	
	//*********************************************************************
	void Do_Step_Eulerian_trapezoidal () {
		Save_X0_Vx0_Simulation();
		Set_Vx_interplotion(vx); Poisson_Solver_Simulation();
		
		Vpush_Simulation(P[0]->d_Time);
		
		Xpush_Simulation();      
		
		Set_Vx_interplotion(vx); Poisson_Solver_Simulation();
		Vpush_Simulation_over_Vx0(P[0]->d_Time);
		
		Simulation_Calculate_Vxhalf_Store_in_Vx0();
		Xpush_Simulation_over_X0_by_Vxhalf();
		#ifdef Co_moving_frame
		for (i_species=0; i_species<=N_species; i_species++){
			if ((P[0]->i_Time)%(P[0]->nTime_Co_Moving) == 0){
				Species[i_species]->Moving_box_simulation(P[i_species], P[0]->d_Time);
			}
		}
		#endif
	
	}
	//*********************************************************************
	
	//*********************************************************************
	void Save_X0_Vx0_Simulation(){
		for (i_species=0; i_species<=N_species; i_species++){
			for(Storage<PhasePointClass>::Iterator particle_iter = Species[i_species]->PhasePoints->begin(); particle_iter != Species[i_species]->PhasePoints->end(); ++particle_iter) {
				(*particle_iter).Save_X0_Vx0();
			}
		}
	}
	//*********************************************************************
	//*********************************************************************
	void Simulation_Calculate_Vxhalf_Store_in_Vx0(){
		#if flag_timing_of_run > 0
		TimingClass Timer(P[0]);
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
		
		
		//-----------------------------------------------------------------------
		for (i_species=0; i_species<=N_species; i_species++){
			for(Storage<PhasePointClass>::Iterator particle_iter = Species[i_species]->PhasePoints->begin(); particle_iter != Species[i_species]->PhasePoints->end(); ++particle_iter) {
   				(*particle_iter).Calculate_Vxhalf_Store_in_Vx0();
			}
		}
		//-----------------------------------------------------------------------
		
		
		#if flag_timing_of_run>0	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("Store");
		#endif
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (1): after storing");
		#endif
	}
	//*********************************************************************
	//*********************************************************************
	void Xpush_Simulation_over_X0_by_Vxhalf(){
		#if flag_timing_of_run==1
		TimingClass Timer(P[0]);
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
		
		
		//-----------------------------------------------------------------------
		for (i_species=0; i_species<=N_species; i_species++){
			for(Storage<PhasePointClass>::Iterator particle_iter = Species[i_species]->PhasePoints->begin(); particle_iter != Species[i_species]->PhasePoints->end(); ++particle_iter) {
   				(*particle_iter).Xpush_over_X0_by_Vxhalf(P[i_species], P[0]->d_Time);
			}
		}
		//-----------------------------------------------------------------------
		
		
		#if flag_timing_of_run==1	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("Xpush");
		#endif
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (1): after Xpush of Do Step");
		#endif
		
		#if flag_timing_of_run==1
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
		
		//-----------------------------------------------------------------------
		for (i_species=0; i_species<=N_species; i_species++){
			Species[i_species]->RELOCATION_CPU();
		}
		//-----------------------------------------------------------------------
		
		#if flag_timing_of_run==1	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("Xrel");
		#endif
		
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (2): after Relocation of Do Step");
		#endif      
	}
	//*********************************************************************
  
  
	//*********************************************************************
	void Xpush_Simulation(){
		#if flag_timing_of_run==1
		TimingClass Timer(P[0]);
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
		
		
		//-----------------------------------------------------------------------
		for (i_species=0; i_species<=N_species; i_species++){
			for(Storage<PhasePointClass>::Iterator particle_iter = Species[i_species]->PhasePoints->begin(); particle_iter != Species[i_species]->PhasePoints->end(); ++particle_iter) {
   				(*particle_iter).Xpush_trapezoidal(P[i_species], P[0]->d_Time);
//   				(*particle_iter).Xpush(P[i_species], P[0]->d_Time);
			}
		}
		//-----------------------------------------------------------------------
		
		
		#if flag_timing_of_run==1	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("Xpush");
		#endif
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (1): after Xpush of Do Step");
		#endif
		
		#if flag_timing_of_run==1
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
		
		//-----------------------------------------------------------------------
		for (i_species=0; i_species<=N_species; i_species++){
			Species[i_species]->RELOCATION_CPU();
		}
		//-----------------------------------------------------------------------
		
		#if flag_timing_of_run==1	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("Xrel");
		#endif
		
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (2): after Relocation of Do Step");
		#endif      
	}
	//*********************************************************************
	
	//*********************************************************************
	void Xpush_Simulation_new(int P_x_out, int P_x_in, int P_v_in,double & dt_in){
		#if flag_timing_of_run==1
		TimingClass Timer(P[0]);
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
		
		
		//-----------------------------------------------------------------------
		for (i_species=0; i_species<=N_species; i_species++){
			for(Storage<PhasePointClass>::Iterator particle_iter = Species[i_species]->PhasePoints->begin(); particle_iter != Species[i_species]->PhasePoints->end(); ++particle_iter) {
   				(*particle_iter).Xpush_new(P_x_out, P_x_in, P_v_in, dt_in, P[i_species]);
			}
		}
		//-----------------------------------------------------------------------
		
		
		#if flag_timing_of_run==1	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("Xpush");
		#endif
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (1): after Xpush of Do Step");
		#endif
		
		#if flag_timing_of_run==1
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
		
		//-----------------------------------------------------------------------
		for (i_species=0; i_species<=N_species; i_species++){
			Species[i_species]->RELOCATION_CPU();
		}
		//-----------------------------------------------------------------------
		
		#if flag_timing_of_run==1	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("Xrel");
		#endif
		
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (2): after Relocation of Do Step");
		#endif      
	}
		
		
		
		
	
	//*********************************************************************
	//*********************************************************************
	void Simulation_average_V_V0(int& P_vx_out, int& P_v_in_1, int& P_v_in_2){
		#if flag_timing_of_run==1
		TimingClass Timer(P[0]);
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
		
		//-----------------------------------------------------------------------
		for (i_species=0; i_species<=N_species; i_species++){
			for(Storage<PhasePointClass>::Iterator  particle_iter  = Species[i_species]->PhasePoints->begin(); particle_iter != Species[i_species]->PhasePoints->end(); ++particle_iter) {
   				(*particle_iter).average_V_V0(P_vx_out, P_v_in_1, P_v_in_2 ,P[i_species]);
			}
		}
		//-----------------------------------------------------------------------  
		#if flag_timing_of_run==1	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("average V");
		#endif
		
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (4): after average V of Do Step");
		#endif
	}
	//*********************************************************************
	
	

	//*********************************************************************
	void Vpush_Simulation_new(int P_vx_out, int P_vx_in, double & dt_in){
		#if flag_timing_of_run==1
		TimingClass Timer(P[0]);
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
		
		//-----------------------------------------------------------------------
		Grid->Ex_Ghost_Talk();      
		for (i_species=0; i_species<=N_species; i_species++){
			for(Storage<PhasePointClass>::Iterator  particle_iter  = Species[i_species]->PhasePoints->begin(); particle_iter != Species[i_species]->PhasePoints->end(); ++particle_iter) {
				Grid->find_Ep((*particle_iter).Data_PP->at(P[i_species]->x), Ex);
   				(*particle_iter).Vpush_new( P_vx_out, P_vx_in, P[i_species], particle_iter,dt_in,Ex);
			}
		}
		//-----------------------------------------------------------------------  
		#if flag_timing_of_run==1	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("Vpush");
		#endif
		
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (4): after vpush of Do Step");
		#endif
	}
	//*********************************************************************
	
	
	//*********************************************************************
	void Vpush_Simulation(double & d_Time_S){
		#if flag_timing_of_run==1
		TimingClass Timer(P[0]);
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
		
		//-----------------------------------------------------------------------
		Grid->Ex_Ghost_Talk();      
		for (i_species=0; i_species<=N_species; i_species++){
			for(Storage<PhasePointClass>::Iterator  particle_iter  = Species[i_species]->PhasePoints->begin(); particle_iter != Species[i_species]->PhasePoints->end(); ++particle_iter) {
				Grid->find_Ep((*particle_iter).Data_PP->at(P[i_species]->x), Ex);
   				(*particle_iter).Vpush_trapezoidal( P[i_species], particle_iter,d_Time_S,Ex);
//   				(*particle_iter).Vpush( P[i_species], particle_iter,d_Time_S,Ex);
			}
		}
		//-----------------------------------------------------------------------  
		#if flag_timing_of_run==1	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("Vpush");
		#endif
		
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (4): after vpush of Do Step");
		#endif
	}
	//*********************************************************************
	
	//*********************************************************************
	void Vpush_Simulation_over_Vx0(double & d_Time_S){
	
		#if flag_timing_of_run==1
		TimingClass Timer(P[0]);
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
		
		
		//-----------------------------------------------------------------------
		Grid->Ex_Ghost_Talk();      
		for (i_species=0; i_species<=N_species; i_species++){
			for(Storage<PhasePointClass>::Iterator  particle_iter  = Species[i_species]->PhasePoints->begin(); particle_iter != Species[i_species]->PhasePoints->end(); ++particle_iter) {
				Grid->find_Ep((*particle_iter).Data_PP->at(P[i_species]->x), Ex);
				(*particle_iter).Vpush_over_Vx0( P[i_species], particle_iter,d_Time_S,Ex);
			}
		}
		//-----------------------------------------------------------------------  
			
		#if flag_timing_of_run==1	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("Vpush_Vx0");
		#endif
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (4): after vpush of Do Step");
		#endif
	}
	//*********************************************************************

	//*********************************************************************
	void Set_Vx_interplotion(int & i_vx_in){
		for (i_species=0; i_species<=N_species; i_species++){
			P[i_species]->vx_interpol = i_vx_in;
		}
		
	}
  
	//*********************************************************************
	void Poisson_Solver_Simulation(){
		
		#if flag_timing_of_run==1
		TimingClass Timer(P[0]);
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
			
		Grid->RO_Zero();
		
		
		#ifdef debug_Poisson
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (2): Start of Poisson Solver Simulation");
		#endif 
		//-----------------------------------------------------------------------
		for (i_species=0; i_species<=N_species; i_species++){
			fill(Den_g_All.begin(),Den_g_All.end(),0.0);
			Species[i_species]->Interpolate_Integrate (Den_g_All);
			if (P[i_species]->Boltzmann_Relation == 0) Grid->RO_Calculator(Den_g_All);
		}
		//        Grid->RO_Zero_on_borders();
		//-----------------------------------------------------------------------
		#if flag_timing_of_run==1	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("IngInp");
		#endif
		
		#ifdef debug_Poisson
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (2): Start of Poisson Solver ");
		#endif 
		
		#if flag_timing_of_run==1
		if (rank==rank_timing) Timer.get_initial_time();		
		#endif
		
		
		//-----------------------------------------------------------------------
		Grid->Ex_Zero();

		if (N_species == 0) Grid->Quasineutrality(P[0]->Charge);
		Grid->Poisson_Solver();
		
		//-----------------------------------------------------------------------
		#if flag_timing_of_run==1	
		if (rank==rank_timing) Timer.get_final_time();
		if (rank==rank_timing) Timer.print_timing("Poisson");
		#endif
		
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (3): after electric field");
		#endif 
	}
	//*********************************************************************


  

  
  
	//***********************************************************************************************
	template<class T>
	void RecursiveCommas(std::ostream& os, T n){

		T rest = n % 1000; //"last 3 digits"
		n /= 1000;         //"begining"

		if (n > 0) {
			RecursiveCommas(os, n); //printing "begining"
			//and last chunk
			os << ',' << std::setfill('0') << std::setw(3) << rest;
		}
		else
			os << rest; //first chunk of the number
	}

	//***********************************************************************************************



	//***********************************************************************************************
	void Run_Simulation(){
		#ifdef debug
		P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (0): Start Running the simulation");
		#endif 
		if (rank==0) {P[0]->file_Info << fixed << setprecision(2);}

		TimingClass Time_run (P[0]);
		#if flag_timing_of_run==1
		TimingClass Time_momentum (P[0]);
		#endif
		

		if (rank==rank_timing) Time_run.get_initial_time();
			
			
		P[0]->i_Time=0;
		
		#ifdef print_Spatial_HDF5
		i_HDF5_spatial = 1;
		Create_File_HDF5_Spatial(i_HDF5_spatial);
		#endif
			
			
		#ifdef print_Contour_HDF5	
		i_HDF5_contour = 1;
		Create_File_HDF5_Contour(i_HDF5_contour);
		#endif
		
			
		//+++++++++++++++++++++++++++++++++++++++++++++++++
		if (P[0]->convergance_step == 1) Convergance_Step();
		//+++++++++++++++++++++++++++++++++++++++++++++++++
			
		//+++++++++++++++++++++
// 		Initial_Step();
		Initial_Step_Leapfrog_Trapezoidal();
		//+++++++++++++++++++++

		if (rank==rank_timing) Time_run.get_final_time();
		
		Number_PhasePoints();      
		if (rank==0) {
		P[0]->file_Info << "Initial Step"<< "   ";
			for (i_species=0; i_species<=N_species; i_species++){
				P[0]->file_Info << "number of "<<P[i_species]->Name <<" = " ; RecursiveCommas(P[0]->file_Info, Num_PhasePoints_Tot[i_species]);P[0]->file_Info << "   ";
			}
		}
		if (rank==rank_timing) Time_run.print_timing("initial loop");
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (rank==0) {P[0]->file_Info<<endl;}
		
		///////////////////////////////////////////////////////////////////////////////////
		for (P[0]->i_Time = 1; P[0]->i_Time<=P[0]->N_Time; P[0]->i_Time++) { 
			if (rank==rank_timing) Time_run.get_initial_time();
			
			#ifdef long_time_run
			if ((P[0]->i_Time)%(P[0]->N_Time/10) == 0){
				if (P[0]->i_Time != P[0]->N_Time){
					#ifdef print_Spatial_HDF5
					Close_File_HDF5_Spatial();
					i_HDF5_spatial += 1;
					Create_File_HDF5_Spatial(i_HDF5_spatial);
					#endif
				
					#ifdef print_Contour_HDF5
					Close_File_HDF5_Contour();
					i_HDF5_contour += 1;
					Create_File_HDF5_Contour(i_HDF5_contour);
					#endif
				}
			}
			#endif
			
			for (i_species=1; i_species<=N_species; i_species++){P[i_species]->i_Time = P[0]->i_Time;}
			for (i_species=0; i_species<=N_species; i_species++){P[i_species]->R_Time = P[i_species]->i_Time*P[i_species]->d_Time;}
			
			if (rank==0) {P[0]->file_Info << "|Time = " << fixed <<P[0]->R_Time << "  | ";}
			
			#ifdef debug
			P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (0): begining of Do Step");
			#endif  
			//+++++++++++++++++++++
			//Do_Step();
			//Do_Step_Eulerian_trapezoidal();
			Do_Step_Leapfrog_trapezoidal();
			//+++++++++++++++++++++

			#if flag_timing_of_run==1			
			if (rank==rank_timing) Time_momentum.get_initial_time();
			#endif
			
			#ifdef debug
			P[0]->Print_Debug(P[0]->line_NO=__LINE__,P[0]->Debug_mes="--- (1): end of Do Step");
			#endif 
			if (rank==0) Ene_Tot = 0;
			
			#ifdef print_Temporal_Techplot
			for (i_species=0; i_species<=N_species; i_species++){
				Species[i_species]->Entropy_F2_Species();
				Species[i_species]->Entropy_FlnF_Species();
				Species[i_species]->Ene_kin_Tot = 0.0;
				if (P[i_species]->Boltzmann_Relation == 0 ) Species[i_species]->Energy_kin_Species();
				if (rank==0) Ene_Tot += Species[i_species]->Ene_kin_Tot;
				if (rank==0) Species[i_species]->Write_Entropy_Temporal_Techplot();
			}
			#endif
			
			if (rank==0) Ene_Tot += Grid->Ene_ele_Tot;
			if (rank==0) file_Time_Ene_Tot <<P[0]->R_Time<<"   "<<Ene_Tot<<endl;
			Number_PhasePoints();
			#if flag_timing_of_run==1
			if (rank==rank_timing) Time_momentum.get_final_time();
			if (rank==rank_timing) Time_momentum.print_timing("Moment");
			#endif
			if (rank==0) {
				for (i_species=0; i_species<=N_species; i_species++){
					P[0]->file_Info << " "<< P[i_species]->Name <<" = "; RecursiveCommas(P[0]->file_Info, Num_PhasePoints_Tot[i_species]);P[0]->file_Info << "| ";
				}
			}
			if (rank==rank_timing) Time_run.get_final_time();
			if (rank==rank_timing) Time_run.print_timing("loop");
			
			#ifdef print_Temporal_Techplot
			if (size>=4){
				if (rank%(size/4)==0) Grid->Write_Temporal_Ex_Field_Techplot();
			}else{
				Grid->Write_Temporal_Ex_Field_Techplot();
			}
			if (rank==0) Grid->Write_Temporal_Ene_ele_Techplot();
			#endif 
			
			#ifdef print_Contour_HDF5
			if (P[0]->i_Time%P[0]->i_Contour==0) Write_Contour_HDF5();
			#endif
			
			#ifdef print_Spatial_HDF5 
			if (P[0]->i_Time%P[0]->i_Spatial_Output==0) Write_Spatial_HDF5();  
			#endif
			if (rank==0) {P[0]->file_Info<<endl;}
		}
		///////////////////////////////////////////////////////////////////////////////////


		#ifdef print_Spatial_HDF5
		Close_File_HDF5_Spatial();
		#endif
		
		#ifdef print_Contour_HDF5
		Close_File_HDF5_Contour();
		#endif
	}
	//*********************************************************************
    
	//*********************************************************************
	void Write_Spatial_HDF5() {
		Create_Group_HDF5_Spatial();
		Grid->Write_Dset_RO_Ex_PHI_HDF5();
		for (i_species=0; i_species<=N_species; i_species++){
			Species[i_species]->Write_Dset_Num_Den_HDF5(P[0]->Group_ID_S);
			#ifdef print_energy_kinetic_in_X_dimesion
			Species[i_species]->Write_Dset_Ene_kin_HDF5(P[0]->Group_ID_S);
			#endif    
		}
		Close_Group_HDF5_Spatial();
	}       
	//*********************************************************************
    
    
    
	//*********************************************************************
	void Create_Group_HDF5_Spatial() {
		// 1) creating the group----------------------------------------------------
		stringstream Group_Name_X;
		Group_Name_X << P[0]->i_Time;
		P[0]->Group_ID_S = H5Gcreate(P[0]->File_ID_S, Group_Name_X.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
		//-----------------------------------------------------------------------------------------------
		void Close_Group_HDF5_Spatial() {
		// 1) close the group----------------------------------------------------
		P[0]->status_hdf5 = H5Gclose(P[0]->Group_ID_S);
	}
	//*********************************************************************
    
	//*********************************************************************
	void Write_Contour_HDF5() {
				Create_Group_Contour_HDF5();
				for (i_species=0; i_species<=N_species; i_species++){
				#ifdef print_FPP
					Species[i_species]->Write_Dset_Contour_HDF5(P[0]->Group_ID_DF);
				#endif	
				#ifdef print_EPP
					Species[i_species]->Write_Dset_EPP_HDF5(P[0]->Group_ID_DF);
				#endif
				#ifdef print_MPP
					Species[i_species]->Write_Dset_MPP_HDF5(P[0]->Group_ID_DF);
				#endif
				//   	  if ((P[0]->i_Time)%(P[0]->N_Time/20) == 0){
	// 			if (P[0]->i_Time == P[0]->N_Time ){
	// 				Species[i_species]->PhasePoints_Sorting_Storing(P[0]->Group_ID_DF);}
		
	}
	Close_Group_HDF5_Contour();
	}
	//*********************************************************************

	//*********************************************************************
	void Create_Group_Contour_HDF5() {
		stringstream Group_Name_DF;
		Group_Name_DF << "Timestep_"<<P[0]->i_Time;
		P[0]->Group_ID_DF = H5Gcreate(P[0]->File_ID_DF, Group_Name_DF.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
	
	//---------------------------------------------------------------------
	void Close_Group_HDF5_Contour() {
		P[0]->status_hdf5 = H5Gclose(P[0]->Group_ID_DF);
	}
	//*********************************************************************
  
   
	//***********************************************************************************************
	void Create_File_HDF5_Spatial(int & i_HDF5_spatial){
		#define H5FILE_NAME_Spatial_1     "Spatial_1.h5"
		#define H5FILE_NAME_Spatial_2     "Spatial_2.h5"
		#define H5FILE_NAME_Spatial_3     "Spatial_3.h5"
		#define H5FILE_NAME_Spatial_4     "Spatial_4.h5"
		#define H5FILE_NAME_Spatial_5     "Spatial_5.h5"
		#define H5FILE_NAME_Spatial_6     "Spatial_6.h5"
		#define H5FILE_NAME_Spatial_7     "Spatial_7.h5"
		#define H5FILE_NAME_Spatial_8     "Spatial_8.h5"
		#define H5FILE_NAME_Spatial_9     "Spatial_9.h5"
		#define H5FILE_NAME_Spatial_10    "Spatial_10.h5"
		
		
		
		
		MPI_Info info  = MPI_INFO_NULL;  
		

		// Set up file access property list with parallel I/O access     
		P[0]->PList_ID = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(P[0]->PList_ID, MPI_COMM_WORLD, info);

		
		//Create a new file collectively and release property list identifier     
		if (i_HDF5_spatial == 1){
		P[0]->File_ID_S = H5Fcreate(H5FILE_NAME_Spatial_1, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_spatial == 2){
		P[0]->File_ID_S = H5Fcreate(H5FILE_NAME_Spatial_2, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_spatial == 3){
		P[0]->File_ID_S = H5Fcreate(H5FILE_NAME_Spatial_3, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_spatial == 4){
		P[0]->File_ID_S = H5Fcreate(H5FILE_NAME_Spatial_4, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_spatial == 5){
		P[0]->File_ID_S = H5Fcreate(H5FILE_NAME_Spatial_5, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_spatial == 6){
		P[0]->File_ID_S = H5Fcreate(H5FILE_NAME_Spatial_6, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_spatial == 7){
		P[0]->File_ID_S = H5Fcreate(H5FILE_NAME_Spatial_7, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_spatial == 8){
		P[0]->File_ID_S = H5Fcreate(H5FILE_NAME_Spatial_8, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_spatial == 9){
		P[0]->File_ID_S = H5Fcreate(H5FILE_NAME_Spatial_9, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_spatial == 10){
		P[0]->File_ID_S = H5Fcreate(H5FILE_NAME_Spatial_10, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}
		
		
		
		//     P[0]->File_ID_S = H5Fcreate(H5FILE_NAME_Spatial, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		H5Pclose(P[0]->PList_ID);    
	}
	//-----------------------------------------------------------------------------------------------
	void Close_File_HDF5_Spatial(){
		H5Fclose(P[0]->File_ID_S);
	}
	//***********************************************************************************************

	//***********************************************************************************************
	void Create_File_HDF5_Contour(int & i_HDF5_contour){
		
		//     #define H5FILE_NAME_contour     "Contour.h5"
		#define H5FILE_NAME_contour_1     "Contour_1.h5"
		#define H5FILE_NAME_contour_2     "Contour_2.h5"
		#define H5FILE_NAME_contour_3     "Contour_3.h5"
		#define H5FILE_NAME_contour_4     "Contour_4.h5"
		#define H5FILE_NAME_contour_5     "Contour_5.h5"
		#define H5FILE_NAME_contour_6     "Contour_6.h5"
		#define H5FILE_NAME_contour_7     "Contour_7.h5"
		#define H5FILE_NAME_contour_8     "Contour_8.h5"
		#define H5FILE_NAME_contour_9     "Contour_9.h5"
		#define H5FILE_NAME_contour_10    "Contour_10.h5"
		
		MPI_Info info  = MPI_INFO_NULL;

		// Set up file access property list with parallel I/O access     
		P[0]->PList_ID = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(P[0]->PList_ID, MPI_COMM_WORLD, info);

		
		//Create a new file collectively and release property list identifier     
		if (i_HDF5_contour == 1){
		P[0]->File_ID_DF = H5Fcreate(H5FILE_NAME_contour_1, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_contour == 2){
		P[0]->File_ID_DF = H5Fcreate(H5FILE_NAME_contour_2, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_contour == 3){
		P[0]->File_ID_DF = H5Fcreate(H5FILE_NAME_contour_3, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_contour == 4){
		P[0]->File_ID_DF = H5Fcreate(H5FILE_NAME_contour_4, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_contour == 5){
		P[0]->File_ID_DF = H5Fcreate(H5FILE_NAME_contour_5, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_contour == 6){
		P[0]->File_ID_DF = H5Fcreate(H5FILE_NAME_contour_6, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_contour == 7){
		P[0]->File_ID_DF = H5Fcreate(H5FILE_NAME_contour_7, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_contour == 8){
		P[0]->File_ID_DF = H5Fcreate(H5FILE_NAME_contour_8, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_contour == 9){
		P[0]->File_ID_DF = H5Fcreate(H5FILE_NAME_contour_9, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}else if (i_HDF5_contour == 10){
		P[0]->File_ID_DF = H5Fcreate(H5FILE_NAME_contour_10, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		}
		
		
		
		//     P[0]->File_ID_DF = H5Fcreate(H5FILE_NAME_contour, H5F_ACC_TRUNC, H5P_DEFAULT, P[0]->PList_ID);
		H5Pclose(P[0]->PList_ID); 
	
	}
	//------------------------------------------------------------------------------------------------
	void Close_File_HDF5_Contour(){
		H5Fclose(P[0]->File_ID_DF);
	}
	//***********************************************************************************************
  
	//***********************************************************************************************
	void Number_PhasePoints (){
		vector<double> Num_PhasePoints_array(size);
		for (i_species=0; i_species<=N_species; i_species++){
			double Num_PhasePoints;			
			Num_PhasePoints = Species[i_species]->PhasePoints->getFill();			
			//-----------------------------------------------------------------------------------------------
			MPI_Gather(&Num_PhasePoints, 1, MPI_DOUBLE,&Num_PhasePoints_array[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank==0) {
				Num_PhasePoints_Tot[i_species]=0;
				for (int i=0; (i<=NODE); ++i){
					Num_PhasePoints_Tot[i_species] +=  Num_PhasePoints_array[i];
				} 
			}
		}
	}
	//***********************************************************************************************
  
	//***********************************************************************************************
	void Print_Time_Date() {
		time_t rawtime;
		struct tm * timeinfo;
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		cout<< endl;
		cout<< endl;
		cout<< "------------------------------------------------------"<<endl;
		printf ( "The current date/time is: %s", asctime (timeinfo) ); 
		cout<< "------------------------------------------------------"<<endl;
		cout<< endl;
		cout<< endl;
	}
	//***********************************************************************************************

	//***********************************************************************************************
	void Create_Folders_Techplot() {      
		#ifdef print_Temporal_Techplot
		mkdir("./temporal", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		#endif
	}
	//***********************************************************************************************
       

	//***********************************************************************************************
	void Read_Num_Species (){
	
		size_t flags = 0;
		json_error_t error;
		//       json_t* file_json, config_json,N_species_json;
		
		if (rank==0){cout << "-----------------------------------------------------"<<endl;}
		
		//--- loading the JSON file ------------------------------------------
		json_t* file_json = json_load_file("./config.json", flags, &error);
		if(!file_json){
			cerr << "in line " << error.line << ": " << error.text << endl;
			exit(1);
		}
		//--------------------------------------------------------------------
		
		//--- returning the config object ------------------------------------  
		json_t* config_json = json_object_get(file_json, "config");
		if(!config_json){
			cerr << "in line " << error.line << ": " << error.text << endl;
			exit(1);
		}
		//--------------------------------------------------------------------
		
		
		//--- returning the N_species object --------------------------------  
		json_t* N_species_json = json_array_get(config_json, 0);
		if (!N_species_json) cout<<"can't catch 0 member of config=N_species" <<endl;
		N_species_json = json_object_get(N_species_json, "N_species");
		if (!N_species_json) cout<<"can't catch N_species object" <<endl;
		N_species = json_integer_value(N_species_json);
		if (rank==0) cout << "Number of species " <<" = " << N_species << endl;
		//--------------------------------------------------------------------    
	}
	//***********************************************************************************************
  
  
	void Finalize_MPI(){
		MPI->Finalize_MPI();
	} 
	
  };  
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  

  
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int main () {
  

  // --------------------------------------------------------------------------
  SimulationClass SimulationObject;
  SimulationObject.Initialization_SimulationClass();
  SimulationObject.Run_Simulation();
  

 
  // --------------------------------------------------------------------------
  // Finalize the MPI environment.
  SimulationObject.Finalize_MPI();
  
  // --------------------------------------------------------------------------
  //MPI_Barrier(MPI_COMM_WORLD);
  //cout << "===============================================================" <<rank <<P[0]->Name<<endl;
  //exit(0);
  
  return 0;
  }
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
