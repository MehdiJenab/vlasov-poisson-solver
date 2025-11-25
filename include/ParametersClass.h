#ifndef ParametersClass_Included
#define ParametersClass_Included


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
#include "MpiClass.h"
#include "naming.h"


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
class ParametersClass {

public:
  
	//1) grid information 
		//1-1) spatial direction
		//1-1-1) grid on spatial direction
			int nX;
			int nX_cpu;
		//1-1-2) box on spatial direction
			double X_min_Tot,X_max_Tot,X_lenght;
			double X_min,X_max;
			double dX;
			int iX_min;
		
		//1-2) velocity direction 
		//1-2-1) grid on velocity direction
			int nVx;
		//1-2-2) box on velocity direction
			double Vx_min, Vx_max, velocity_box_simulation;
			double dVx;
			double Kappa_DF;
			int DF_X_NMRE,DF_Vx_MK, DF_NotShifted_Shifted, Boltzmann_Relation;
			int iX_switching_EPP;
			int nTime_Co_Moving;
	

	//2) distribution fucntion information
		//2-1) general information
			char const * Name; double Mass, Temp, Charge, norm_factor;
		//2-2) perturbation information 	
			double Alpha; int NoW; 
			double  Kappa;
		//2-3) jet information
			double Den_Ratio_b;
		//2-4) Quasineutrality information 
			double Density;
		//2-4) per cell information
			int pVx, pX;
		
			int number_of_data;
			int size_shadow_Vx;
		

	//3) Time grid
		double R_Time; int i_Time;
		double d_Time; int N_Time;
		int Num_Spatial_Output, Num_Contour;
		int i_Spatial_Output, i_Contour;
		int Plasma_approximation; // to use Boltzmann Relation
		int convergance_step; // for introducing solitons to the code
	//4) debugging
		int line_NO; //line number of the debugging 
		const char* Debug_mes; //message passed to terminal by debug
		ofstream file_Info;
	
	// data structure 
		/*
		int fp = 0, x = 1, x0 = 2, vx=3, vx0=4;
		#ifdef print_EPP
		int epp=5; 
		#endif
		#ifdef print_MPP
		int mpp=6; 
		#endif
		*/
	
	// data structure of "leapfrog trapezoidal"
		int fp = 0, x = 1, vx=2, vx0=3, vx1=4;
		#ifdef print_EPP
		int epp=5; 
		#endif
		#ifdef print_MPP
		int mpp=6; 
		#endif
		int vx_interpol = vx;
	//5) HDF5 parameters --------------------------------------------------
		hid_t File_ID_DF,File_ID_S;
		hid_t Group_ID_DF,Group_ID_S;
		hid_t PList_ID;      // Property list identifier
		
		hid_t DSet_ID_DF, DSet_ID_DFp_x, DSet_ID_DFp_Num, DSet_ID_E, DSet_ID_PHI, DSet_ID_RO, DSet_ID_Num_Den, DSet_ID_Ene_kin;      // Dataset identifier
		hid_t Space_ID_DF, Space_ID_DFp_x, Space_ID_DFp_Num,Space_ID_E, Space_ID_PHI, Space_ID_RO, Space_ID_Num_Den, Space_ID_Ene_kin;  // Dataspace identifier in file
		hid_t Memory_ID_DF, Memory_ID_DFp_x, Memory_ID_DFp_Num, Memory_ID_E, Memory_ID_PHI, Memory_ID_RO, Memory_ID_Num_Den, Memory_ID_Ene_kin; // Dataspace identifier in memory
		
		hsize_t Dim_X[1], Dim_X_Chunk[1]; //Dataset dimensions
		hsize_t Dim_X_Vx[2],Dim_X_Vx_Chunk[2]; //Dim_X_Vx = (/nX_Tot,nVx_g/), Dim_X_Vx_Chunk = (/nX_cpu,nVx_g/)
		
		
		hsize_t Count_X_Vx[2], OffSet_X_Vx[2],Stride_X_Vx[2], Block_X_Vx[2];
		hsize_t Count_X[1], OffSet_X[1],Stride_X[1], Block_X[1];
		
		
		hsize_t Dim_3D[3], Dim_3D_Chunk[3]; //Dataset dimensions
		hsize_t Count_3D[3], OffSet_3D[3],Stride_3D[3], Block_3D[3];

		int Rank_3D;
		int Rank_DFp;
		int Rank_X_Vx,Rank_X;
		herr_t	status_hdf5;
	//---------------------------------------------------------------------
	//---- MPI variables ------
		MpiClass * MPI;
		int size, rank,NODE;
		int * Left_cpu;
		int * Right_cpu;
	//---------------------------------------------------------------------
		
	vector<int> VxUp;
	double minus_2,minus_1,ghost_nX_cpu,nX_cpu_plus_1; // for controling ghosts cells
		
	int multigrid_call_num, size_MG, Iteration_cpu_rank; //to control the number of multigrid call
	int n_relax_up=100000000;
	double Max_Error_GS=0.00000001;
	int n_relax_down=10;
    
	//===============================================================================================
	ParametersClass(int & i_species, MpiClass * MPI_in ){
		MPI = MPI_in;
		size = MPI->size;
		rank = MPI->rank;
		NODE = MPI->NODE;
		Left_cpu = new int [size];
		Right_cpu = new int [size];
		for (int i=0; i<size; ++i){
			Left_cpu[i]     =  MPI->Left_cpu[i];
			Right_cpu[i]    =  MPI->Right_cpu[i];	
		}
		number_of_data = 5;
		#ifdef print_EPP
		++number_of_data;
		#endif
		#ifdef print_MPP
		++number_of_data;
		#endif
		
			
		Parameters_Initialization_Grid();
		Parameters_Initialization_PhasePoints(i_species);
		
		#ifdef HDF5_Usage
		HDF5_INITIAL();
		#endif
		
		R_Time = 0; 
		
		#ifdef print_Temporal_Techplot
		if (i_species ==0){
			stringstream pathname;
			pathname.str("");
			pathname << "./temporal/Info.dat";
			file_Info.open(pathname.str().c_str());
		}
		#endif
	};
	//===============================================================================================
  
	//*********************************************************************
	void Parameters_Initialization_Grid (){      
		
		multigrid_call_num=0;  
		size_MG=size;
		Iteration_cpu_rank=0;

		
		Read_Grid_info();      
		X_lenght = (  X_max_Tot - X_min_Tot  );
		dX= (  X_max_Tot - X_min_Tot  ) / ( double(this->nX) );
		//       cout<<"++++ Rank="<<rank<<" dX="<<dX<<endl;
		
		if (Num_Spatial_Output <= N_Time) {i_Spatial_Output = N_Time/Num_Spatial_Output;} else{i_Spatial_Output = N_Time;}
		if (Num_Contour <= N_Time ){i_Contour = N_Time/Num_Contour;} else {i_Contour = N_Time;}
		
		
		//Parallel initilization
		nX_cpu = nX/size;
			
		//       X_min= X_min_Tot + (double)rank * ( ( X_max_Tot -X_min_Tot  )/ (double)size);
		//       X_max= X_min + ( ( X_max_Tot - X_min_Tot  )/ (double) size);
		//       iX_min=X_min/dX;
		
		iX_min = rank * nX_cpu;
		X_min  = (double)iX_min * dX; 
		X_max = (double) (iX_min + nX_cpu) * dX;
		velocity_box_simulation = 0.0;
		iX_switching_EPP =int( 310.0/dX);
		nTime_Co_Moving = 1;
		#ifdef Co_moving_frame
		velocity_box_simulation = 6.25;
		nTime_Co_Moving = int(dX/(velocity_box_simulation*d_Time));
		if (rank==0) cout<<" nTime_Co_Moving = "<<nTime_Co_Moving<<endl;
		if (nTime_Co_Moving==0) nTime_Co_Moving=1;
		#endif
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		

			
			
		#ifdef debug
		printf("minimum = %f, max = %f, rank: %i, iX_min:%i \n", X_min, X_max, rank, iX_min);
		#endif

	}
	//*********************************************************************

	//***********************************************************************************************
	void Parameters_Initialization_PhasePoints (int &i_species){
		//Allocation ();
		Read_DF_info(i_species);
		norm_factor = Mass/Temp;
		
		Vx_min = Vx_min*sqrt(1.0/norm_factor) ;
		Vx_max = Vx_max*sqrt(1.0/norm_factor) ;
		dVx = ( abs( Vx_min - Vx_max ) ) / ( double( nVx ) );

		Kappa = 2.0*M_PI /(( abs( X_min_Tot - X_max_Tot) ) / ( double( NoW ) ));
			
		size_shadow_Vx = 7; //default value in case a value for species is not specified in the below
		if (i_species==0) size_shadow_Vx = 5;	
		if (i_species==1) size_shadow_Vx = 9;
		//boundary condition on velocity direction
		VxUp.resize(nVx);      
		for (int iVx=0; iVx<=nVx-1; iVx++){
			VxUp[iVx]=iVx+1;
		}
		//boundary condition
		VxUp[nVx-1]=nVx-1;      
	}
	//***********************************************************************************************
    
	//***********************************************************************************************
	void HDF5_INITIAL(){	
		Dim_X_Vx_Chunk[0]=nX_cpu;
		Dim_X_Vx_Chunk[1]=nVx;
		
		Dim_X_Vx[0]=nX;   
		Dim_X_Vx[1]=nVx;
		
		Dim_X[0]=nX;
		Dim_X_Chunk[0]=nX_cpu;
		
		
		Stride_X[0] = 1 ;
		Stride_X_Vx[0] = 1;
		Stride_X_Vx[1] = 1;
		
		
		Count_X[0] =  1 ;
		Count_X_Vx[0] =  1;
		Count_X_Vx[1] =  1; 
		
		Block_X[0]=Dim_X_Chunk[0];
		Block_X_Vx[0]=Dim_X_Vx_Chunk[0];
		Block_X_Vx[1]=Dim_X_Vx_Chunk[1];
		
		OffSet_X_Vx[0] = iX_min;
		OffSet_X_Vx[1] =0;
		#ifdef debug
			cout<< "{rank="<<rank<<"} offset="<< OffSet_X_Vx[0]<<" ofset+chunk="<< OffSet_X_Vx[1]+Dim_X_Vx_Chunk[1]<<endl;
		#endif
		OffSet_X[0]=OffSet_X_Vx[0]; 
		
		Rank_DFp = 1;
		Rank_X_Vx= 2;
		Rank_X=1; //#define 
		
		
		Rank_3D=3;
		
		Dim_3D[0]=nX;
		Dim_3D[1]=1;
		Dim_3D[2]=1;   
		
		Dim_3D_Chunk[0]=nX_cpu;
		Dim_3D_Chunk[1]=1;
		Dim_3D_Chunk[2]=1;
		
		
		Stride_3D[0] = 1 ;
		Stride_3D[1] = 1 ;    
		Stride_3D[2] = 1 ;
		
		Count_3D[0] =  1 ;
		Count_3D[1] =  1;
		Count_3D[2] =  1; 
		
		Block_3D[0]=Dim_3D_Chunk[0];
		Block_3D[1]=Dim_3D_Chunk[0];
		Block_3D[2]=Dim_3D_Chunk[1];
		
		
		OffSet_3D[0] = iX_min;
		OffSet_3D[1] = 0;
		OffSet_3D[2] = 0;
	}
	//***********************************************************************************************
      

	
	//***********************************************************************************************
	void Read_DF_info(int& i_species) {
		/*take note: arrays are starting from 0, 
		and if the pointers goes beyond the domain of the array, no error will be shown
		the type of the element has to agree with the type of the function */
		size_t flags = 0;
		json_error_t error;
		int i_element;

		if (rank==0) cout << "-----------------------------------------------------"<<endl;
		
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
		
		
		//--- returning the PhasePoint object --------------------------------  
		json_t* PhasePoint_json = json_array_get(config_json, 1);
		if (!PhasePoint_json) cout<<"can't catch 1 member of config=PhasePoint" <<endl;
		PhasePoint_json = json_object_get(PhasePoint_json, "PhasePoint");
		if (!PhasePoint_json) cout<<"can't catch PhasePoint object" <<endl;
		//--------------------------------------------------------------------

		//--- returning the i_species parameters ------------------------------------------
		json_t* Species_json = json_array_get(PhasePoint_json, i_species);
		if (!Species_json) cout<<"can't catch 0 of PhasePointSpecies object" <<endl;
		Species_json = json_object_get(Species_json, "Species");
		if (!Species_json) cout<<"can't catch Species object" <<endl;
		
		
		//---------------------
		i_element=0;
		json_t* element_json      = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element ("<<i_element<<") of " <<i_species <<endl;
		element_json = json_object_get(element_json, "Name");
		if (!element_json) cout<<"can't catch Name of " <<i_species <<endl;
		Name = json_string_value(element_json);
		if (rank==0) cout << "Name of " <<i_species <<" = " << Name << endl;
		
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "Mass");
		if (!element_json) cout<<"can't catch Mass of " <<i_species <<endl;
		Mass = json_real_value(element_json);
		if (rank==0) cout << "Mass of " <<i_species <<" = " << Mass << endl;
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "Charge");
		if (!element_json) cout<<"can't catch Charge of " <<i_species <<endl;
		Charge = json_real_value(element_json);
		if (rank==0) cout << "Charge of " <<i_species <<" = " << Charge << endl;
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element ("<<i_element<<") of " <<i_species <<endl;
		element_json = json_object_get(element_json, "Temperature");
		if (!element_json) cout<<"can't catch Temperature of " <<i_species <<endl;
		Temp = json_real_value(element_json);
		if (rank==0) cout << "Temperature of " <<i_species <<" = " << Temp << endl;
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element ("<<i_element<<") of " <<i_species <<endl;
		element_json = json_object_get(element_json, "Den_Ratio_b");
		if (!element_json) cout<<"can't catch Den_Ratio_b of " <<i_species <<endl;
		Den_Ratio_b = json_real_value(element_json);
		if (rank==0) cout << "Den_Ratio_b of " <<i_species <<" = " << Den_Ratio_b << endl;
		
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "Density");
		if (!element_json) cout<<"can't catch Density of " <<i_species <<endl;
		Density = json_real_value(element_json);
		if (rank==0) cout << "Density of " <<i_species <<" = " << Density << endl;


		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "Alpha");
		if (!element_json) cout<<"can't catch Alpha of " <<i_species <<endl;
		Alpha = json_real_value(element_json);
		if (rank==0) cout << "Alpha of " <<i_species <<" = " << Alpha << endl;
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "NoW");
		if (!element_json) cout<<"can't catch NoW of " <<i_species <<endl;
		NoW = json_integer_value(element_json);
		if (rank==0) cout << "NoW of " <<i_species <<" = " << NoW << endl;    
		
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "pX");
		if (!element_json) cout<<"can't catch pX of " <<i_species <<endl;
		pX = json_integer_value(element_json);
		if (rank==0) cout << "pX of " <<i_species <<" = " << pX << endl;    
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "pVx");
		if (!element_json) cout<<"can't catch pVx of " <<i_species <<endl;
		pVx = json_integer_value(element_json);
		if (rank==0) cout << "pVx of " <<i_species <<" = " << pVx << endl; 
		
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "nVx");
		if (!element_json) cout<<"can't catch nVx of " <<i_species <<endl;
		nVx = json_integer_value(element_json);
		if (rank==0) cout << "nVx of " <<i_species <<" = " << nVx << endl; 
		
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "Vx_min");
		if (!element_json) cout<<"can't catch Vx_min of " <<i_species <<endl;
		Vx_min = json_real_value(element_json);
		if (rank==0) cout << "Vx_min of " <<i_species <<" = " << Vx_min << endl;
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "Vx_max");
		if (!element_json) cout<<"can't catch Vx_max of " <<i_species <<endl;
		Vx_max = json_real_value(element_json);
		if (rank==0) cout << "Vx_max of " <<i_species <<" = " << Vx_max << endl; 
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "DF_X_NMRE");
		if (!element_json) cout<<"can't catch DF_X_NMRE of " <<i_species <<endl;
		DF_X_NMRE = json_integer_value(element_json);
		if (rank==0) cout << "DF_X_NMRE of " <<i_species <<" = " << DF_X_NMRE << endl;
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "DF_Vx_MK");
		if (!element_json) cout<<"can't catch DF_Vx_MK of " <<i_species <<endl;
		DF_Vx_MK = json_integer_value(element_json);
		if (rank==0) cout << "DF_Vx_MK of " <<i_species <<" = " << DF_Vx_MK << endl;
		
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "Kappa_DF");
		if (!element_json) cout<<"can't catch Kappa_DF of " <<i_species <<endl;
		Kappa_DF = json_real_value(element_json);
		if (rank==0) cout << "Kappa_DF of " <<i_species <<" = " << Kappa_DF << endl;
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "Boltzmann_Relation");
		if (!element_json) cout<<"can't catch Boltzmann_Relation of " <<i_species <<endl;
		Boltzmann_Relation = json_integer_value(element_json);
		if (rank==0) cout << "Boltzmann_Relation of " <<i_species <<" = " << Boltzmann_Relation << endl;
		
		
		//---------------------
		++i_element;
		element_json = json_array_get(Species_json, i_element);
		if (!element_json) cout<<"can't catch element  ("<<i_element<<") of "  <<i_species <<endl;
		element_json = json_object_get(element_json, "DF_NotShifted_Shifted");
		if (!element_json) cout<<"can't catch DF_NotShifted_Shifted of " <<i_species <<endl;
		DF_NotShifted_Shifted = json_integer_value(element_json);
		if (rank==0) cout << "DF_NotShifted_Shifted of " <<i_species <<" = " << DF_NotShifted_Shifted << endl;
		
		if (rank==0) cout << "-----------------------------------------------------"<<endl;
	}
	//***********************************************************************************************

  //***********************************************************************************************
  void Read_Grid_info() {
      size_t flags = 0;
      json_error_t error;
      int i_element;

      if (rank==0)    cout << "------- start of grid info---------------------------"<<endl;
    
      //--- loading the JSON file ------------------------------------------
      json_t* file_json = json_load_file("./config.json", flags, &error);
      if(!file_json) {
	cerr << "in line " << error.line << ": " << error.text << endl;
	exit(1);
      }
      //--------------------------------------------------------------------
    
      //--- returning the config object ------------------------------------  
      json_t* config_json = json_object_get(file_json, "config");
      if(!config_json) {
	cerr << "in line " << error.line << ": " << error.text << endl;
	exit(1);
      }
      //--------------------------------------------------------------------
    

      //--- returning the grid object --------------------------------  
      json_t* Grid_json = json_array_get(config_json, 2);
      if (!Grid_json) cout<<"can't catch second member of config= Grid" <<endl;	    
      Grid_json = json_object_get(Grid_json, "Grid");
      if (!Grid_json) cout<<"can't catch Grid object" <<endl;
      //--------------------------------------------------------------------

    
      //--- returning the grid parameters -------------------------------

      //----------------------
      i_element=0;
      json_t* element_json      = json_array_get(Grid_json, i_element);
      if (!element_json) cout<<"can't catch element of 0 of Grid" <<endl;
      element_json = json_object_get(element_json, "nX");
      if (!element_json) cout<<"can't catch nX of Grid" <<endl;
      nX = json_integer_value(element_json);
      if (rank==0) cout << "nX of Grid" <<" = " << nX << endl;
    
      //----------------------
      ++i_element;
      element_json      = json_array_get(Grid_json, i_element);
      if (!element_json) cout<<"can't catch element of 1 of Grid" <<endl;
      element_json = json_object_get(element_json, "X_min");
      if (!element_json) cout<<"can't catch X_min of Grid" <<endl;
      if (!json_is_real(element_json)) cerr <<"That's not an real, rank="<<rank  <<endl;
      X_min_Tot = json_real_value(element_json);
      if (rank==0) cout << "X_min of Grid" <<" = " << X_min_Tot << endl;
    
      //----------------------
      ++i_element;
      element_json      = json_array_get(Grid_json, i_element);
      if (!element_json) cout<<"can't catch element of 2 of Grid" <<endl;
      element_json = json_object_get(element_json, "X_max");
      if (!element_json) cout<<"can't catch X_max of Grid" <<endl;
      X_max_Tot = json_real_value(element_json);
      if (rank==0) cout << "X_max of Grid" <<" = " << X_max_Tot << endl;
   
      if (rank==0) cout << "-------- end of grid info ---------------------------"<<endl;
      
      //--- returning the Time parameters -------------------------------
      if (rank==0) cout << "------- start reading of time info---------------------------"<<endl;
      
      json_t* Time_json = json_array_get(config_json, 3);
      if (!Time_json) cout<<"can't catch 3rd member of config= Time" <<endl;	    
      Time_json = json_object_get(Time_json, "Time");
      if (!Time_json) cout<<"can't catch Time object" <<endl;
      
      //----------------------
      i_element=0;
      element_json      = json_array_get(Time_json, i_element);
      if (!element_json) cout<<"can't catch element of "<<i_element<<" of Time" <<endl;
      element_json = json_object_get(element_json, "N_Time");
      if (!element_json) cout<<"can't catch N_Time of Time" <<endl;
      N_Time = json_integer_value(element_json);
      if (rank==0) cout << "N_Time of Time" <<" = " << N_Time << endl;
    
      //----------------------
      ++i_element;
      element_json      = json_array_get(Time_json, i_element);
      if (!element_json) cout<<"can't catch element of "<<i_element<<" of Time" <<endl;
      element_json = json_object_get(element_json, "d_Time");
      if (!element_json) cout<<"can't catch d_Time of Time" <<endl;
      d_Time = json_real_value(element_json);
      if (rank==0) cout << "d_Time of Time" <<" = " << d_Time << endl; 
      
      
      //----------------------
      ++i_element;
      element_json      = json_array_get(Time_json, i_element);
      if (!element_json) cout<<"can't catch element of "<<i_element<<" of Time" <<endl;
      element_json = json_object_get(element_json, "Num_Spatial_Output");
      if (!element_json) cout<<"can't catch Num_Spatial_Output of Time" <<endl;
      Num_Spatial_Output = json_integer_value(element_json);
      if (rank==0) cout << "Num_Spatial_Output of Time" <<" = " << Num_Spatial_Output << endl; 
      
      //----------------------
      ++i_element;
      element_json      = json_array_get(Time_json, i_element);
      if (!element_json) cout<<"can't catch element of "<<i_element<<" of Time" <<endl;
      element_json = json_object_get(element_json, "Num_Contour");
      if (!element_json) cout<<"can't catch Num_Contour of Time" <<endl;
      Num_Contour = json_integer_value(element_json);
      if (rank==0) cout << "Num_Contour of Time" <<" = " << Num_Contour << endl; 
      
      //----------------------
      ++i_element;
      element_json      = json_array_get(Time_json, i_element);
      if (!element_json) cout<<"can't catch element of "<<i_element<<" of Time" <<endl;
      element_json = json_object_get(element_json, "Plasma_approximation");
      if (!element_json) cout<<"can't catch Plasma_approximation of Time" <<endl;
      Plasma_approximation = json_integer_value(element_json);
      if (rank==0) cout << "Plasma_approximation " <<" = " << Plasma_approximation << endl; 
      
      //----------------------
      ++i_element;
      element_json      = json_array_get(Time_json, i_element);
      if (!element_json) cout<<"can't catch element of "<<i_element<<" of Time" <<endl;
      element_json = json_object_get(element_json, "convergance_step");
      if (!element_json) cout<<"can't catch convergance_step of Time" <<endl;
      convergance_step = json_integer_value(element_json);
      if (rank==0) cout << "convergance_step " <<" = " << convergance_step << endl; 
      if (rank==0) cout << "-------- end of Time info ---------------------------"<<endl;         
  } 
  //***********************************************************************************************
  
  
  //***********************************************************************************************
  void Ghost_Grid_Points_Vector(vector<double>& array, int & nX_array, int jump_num) {
        
 //uses balck-red approach for send and recieve
    int rank_snd, rank_rcv;
    vector<double> * array_RCV, * array_SND;
    array_RCV  = new vector<double>(2);
    array_SND  = new vector<double>(2);
    
    
    array_SND->at(0) = array[nX_array-1];
    array_SND->at(1) = array[nX_array-2];
    
    
    //---- to the right -----------------------------------------------------------------------
    if (rank%(jump_num*2)==0) {
	rank_snd = rank + jump_num;
	if (rank_snd > NODE) rank_snd=0; //boundary condition over cpus
	MPI_Send(&array_SND->at(0), 2, MPI_DOUBLE,rank_snd,233, MPI_COMM_WORLD);
    }else{
	rank_rcv= rank - jump_num;
	if (rank_rcv<0) rank_rcv= size-jump_num; //boundary condition over cpus
	MPI_Recv(&array_RCV->at(0), 2, MPI_DOUBLE,rank_rcv,233, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    if (rank%(jump_num*2)==0) {    
	rank_rcv= rank - jump_num;
	if (rank_rcv<0) rank_rcv= size-jump_num; //boundary condition over cpus
	MPI_Recv(&array_RCV->at(0), 2, MPI_DOUBLE,rank_rcv,233, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   
    }else{
	rank_snd = rank + jump_num;
	if (rank_snd > NODE) rank_snd=0; //boundary condition over cpus
	MPI_Send(&array_SND->at(0), 2, MPI_DOUBLE,rank_snd,233, MPI_COMM_WORLD);
    }


    
    minus_1 = array_RCV->at(0);
    minus_2 = array_RCV->at(1);
    
    
    array_SND->at(0) = array[0];
    array_SND->at(1) = array[1];
    
    
    // ---- to the left ------------------------------------------------------
    if (rank%(jump_num*2)==0) {
	rank_snd = rank - jump_num;
	if (rank_snd<0) rank_snd= size-jump_num; //boundary condition over cpus
	MPI_Send(&array_SND->at(0), 2, MPI_DOUBLE,rank_snd,377, MPI_COMM_WORLD);
    } else {
	rank_rcv= rank + jump_num;
	if (rank_rcv > NODE) rank_rcv=0; //boundary condition over cpus
	MPI_Recv(&array_RCV->at(0), 2, MPI_DOUBLE,rank_rcv,377, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    
    if (rank%(jump_num*2)==0) {
        rank_rcv= rank + jump_num;
	if (rank_rcv > NODE) rank_rcv=0; //boundary condition over cpus
	MPI_Recv(&array_RCV->at(0), 2, MPI_DOUBLE,rank_rcv,377, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
	rank_snd = rank - jump_num;
	if (rank_snd<0) rank_snd= size-jump_num; //boundary condition over cpus
	MPI_Send(&array_SND->at(0), 2, MPI_DOUBLE,rank_snd,377, MPI_COMM_WORLD);
    }    
    
    
    ghost_nX_cpu  = array_RCV->at(0);
    nX_cpu_plus_1 = array_RCV->at(1);
    
    delete array_RCV;
    delete array_SND;
    
    
    }
  //***********************************************************************************************
    
     
  //***********************************************************************************************
  void Sum_up_Double(double & var, double& var_tot  ,int jump_num) {
        
    int rank_rcv, rank_snd;
    double var_rcv;
    if (rank%jump_num==0) {
	if (rank%(jump_num*2)==0) {
	  rank_rcv = rank + jump_num;  
// 	  printf("%i<--- %i,  level = %i \n",rank,rank_rcv,jump_num);
	  MPI_Recv(&var_rcv, 1, MPI_DOUBLE,rank_rcv,233, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  var = var + var_rcv;
	} else {
	  rank_snd = rank - jump_num;  
	  MPI_Send(&var, 1, MPI_DOUBLE,rank_snd,233, MPI_COMM_WORLD);
	  
	}
	jump_num = jump_num * 2;
	if(jump_num == size){
	  if (rank==0) var_tot = var;
	}else{
	  Sum_up_Double(var, var_tot, jump_num);
	}
    }
  }
  //***********************************************************************************************
  
  //***********************************************************************************************
  void Sum_up_INT(int & var, int& var_tot  ,int jump_num) {
        
    int rank_rcv, rank_snd;
    int var_rcv;
    if (rank%jump_num==0) {
	if (rank%(jump_num*2)==0) {
	  rank_rcv = rank + jump_num;  
// 	  printf("%i<--- %i,  level = %i \n",rank,rank_rcv,jump_num);
	  MPI_Recv(&var_rcv, 1, MPI_DOUBLE,rank_rcv,233, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  var = var + var_rcv;
	} else {
	  rank_snd = rank - jump_num;  
	  MPI_Send(&var, 1, MPI_DOUBLE,rank_snd,233, MPI_COMM_WORLD);
	  
	}
	jump_num = jump_num * 2;
	if(jump_num == size){
	  if (rank==0) var_tot = var;
	}else{
	  Sum_up_INT(var, var_tot, jump_num);
	}
    }
  }
  //***********************************************************************************************
  
  //***********************************************************************************************
  void Gather_Compare_Cast_INT(int & var, int& var_tot  ,int jump_num, const char*& Debug_mes) {
    //Recrusive Tree Based function to compare var and get the var_tot and cast it 
    //to all the cpus existing inside the comm world with rank=i*jump_num, i=0,1,2, ... jump_num here refers to the initial value of jump_num
    int rank_rcv, rank_snd;
    int var_rcv;
    if (rank%jump_num==0) {
	if (rank%(jump_num*2)==0) {
	  rank_rcv = rank + jump_num;  
	  //printf("%i<--- %i,  level = %i \n",rank,rank_rcv,jump_num);
	  MPI_Recv(&var_rcv, 1, MPI_INT,rank_rcv,233, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  var = var + var_rcv;
	  if (var>=1) var=1; //important part which make sure that the comparison always ends up returning 1
	} else {
	  rank_snd = rank - jump_num;  
	  MPI_Send(&var, 1, MPI_INT,rank_snd,233, MPI_COMM_WORLD);
	  
	}
	jump_num = jump_num * 2;
	
	if(jump_num == size){
	  if (rank==0) var_tot = var;
	}else{
	  Gather_Compare_Cast_INT(var, var_tot, jump_num, Debug_mes);
	  
	}
	jump_num = jump_num / 2;
	
	if (rank%(jump_num*2)==0) {
	  var = var_tot;
	  rank_snd = rank + jump_num;  
//   	  printf("%i---> %i,  level = %i \n",rank,rank_snd,jump_num);
	  MPI_Send(&var, 1, MPI_INT,rank_snd,377, MPI_COMM_WORLD);  
	} else {
	  rank_rcv = rank - jump_num;  
	  MPI_Recv(&var_rcv, 1, MPI_INT,rank_rcv,377, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  var_tot = var_rcv;
	}
//  	cout <<"finished, level"<<jump_num<<"  rank="<<rank<<"  "<<R_Time<<"   "<<Debug_mes <<endl;
    }
  }
  //***********************************************************************************************
  
  
    
  //***********************************************************************************************
  void Print_Debug(int& Line_NO, const char*& Debug_mes) {
             
      MPI_Barrier(MPI_COMM_WORLD);
	if (rank==0){
	  cout<<endl;
	  cout<< "------------------------------------------------------"<<endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	#ifdef debug_rank_All
	    int flag_debug=1;
	    if (rank==0) {
	      cout<<Debug_mes <<" {"<<rank<<"} "<<Line_NO <<endl;
	      for (int i=1; i<=NODE; i++){
		MPI_Send(&flag_debug, 1, MPI_INT,i,2, MPI_COMM_WORLD);
		MPI_Recv(&flag_debug, 1, MPI_INT,i,3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      }
	    }else{
	      MPI_Recv(&flag_debug, 1, MPI_INT,0,2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
	      cout<<Debug_mes <<" {"<<rank<<"} "<<Line_NO <<endl;
	      MPI_Send(&flag_debug, 1, MPI_INT,0,3, MPI_COMM_WORLD);
	    }
	#else
	  if (rank==0) cout<<Debug_mes <<" {"<<rank<<"} "<<Line_NO <<endl;
	#endif
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank==0){
	  cout<< "------------------------------------------------------"<<endl;
	  cout<<endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
    }
  //***********************************************************************************************
    
  //***********************************************************************************************
  void Print_Debug_Single(int& Line_NO, const char*& Debug_mes) {
      cout<<endl;
      cout<< "==========================================="<<endl;
      cout<<Debug_mes <<" {"<<rank<<"} "<<Line_NO <<endl;
    }
  //***********************************************************************************************
      
  };
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#endif //ParametersClass_Included
