#ifndef GridClass_Included
#define GridClass_Included
#include "config.h"


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
class GridClass {

public:
  
  ParametersClass * P;
	//---- MPI variables ------
  	int size, rank,NODE;
	int * Left_cpu;
	int * Right_cpu;
  
  ofstream file_RO,file_Ex,file_PHI,file_Time_Ex,file_Time_Ene_ele;
  
  int jump_num; //a variable checking the number of cpus living in each level of Tree-Based-Multigrid

  // electric energy variables
  double Ene_ele_Tot, Ene_ele_cpu;  

  // vectors 
   vector<double> RO_X,RO_X_Old, PHI_X;

  vector<double> E_X;
  vector<double> E_X_0;
  vector<double> E_X_ghost;
  
#if flag_Poisson_Solver==0
	vector<double> Matrix_XX_cpu;
  
#endif
  
  //===============================================================================================
  GridClass (ParametersClass * P_in, MpiClass * MPI){
    P = P_in;
    
    size = MPI->size;
    rank = MPI->rank;
    NODE = MPI->NODE;
    Left_cpu = new int [size];
    Right_cpu = new int [size];
    for (int i=0; i<size; ++i){
	Left_cpu[i]     =  MPI->Left_cpu[i];
	Right_cpu[i]    =  MPI->Right_cpu[i];	
    }


 
  };
  //===============================================================================================

  
  
	//***********************************************************************************************
	void Read_Matrix_from_json_file() {
		/*take note: arrays are starting from 0, 
		* and if the pointers goes beyond the domain of the array, no error will be shown
		* the type of the element has to agree with the type of the function */
		size_t flags = 0;
		json_error_t error;


		if (rank==0) cout << "-----------------------------------------------------"<<endl;
		
		//--- loading the JSON file ------------------------------------------
		json_t* file_json = json_load_file("./Inversed_Matrix.json", flags, &error);
		if(!file_json){
			cerr << "in line " << error.line << ": " << error.text << endl;
			exit(1);
		}
		
		//--- returning the phi object ------------------------------------  
		json_t* matrix_element_json = json_object_get(file_json, "Matrix_elements");
		if(!matrix_element_json){
			cerr << "in line " << error.line << ": " << error.text << endl;
			exit(1);
		}
		//--------------------------------------------------------------------
		
		
		//------reading the matrix--------------------------------------------
		int iRow_cpu = 0;
		for (int iRow = rank*P->nX_cpu; iRow<(rank + 1)*P->nX_cpu; ++iRow) {
			json_t* Row_json = json_array_get(matrix_element_json, iRow);
			if (!Row_json) {
				cout<<"can't catch element = ("<<iRow<<") Rank="<<rank <<endl;
			}else{
				
				for (int iColumn=0; iColumn<P->nX; ++iColumn) {
					json_t* element_json = json_array_get(Row_json, iColumn);
					if (!element_json) {
						cout<<"can't catch element IRow= ("<<iRow<<")"<< " iColumn=("<<iColumn<<") Rank="<<rank <<endl;
					}else{
						Matrix_XX_cpu[iColumn + (iRow_cpu*P->nX)] = json_real_value(element_json);
					}
				}
				++iRow_cpu;
				
			}
		}
	}
      
      
      
  
  //***********************************************************************************************
  void Initialization_GridClass (){
    #ifdef debug_initial
      P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="--- Initial Step (0): initialization of grid class");
    #endif
      
      
    E_X.resize(P->nX_cpu);
    E_X_0.resize(P->nX_cpu); // storing for leapfrog trapezoidal method
    E_X_ghost.resize(P->nX_cpu+3); // three points are added -1, nX_cpu, nX_cpu+1

    RO_X.resize (P->nX_cpu);
    RO_X_Old.resize (P->nX_cpu);
    PHI_X.resize (P->nX_cpu);
    fill(PHI_X.begin(), PHI_X.end(), 0.0);
    fill(E_X_0.begin(), E_X_0.end(), 0.0);
	#if flag_Poisson_Solver==0
		Matrix_XX_cpu.resize (P->nX * P->nX_cpu);
		Read_Matrix_from_json_file();
	#endif
    stringstream pathname;
    
      
    #ifdef print_Temporal_Techplot
      MPI_Barrier(MPI_COMM_WORLD);
      if (size>=4){
        if (rank%(size/4)==0) {
            //int iX_half=int((P->nX_cpu/2));
            pathname.str("");
            pathname << "./temporal/Ex_T_"<<rank<<".dat";
            file_Time_Ex.open(pathname.str().c_str());

            if ( (file_Time_Ex.rdstate() & std::ifstream::failbit ) != 0 )
                cout << "Error opening Temporal E_X for rank="<<rank<<endl;
        }
      }else{
            pathname.str("");
            pathname << "./temporal/Ex_T_"<<rank<<".dat";
            file_Time_Ex.open(pathname.str().c_str());

            if ( (file_Time_Ex.rdstate() & std::ifstream::failbit ) != 0 )
                cout << "Error opening Temporal E_X for rank="<<rank<<endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      
    #endif 
      
    #ifdef print_Temporal_Techplot  
      if (rank==0) {
	pathname.str("");
	pathname << "./temporal/Ene_ele.dat";
	file_Time_Ene_ele.open(pathname.str().c_str());
      }
    #endif 
      
    }
  //***********************************************************************************************
  
  //***********************************************************************************************
  void RO_Zero (){
      fill(RO_X.begin(), RO_X.end(), 0.0);
    }
  //***********************************************************************************************
  
  //***********************************************************************************************
  void RO_Zero_on_borders (){ // open boundary condition
		if (rank == 0){
			for (int iX = 0; iX<P->nX_cpu; ++iX) {
				RO_X[iX] = 0.0;
			}
	
		}
		if  (rank == NODE){
			for (int iX = 0; iX<P->nX_cpu; ++iX) {
				RO_X[iX] = 0.0;
			}
		}	

    }
    
  //***********************************************************************************************
  
  
  //***********************************************************************************************
  void Ex_Zero (){
      fill(E_X.begin(), E_X.end(), 0.0);
      fill(E_X_ghost.begin(), E_X_ghost.end(), 0.0);
    }
  //***********************************************************************************************
  //***********************************************************************************************
  void RO_Calculator (vector<double>& Den_g){
      for (int iX = 0; iX<=P->nX_cpu-1; iX++) {
	RO_X[iX] +=  Den_g[iX]; 
    }
  }
  //***********************************************************************************************
  
    //***********************************************************************************************
  void RO_Old_Store(){
      for (int iX = 0; iX<=P->nX_cpu-1; iX++) {
	RO_X_Old [iX] = RO_X[iX];
    }
  }
  //***********************************************************************************************

  //***********************************************************************************************
  void Quasineutrality (double & Charge){ 
      for (int iX = 0; iX<=P->nX_cpu-1; iX++) {
	RO_X[iX] += -Charge ; // to follow quasineutrality when there is just one species
      }
    }
  //***********************************************************************************************
 
 
  //***********************************************************************************************
  void Compare_RO (double & RO_diff){ 
    RO_diff = 100.0  ;
    for (int iX = 0; iX<=P->nX_cpu-1; iX++) {  
      if (abs(RO_X[iX]-RO_X_Old[iX])<RO_diff) {
	RO_diff = abs(RO_X[iX]-RO_X_Old[iX]);
      }
     }
   }
  //***********************************************************************************************
  
  //***********************************************************************************************
  void PHI_Static (int & i_converage){ 
    double X_static;
    double Lenght = abs( P->X_min_Tot - P->X_max_Tot);
    double veclocity = 1.6;
    double X_minimum = Lenght/double(P->NoW)*(5.0);//, X_maximum = Lenght/double(P->NoW)*(3.75);
    for (int iX = 0; iX<=P->nX_cpu-1; iX++)  { //to creat soliton in the code
	X_static = P->X_min  + double(iX)  * P->dX;
	PHI_X[iX] = P->Alpha * cos(P->Kappa * (X_static- veclocity*i_converage*P->d_Time)) * exp(-pow(abs(X_static-X_minimum)/5.0,2));
//  	if (X_static>X_minimum && X_static<X_maximum) PHI_X[iX]  = P->Alpha * cos(P->Kappa * (X_static- veclocity*i_converage*P->d_Time)  );
// 	PHI_X[iX]  = P->Alpha * exp( -pow ( ((X_static- veclocity*i_converage*P->d_Time)  - Lenght/2.0)*10.0 ,2) );
	
    }
  }
	//*********************************************************************************************** 
	void test_poisson_solver(){
		double X0;
		for (int iX = 0; iX<=P->nX_cpu-1; iX++)  { 
			X0 = P->X_min  + double(iX)  * P->dX;
			RO_X[iX]  = -0.1;//P->Alpha * cos(P->Kappa * (X0)  );
		}
	}
	//*********************************************************************************************** 
  
  

  
	//***********************************************************************************************
	void Poisson_Solver(){
		fill(PHI_X.begin(), PHI_X.end(), 0.0);
//  		test_poisson_solver();
		#ifdef debug_Poisson
			P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="------  Start of Poisson Solver ---");
		#endif
			
		#ifdef flag_boundary_condition_open
			RO_Zero_on_borders();
		#endif
			
// 		Zalesak_limiter_fourth_order(RO_X);
		if (P->Plasma_approximation == 0) { 
		Minus_RO_for_Poisson (RO_X,P->nX_cpu);
		#if flag_Poisson_Solver==0
			Poisson_solver_inverse_matrix(PHI_X,RO_X,Matrix_XX_cpu);
		#else
			if (size ==1) {
				Multigrid_Solver_Single(PHI_X, RO_X, P->nX, P->dX);
			}else{
			//  	  Multigrid_Solver_Parallel(PHI_X, RO_X, P->nX_cpu, P->dX);
				int jump_num = 1;
				Poisson_Solver_GS_Parallel(PHI_X,RO_X, P->nX_cpu, P->dX ,P->n_relax_up,P->Max_Error_GS,jump_num);
			}
		#endif
// 		Zalesak_limiter_fourth_order(PHI_X);
		}else if (P->Plasma_approximation == 1){
			PHI_calculator(PHI_X, RO_X, P->nX_cpu);
		}
 		E_x_Calculator();
// 		E_x_Calculator_2points();
		Ene_ele_Calculator();
	}
	//***********************************************************************************************
	
	//***********************************************************************************************
	void Poisson_solver_inverse_matrix (vector<double>& PHI_X, vector<double>& RO_X, vector<double>& Matrix_XX_cpu){
		vector<double> RO_X_total(P->nX);
		double multiplication_value;
		double average_value, X0;
		
		// gather and broadcast Ro, so each cpu has a copy of all elements of RO
		MPI_Gather(&RO_X.at(0),       	P->nX_cpu, MPI_DOUBLE,
		           &RO_X_total.at(0), 	P->nX_cpu, MPI_DOUBLE,
					   0, MPI_COMM_WORLD);
		#ifdef debug_Poisson
			P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="------  MPI_Gather for Ro_X is finsished ---");
		#endif
		/*	
		if (rank==0){
			// restrictive condition on the matrix, boundaries = 0.0
// 	 		RO_X_total.at(0) = 0.0;
			// Correction to Ro: adjusting average to zero	
			for (int i =0; i<3; ++i){
				average_value = 0.0;
				int i_counter =0;
				for (int iX = 0; iX<25; ++iX) {
					average_value = average_value + RO_X_total[iX];
					++i_counter;
				}
				for (int iX = P->nX-25; iX<P->nX; ++iX) {
					average_value = average_value + RO_X_total[iX];
					++i_counter;
				}
				average_value = average_value/ float(i_counter);//;/ float(P->nX);
				for (int iX = 0; iX<P->nX; ++iX) {
					RO_X_total[iX] = RO_X_total[iX] - average_value;
				}
	// 			cout<<"+++ average value="<<average_value<<std::flush; ;
			}
// 			average_value = 0.0;//RO_X_total[0];
// 			for (int iX = 0; iX<P->nX; ++iX) {
// 				if (abs(RO_X_total[iX])>average_value){
// 					average_value=RO_X_total[iX];
// 				}
// 			for (int iX = 0; iX<P->nX; ++iX) {
// 				X0 = iX*P->dX +P->X_min;
// 				RO_X_total[iX] = average_value*cos(P->Kappa*X0);
// 			}
// 			}
					
			
		}
		*/
		
			
		MPI_Bcast(&RO_X_total.at(0), 	P->nX, MPI_DOUBLE,
		         			0, MPI_COMM_WORLD);
		#ifdef debug_Poisson
			P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="------  MPI_Bcast for Ro_X is finsished ---");
		#endif

			
		// multiplication of inverse matrix with Ro_total to get the PHI_X on each cpu
		for (int iRow_cpu = 0; iRow_cpu<P->nX_cpu; ++iRow_cpu) {
			multiplication_value = 0.0;
			for (int iColumn=0; iColumn<P->nX; ++iColumn) {
					multiplication_value +=	Matrix_XX_cpu[iColumn + (iRow_cpu*P->nX)] * RO_X_total[iColumn] ;
			}
			PHI_X[iRow_cpu]=multiplication_value;
		}
		#ifdef debug_Poisson
			P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="------  multiplication of inverse matrix is finsished ---");
		#endif
	}
	//***********************************************************************************************
    
    
  //***********************************************************************************************
  void Minus_RO_for_Poisson (vector<double>& RO_X, int & nX_cpu)  {
      for (int iX = 0; iX<=nX_cpu-1; iX++) {  
	RO_X[iX] = - RO_X[iX];
      }   
  }
    
    
  //***********************************************************************************************
	//***********************************************************************************************
	void Save_E_in_E0 ()  {
		for (int iX = 0; iX<P->nX_cpu; ++iX) {  
			E_X_0[iX] = E_X[iX];
		}
	}
	
	//***********************************************************************************************
	void Average_E_E0_store_in_E ()  {
		/* average E and E0 into ---> E, and return E0= old E */
		double E_X_aux;
		for (int iX = 0; iX<P->nX_cpu; ++iX) {  
			E_X_aux 	= E_X[iX];
			E_X[iX] 	= (E_X_0[iX] + E_X[iX])/2.0;
			E_X_0[iX] 	= E_X_aux;
		}
	}	
	//***********************************************************************************************
  
  //***********************************************************************************************
  void PHI_calculator (vector<double>& PHI_X, vector<double>& RO_X, int & nX_cpu)  {
      for (int iX = 0; iX<=nX_cpu-1; iX++) {	
	PHI_X[iX] = - (P->Temp/P->Charge) *  log(RO_X[iX]); //CAUTION about usage of P-> here which refers just to the first element (must be electrons)
							    // minus in the eqation comes from the fact that \E_x = - \deravative \phi
							    // apparently ln(RO_X) can be repalced by RO_X
      }   
  }
    
    
  //***********************************************************************************************
   
  //***********************************************************************************************
  void Multigrid_Solver_Parallel (vector<double>& Error_S_ini, vector<double>& Resid_S_ini,
				int & nX_S_ini, double & dX_S_ini )  {
// 	int P->size_MG = size; // this is being taking care of inside this function, and at the end of each time step, size_MG = Size and Iteration_cpu_rank=0
    //WARNING print_debug should not be used here since on the second level tree onward there, all the cpus are not present and it produce deadlock!!
    // #ifdef debug_Poisson
    //P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="------ (0): Start of Multigrid_Solver_Parallel ");
    //#endif 
    

    P->multigrid_call_num +=1;
    

    
//     int n_relax_down=10, n_relax_up=1000000; 
    int Direct_Grid_Number=4;
    int nX_S_out=nX_S_ini/2;
    double dX_T_s_1 = 2.0*dX_S_ini;
//     double Max_Error_GS=0.000001;
    
    jump_num = pow(2,P->Iteration_cpu_rank);
    
    #ifdef debug_Poisson
	if (rank == 0){
	  for (int i = 0; i< (size/jump_num); i++){
	    cout <<"|"<<nX_S_ini;
	    for (int j = 0; j<= jump_num; j++) {
	      cout <<" ";
	    }
	  }
	  cout<<endl;	    
	}
    #endif
    
    vector<double> Resid_S_mid(nX_S_ini),Error_S_cor(nX_S_ini);
    vector<double> Error_S_out(nX_S_out),Resid_S_out(nX_S_out);


     fill(Resid_S_out.begin(), Resid_S_out.end(), 0.0);
     fill(Resid_S_mid.begin(), Resid_S_mid.end(), 0.0);

       
       //--- (0) pre relaxation -----
       //cout <<"------ pre relaxtion rank="<<rank<<"   number of multigrid calls="<<P->multigrid_call_num<<endl;
      Poisson_Solver_GS_Parallel(Error_S_ini,Resid_S_ini,nX_S_ini,dX_S_ini,P->n_relax_down,P->Max_Error_GS,jump_num);
 
       //--- (1) computing the residual ----
       //cout <<"------ (1) computing the residual rank="<<rank<<"   number of multigrid calls="<<P->multigrid_call_num<<endl;
      Residual_Computation_Parallel (Error_S_ini,Resid_S_ini,nX_S_ini,dX_S_ini,Resid_S_mid);
  
       //--- (2) restricting the residual ---- 
       //cout <<"------ (2) restricting the residual rank="<<rank<<"   number of multigrid calls="<<P->multigrid_call_num<<endl;
      Restrict_RO_Parallel (Resid_S_mid,nX_S_ini,Resid_S_out,nX_S_out);
 
       
       //--- (3) guess zero for the next level initial value ------
//       fill(Error_S_out.begin(), Error_S_out.end(), 0.0);

      //==== (4) next iteration ===============================================
      if(nX_S_ini > Direct_Grid_Number) {
	  Multigrid_Solver_Parallel(Error_S_out, Resid_S_out, nX_S_out, dX_T_s_1);
      } else {
	  vector<double> Error_S_out_array(Direct_Grid_Number+1);
	  vector<double> Resid_S_out_array(Direct_Grid_Number+1);
	  
	  P->Iteration_cpu_rank +=1;
	  P->size_MG = P->size_MG/2;
	  int num_level=(int) pow(2,P->Iteration_cpu_rank); //2, 4, 8
	  int num_level_minus_1= num_level-(int) pow(2,P->Iteration_cpu_rank-1); //1, 2, 4

	
 	  if (rank% num_level==0 ){
 	      int rank_rcv = rank+num_level_minus_1;
 	      MPI_Recv(&Error_S_out_array[nX_S_out], nX_S_out, MPI_DOUBLE,rank_rcv,233, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 	      MPI_Recv(&Resid_S_out_array[nX_S_out], nX_S_out, MPI_DOUBLE,rank_rcv,377, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  	      for(int i=0;i<=nX_S_out-1;i++) Error_S_out_array[i]=Error_S_out[i];
 	      for(int i=0;i<=nX_S_out-1;i++) Resid_S_out_array[i]=Resid_S_out[i];
 	  }else if (rank% num_level_minus_1==0 ) {
 		int rank_snd= rank-num_level_minus_1;
		MPI_Send(&Error_S_out[0], nX_S_out, MPI_DOUBLE,rank_snd,233, MPI_COMM_WORLD);
 		MPI_Send(&Resid_S_out[0], nX_S_out, MPI_DOUBLE,rank_snd,377, MPI_COMM_WORLD);
	  }
 	 
	 if (P->size_MG == 1) {
	    if (rank==0) Multigrid_Solver_Single(Error_S_out_array,Resid_S_out_array,Direct_Grid_Number,dX_T_s_1);
	 }else {
	    if (rank% num_level==0 ) Multigrid_Solver_Parallel(Error_S_out_array, Resid_S_out_array, Direct_Grid_Number, dX_T_s_1);	      
	 }
 	/*if (P->multigrid_call_num==5) {
 	  cout <<"Start (MG Parallel) number of multigrid calls="<<P->multigrid_call_num
   	<<"  level of cpu reduction="<<P->Iteration_cpu_rank
   	<<"   num_level="<<num_level
   	<<"  rank ="<<rank<<endl;
 	}*/
	      
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
 	if (rank% num_level==0 ){
 		int rank_snd = rank+num_level_minus_1;
 		MPI_Send(&Error_S_out_array[nX_S_out] , nX_S_out, MPI_DOUBLE,rank_snd,610, MPI_COMM_WORLD);
 		MPI_Send(&Resid_S_out_array[nX_S_out], nX_S_out, MPI_DOUBLE,rank_snd,987, MPI_COMM_WORLD);
 		for(int i=0;i<=nX_S_out-1;i++) Error_S_out[i]=Error_S_out_array[i];
 		for(int i=0;i<=nX_S_out-1;i++) Resid_S_out[i]=Resid_S_out_array[i];
 	  }else if (rank% num_level_minus_1==0 ) {
 		int rank_rcv= rank-num_level_minus_1;
		MPI_Recv(&Error_S_out[0], nX_S_out, MPI_DOUBLE,rank_rcv,610, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&Resid_S_out[0], nX_S_out, MPI_DOUBLE,rank_rcv,987, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 	  }
	 P->size_MG = P->size_MG*2;
	 P->Iteration_cpu_rank -=1;
	 jump_num = pow(2,P->Iteration_cpu_rank);
	 ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      }
      //=======================================================================
      
      
      
	//--- (5) prolon the Error or PHI -----------------
	fill(Error_S_cor.begin(), Error_S_cor.end(), 0.0);
	Prolong_Parallel (Error_S_out,nX_S_out,Error_S_cor,nX_S_ini);
    
	//--- (6) adding the correction to the Error or PHI --------
	for (int iX = 0; iX<=nX_S_ini-1; iX++) {
	  Error_S_ini[iX]=-Error_S_ini[iX] + Error_S_cor[iX];
	}
	//Multigrid_PHI_Printing_Single(Error_S_ini,nX_S_ini,dX_S_ini,Error_Tot,nX_S_ini)
  
	//--- (7) post relaxation -------------------------------------------------------
	Poisson_Solver_GS_Parallel(Error_S_ini,Resid_S_ini,nX_S_ini,dX_S_ini,P->n_relax_up,P->Max_Error_GS,jump_num);
	P->multigrid_call_num -=1;
	
	#ifdef debug_Poisson
	if (rank == 0){
	  for (int i = 0; i< (size/jump_num); i++){
	    cout <<"|"<<nX_S_ini;
	    for (int j = 0; j<= jump_num; j++) {
	      cout <<" ";
	    }
	  }
	  cout<<endl;	    
	}
	#endif
    

  } 
  //-----------------------------------------------------------------------------------------------

  //-----------------------------------------------------------------------------------------------
  void Poisson_Solver_GS_Parallel (vector<double>& PHI_X_ini,vector<double>& RO_X_ini
			      ,int& nX_GS,double& dX_GS,int& Max_Loop,double& Max_Error_GS, int& jump_num){
  
      #ifdef debug_Poisson
	if (rank == 0)cout <<" ((rank ="<<rank<<", start of Poisson_Solver_GS_Parallel))  "<<std::flush; 
      #endif
      double Resid_GS_2,Resid_GS_1;
      int Flag_GS_Single,Flag_GS_Total;
      vector<int> Flag_GS_cpu_array(size);
      int i_GS;
      vector<double> Delta_PHI(nX_GS);

      
//       double ghost_minus_2,ghost_minus_1,ghost_nX_cpu,ghost_nX_cpu_plus_1;
   
      fill(Delta_PHI.begin(), Delta_PHI.end(), 0.0);
      Resid_GS_1=0.;
      Resid_GS_2=0.;
   
      i_GS=0;
   
   
     do {

       Flag_GS_Single=1;
       Flag_GS_Total=1;
       
       //-------------------------- Jacobi method of solving Poisson equation ----------------
	  i_GS += 1;
	  P->Ghost_Grid_Points_Vector(PHI_X_ini,nX_GS,jump_num);
	  
 	  PHI_X_ini[0] = P->minus_1 + PHI_X_ini[0+1];//CPU Boundary Correction
 	  PHI_X_ini[0] = PHI_X_ini[0] - (RO_X_ini[0]*( pow(dX_GS,2)));
 	  PHI_X_ini[0] = 0.5 * PHI_X_ini[0];
	      
	  PHI_X_ini[nX_GS-1] = PHI_X_ini[nX_GS-1-1]+ P->ghost_nX_cpu;//CPU Boundary Correction
	  PHI_X_ini[nX_GS-1] = PHI_X_ini[nX_GS-1] - (RO_X_ini[nX_GS-1]*( pow(dX_GS,2)));
	  PHI_X_ini[nX_GS-1] = 0.5 * PHI_X_ini[nX_GS-1];
	  
	  for (int iX = 1; iX<=nX_GS-2; ++iX) {
	      PHI_X_ini[iX] = PHI_X_ini[iX-1]+ PHI_X_ini[iX+1];
	      PHI_X_ini[iX] = PHI_X_ini[iX] - (RO_X_ini[iX]*( pow(dX_GS,2)));
	      PHI_X_ini[iX] = 0.5 * PHI_X_ini[iX];
	  }
	  
	  
	  //-------------------------------------------------------------------------------------
	  

	  //----- calculation of Delta_PHI ------------------------------------------------------
	  P->Ghost_Grid_Points_Vector(PHI_X_ini,nX_GS,jump_num);
	  
	  for (int iX = 1; iX<=nX_GS-2; iX++) {
	      Delta_PHI[iX] = PHI_X_ini[iX-1]-2.0*PHI_X_ini[iX]+PHI_X_ini[iX+1];
	  }
	  Delta_PHI[0] =  P->minus_1 - 2.0*PHI_X_ini[0] + PHI_X_ini[1];
	  Delta_PHI[nX_GS-1] = PHI_X_ini[nX_GS-2]-2.0*PHI_X_ini[nX_GS-1]+P->ghost_nX_cpu;
	  
	  for (int iX = 0; iX<=nX_GS-1; iX++) {
	      Delta_PHI[iX] = Delta_PHI[iX]/( pow(dX_GS,2));      
	  }
	  //-------------------------------------------------------------------------------------
     
	  //---- Finding the maximum residual -------------
	  Resid_GS_2 = 0.0;
	  for (int iX = 0; iX<=nX_GS-1; iX++)  {
	      if ( abs( RO_X_ini[iX] + Delta_PHI[iX] )> Resid_GS_2 )
	      Resid_GS_2  = abs( RO_X_ini[iX] + Delta_PHI[iX] );
	  }
	  //--------------------------------------------
      
	  //---- comparing the maximum error with the acceptable error -----------------------------
	  if ( i_GS > 2) { // this if makes sure that the initial stage is not going to be checked which resutl in division by zero
	      if ( abs((Resid_GS_2-Resid_GS_1)/Resid_GS_1) < Max_Error_GS) {
		  Flag_GS_Single = 0;
		  //Print*,"converged",i_GS,rank,Time
	      } else if(i_GS>Max_Loop) {
		  Flag_GS_Single = 0;
		  //  Print*,"Not converging",rank,Time
	      }	  	
	  }
	  Resid_GS_1 = Resid_GS_2;
	  
	  //-------------------------------- set the flag for each cpu
	  // gather, sum up, broadcast the flags
//  	  MPI_Gather(&Flag_GS_Single,          1, MPI_INT,
//  		     &Flag_GS_cpu_array.at(0), 1, MPI_INT,
//  		     0, MPI_COMM_WORLD);
// 	  if (rank==0) {
// 	    Flag_GS_Total= 0;
// 	    for (int i=0; (i<=NODE); ++i){
// 	      Flag_GS_Total += Flag_GS_cpu_array.at(i);
// 	      //cout <<Flag_GS_Total<<endl;
// 	    }
// 	  }
//    	  MPI_Bcast(&Flag_GS_Total, 1,  MPI_INT, 0, MPI_COMM_WORLD);


   	  P->Gather_Compare_Cast_INT(Flag_GS_Single,Flag_GS_Total,jump_num,P->Debug_mes="--- (1): Poisson Solver"); //jenab removed to test MPI_Gather

      } 
      while (Flag_GS_Total >= 1);
      #ifdef debug_Poisson
	if (rank == 0)cout <<" ((rank ="<<rank<<", i_GS="<<i_GS<<"))  "<<std::flush; 
      #endif
  }// End Poisson_Solver_GS_Parallel  
  //-----------------------------------------------------------------------------------------------
  
  //-----------------------------------------------------------------------------------------------
  void Prolong_Parallel (vector<double>& Func_in,  int& nX_Func_in,
		       vector<double>& Func_out, int& nX_Func_out) {
      
    fill(Func_out.begin(),Func_out.end(),0.0);
    
    P->Ghost_Grid_Points_Vector(Func_in,nX_Func_in,jump_num);
    
    Func_out[2*(nX_Func_in-1)+1] = 0.5 * Func_in[nX_Func_in-1] 
		       + 0.5 * P->ghost_nX_cpu;//CPU Boundary Correction
    
    for (int iX = 0; iX<=nX_Func_in-2; iX++){
      Func_out[2*iX+1] = 0.5 * Func_in[iX] 
		       + 0.5 * Func_in[iX+1];
      }

      for (int iX = 0; iX<=nX_Func_in-1; iX++){
	  Func_out[2*iX]=Func_in[iX];
      }
  }
  //-----------------------------------------------------------------------------------------------
  
  //-----------------------------------------------------------------------------------------------
  void Restrict_RO_Parallel (vector<double>& Func_in,  int& nX_Func_in,
			   vector<double>& Func_out, int& nX_Func_out) {
      // cout <<"------------------ Restrict parallel -----------------multigrid calls="<<P->multigrid_call_num<<"  rank="<<rank<<endl;
      fill(Func_out.begin(),Func_out.end(),0.0);
      for (int iX = 1; iX<=nX_Func_out-2; iX++) {
	  Func_out[iX] = 0.25 * Func_in[iX*2-1]
		       + 0.5  * Func_in[(iX*2)]
		       + 0.25 * Func_in[(iX*2)+1]; 
      }
      
      P->Ghost_Grid_Points_Vector(Func_in,nX_Func_in,jump_num);
      
      Func_out[0] = 0.25 * P->minus_1
		  + 0.5  * Func_in[0] 
		  + 0.25 * Func_in[+1];//CPU Boundary Correction 
   
	       
      Func_out[nX_Func_out-1] = 0.25 * Func_in[nX_Func_in-2]
			      + 0.5  * Func_in[nX_Func_in-1]
			      + 0.25 * P->ghost_nX_cpu;//CPU Boundary Correction

   }
  //-----------------------------------------------------------------------------------------------
 
  //-----------------------------------------------------------------------------------------------
  void Residual_Computation_Parallel(vector<double>& PHI_input,vector<double>& RO_input,
				   int & nX_input, double& dX_input,vector<double>& Resid_out) {
    
    
      // cout <<"-- residual compuation parallel -----multigrid calls="<<P->multigrid_call_num<<"  rank="<<rank<<endl;      
      vector<double> Delta_PHI(nX_input+1);
      fill(Delta_PHI.begin(), Delta_PHI.end(), 0.0);      
  
      P->Ghost_Grid_Points_Vector(PHI_input,nX_input,jump_num);

      
      //----- calculation of Delta_PHI -------------------
      for (int iX = 1; iX<=nX_input-2; iX++) {
	  Delta_PHI[iX] = PHI_input[iX-1]
			- 2.0*PHI_input[iX] 
			+ PHI_input[iX+1];
      }      
      Delta_PHI[0] =  P->minus_1 - 2.0 * PHI_input[0] + PHI_input[1]; //CPU Boundary Correction
      Delta_PHI[nX_input-1] = PHI_input[nX_input-2] - 2.0 * PHI_input[nX_input-1] + P->ghost_nX_cpu; //CPU Boundary Correction
      
      for (int iX = 0; iX<=nX_input-1; iX++){
	  Delta_PHI[iX] = Delta_PHI[iX]/(pow(dX_input,2));
      }
      
      //---------------------------------------------------
      fill(Resid_out.begin(), Resid_out.end(), 0.0);
      for (int iX = 0; iX<=nX_input-1; iX++) {
	  Resid_out[iX] = RO_input[iX] + Delta_PHI[iX];
      }
 

  }//END Residual_Computation_Parallel
  //***********************************************************************************************

   
  
  //***********************************************************************************************
  void Multigrid_Solver_Single (vector<double>& Error_S_ini, vector<double>& Resid_S_ini,
				int & nX_S_ini, double & dX_S_ini )  {
//       int n_relax_down=10, n_relax_up=1000000; 
      int Direct_Grid_Number=4;
      int nX_S_out=nX_S_ini/2;
      double dX_T_s_1 = 2.0*dX_S_ini;
//       double Max_Error_GS=0.000001;
    
      vector<double> Resid_S_mid(nX_S_ini+1),Error_S_cor(nX_S_ini+1);
      vector<double> Error_S_out((nX_S_out)+1),Resid_S_out((nX_S_out)+1);
      vector<double> Error_Tot(nX_S_ini+1); // for printing the 
      
      #ifdef debug_Poisson
	P->Print_Debug_Single(P->line_NO=__LINE__,P->Debug_mes="---------- 0) start of Multigrid_Solver_Single ");
      #endif 
      fill(Resid_S_out.begin(), Resid_S_out.end(), 0.0);
      fill(Error_Tot.begin(), Error_Tot.end(), 0.0);
      fill(Resid_S_mid.begin(), Resid_S_mid.end(), 0.0);
      
      //--- (0) pre relaxation -----
      #ifdef debug_Poisson
	P->Print_Debug_Single(P->line_NO=__LINE__,P->Debug_mes="---------- (0) pre relaxation Poisson_Solver_GS_Single ");
      #endif 
      Poisson_Solver_GS_Single(Error_S_ini,Resid_S_ini,nX_S_ini,dX_S_ini,P->n_relax_down,P->Max_Error_GS);

      //--- (1) computing the residual ----
      
      #ifdef debug_Poisson
	P->Print_Debug_Single(P->line_NO=__LINE__,P->Debug_mes="---------- (1) computing the Residual_Computation_Single ");
      #endif       
      
      Residual_Computation_Single (Error_S_ini,Resid_S_ini,nX_S_ini,dX_S_ini,Resid_S_mid);
   
      //--- (2) restricting the residual ---- 
      Restrict_RO_Single (Resid_S_mid,nX_S_ini,Resid_S_out,nX_S_out);
      
      //--- (3) guess zero for the next level initial value ------
      fill(Error_S_out.begin(), Error_S_out.end(), 0.0);

      //==== (4) next iteration ===============================================
      if(nX_S_ini > Direct_Grid_Number*2) {
	  Multigrid_Solver_Single(Error_S_out,Resid_S_out,nX_S_out,dX_T_s_1);
      } else {
	  //--------------------------------------------------------------------------------------------------------
	  //cout << "--- exact answer on " <<nX_S_out << "-------------------------------"<< endl;
// 	  Eigen::MatrixXd Mat_exact(Direct_Grid_Number,Direct_Grid_Number); // the exact answer is calculated for 4 grid points
	  vector<double> Mat_exact_1D (Direct_Grid_Number*Direct_Grid_Number);
	  #ifdef debug_Poisson
	    P->Print_Debug_Single(P->line_NO=__LINE__,P->Debug_mes="---------- (1) exact answer on ");
	    cout <<Direct_Grid_Number<< " nx of excat  "<<nX_S_ini<<endl;
	  #endif
	  int k;
	  for (int i = 0; i<=nX_S_out-1; i++) {
	      for (int j = 0; j<=nX_S_out-1; j++) {// the exact answer based on documentation
		  k = i * nX_S_out +j;
		  Mat_exact_1D[k] = -1.0 / (pow (M_PI,2)); // elements are -1/pi^2
// 		  Mat_exact(i,j) = -1.0 / (pow (M_PI,2)); // elements are -1/pi^2
		  if (i==j){
// 		    Mat_exact(i,j) = -3.0 * Mat_exact(i,j); // diognal elements are 3/pi^2
		    Mat_exact_1D[k] = -3.0 * Mat_exact_1D[k]; // diognal elements are 3/pi^2		    
		  }
	      }
	  }
	  int k_exact;
	  for (int k = 0; k<=nX_S_out-1; k++) {
	      for (int iX = 0; iX<=nX_S_out-1; iX++) {// the exact answer based on documentation
		  k_exact = k * nX_S_out + iX;
 		  Error_S_out[k] = Mat_exact_1D[k_exact] * Resid_S_out[iX]*dX_T_s_1*dX_T_s_1; 
// 		  Error_S_out[k] = Mat_exact(k,iX) * Resid_S_out[iX]*dX_T_s_1*dX_T_s_1;
	      }
	  }
	  Error_S_out[nX_S_out] = Error_S_out[0]; //bourder correction on the exact level
	  //--------------------------------------------------------------------------------------------------------
	  //Multigrid_PHI_Printing_Single(Error_S_out,nX_S_out,dX_S_ini*2.0D0,Error_Tot,nX_S_out)
      }
      //=======================================================================
      
      
      
      //--- (5) prolon the Error or PHI ------------------
      fill(Error_S_cor.begin(), Error_S_cor.end(), 0.0);
      Prolong_Single (Error_S_out,nX_S_out,Error_S_cor,nX_S_ini);
   
      //--- (6) adding the correction to the Error or PHI --------
      for (int iX = 0; iX<=nX_S_ini; iX++) {
	  Error_S_ini[iX]=-Error_S_ini[iX] + Error_S_cor[iX];
      }
      //Multigrid_PHI_Printing_Single(Error_S_ini,nX_S_ini,dX_S_ini,Error_Tot,nX_S_ini)
 
      //--- (7) post relaxation -------------------------------------------------------
      Poisson_Solver_GS_Single(Error_S_ini,Resid_S_ini,nX_S_ini,dX_S_ini,P->n_relax_up,P->Max_Error_GS);
      //Multigrid_PHI_Printing_Single(Error_S_ini,nX_S_ini,dX_S_ini,Error_Tot,nX_S_ini)*/
      
      #ifdef debug_Poisson
	    P->Print_Debug_Single(P->line_NO=__LINE__,P->Debug_mes="---------- End of Multigrid_Solver_Single ");
      #endif 
            
  }// End Multigrid_Solver_Single 
  //-----------------------------------------------------------------------------------------------

  //-----------------------------------------------------------------------------------------------
  void Poisson_Solver_GS_Single (vector<double>& PHI_X_ini,vector<double>& RO_X_ini
			      ,int& nX_GS,double& dX_GS,int& Max_Loop,double& Max_Error_GS){
  

      double Resid_GS_2,Resid_GS_1;
      bool Flag_GS_Single;
      int i_GS;
      vector<double>  Delta_PHI(nX_GS+1);
      fill(Delta_PHI.begin(), Delta_PHI.end(), 0.0);
   
      Resid_GS_1=0.;
      Resid_GS_2=0.;
   
      i_GS=0;
   
      Flag_GS_Single = true;
   
      while (Flag_GS_Single == true) {
	  //cout<< "-------------------------- Jacobi method of solving Poisson equation ----------------"<<endl;
	  i_GS += 1;
	  PHI_X_ini[0] = PHI_X_ini[nX_GS-1]+ PHI_X_ini[0+1];
	  PHI_X_ini[0] = PHI_X_ini[0] - (RO_X_ini[0]*( pow(dX_GS,2)));
	  PHI_X_ini[0] = 0.5 * PHI_X_ini[0];
	  
	  PHI_X_ini[nX_GS-1] = PHI_X_ini[nX_GS-1-1]+ PHI_X_ini[0]; 
	  PHI_X_ini[nX_GS-1] = PHI_X_ini[nX_GS-1] - (RO_X_ini[nX_GS-1]*( pow(dX_GS,2)));
	  PHI_X_ini[nX_GS-1] = 0.5 * PHI_X_ini[nX_GS-1];
	      
	  
	  for (int iX = 1; iX<=nX_GS-2; iX++) {//jenab changed from  for (int iX = 1; iX<=nX_GS-1; iX++)
	      PHI_X_ini[iX] = PHI_X_ini[iX-1]+ PHI_X_ini[iX+1]; 
	      PHI_X_ini[iX] = PHI_X_ini[iX] - (RO_X_ini[iX]*( pow(dX_GS,2)));
	      PHI_X_ini[iX] = 0.5 * PHI_X_ini[iX];
	  }
//  	  PHI_X_ini(nX_GS) = PHI_X_ini(0);
	  //-------------------------------------------------------------------------------------
	  
	  
	  //cout<<"----- calculation of Delta_PHI ------------------------------------------------------"<<endl;
	  for (int iX = 1; iX<=nX_GS-2; iX++) {
	      Delta_PHI[iX] = PHI_X_ini[iX-1]-2.0*PHI_X_ini[iX]+PHI_X_ini[iX+1];
	  }
	  Delta_PHI[0] =  PHI_X_ini[nX_GS-1] - 2.0*PHI_X_ini[0] + PHI_X_ini[1];
// 	  Delta_PHI(nX_GS) = PHI_X_ini[nX_GS-1]-2.0*PHI_X_ini(nX_GS)+PHI_X_ini[1]; //jenab deleted for new approach in boundary condition 12 Feb
	  Delta_PHI[nX_GS-1] = PHI_X_ini[nX_GS-2]-2.0*PHI_X_ini[nX_GS-1]+PHI_X_ini[0]; //jenab 

	  for (int iX = 0; iX<=nX_GS-1; iX++) {
	      Delta_PHI[iX] = Delta_PHI[iX]/( pow(dX_GS,2));
	  }
	  //-------------------------------------------------------------------------------------
      
      
	  //cout<<"---- Finding the maximum residual -------------"<<endl;
	  Resid_GS_2 = 0.0;
	  for (int iX = 0; iX<=nX_GS; iX++)  {
	      if ( abs( RO_X_ini[iX] + Delta_PHI[iX] )> Resid_GS_2 )
	      Resid_GS_2  = abs( RO_X_ini[iX] + Delta_PHI[iX] );
	  }
	  //--------------------------------------------
      
	  //cout<<"---- comparing the maximum error with the acceptable error -----------------------------"<<endl;
	  if ( i_GS > 2) { // this if makes sure that the initial stage is not going to be checked which resutl in division by zero
	      if ( abs((Resid_GS_2-Resid_GS_1)/Resid_GS_1) < Max_Error_GS) {
		  Flag_GS_Single = false;
		  //Print*,"converged",i_GS,rank,Time
	      } else if(i_GS>Max_Loop) {
		  Flag_GS_Single = false;
		  //  Print*,"Not converging",rank,Time
	      }	  	
	  }
	  Resid_GS_1 = Resid_GS_2;
	  //-----------------------------------------------
	  //cout <<i_GS<<endl;
	  //--------------------------------------------------------------------------------------------------
      }
      //cout << "level "<<nX_GS<< "    number of iteration= "<<i_GS<<endl;
  }// End Poisson_Solver_GS_Single  
  //-----------------------------------------------------------------------------------------------
  
  //-----------------------------------------------------------------------------------------------
  void Prolong_Single (vector<double>& Func_in,  int& nX_Func_in,
		       vector<double>& Func_out, int& nX_Func_out) {
      fill(Func_out.begin(),Func_out.end(),0.0);
      Func_in[nX_Func_in]=Func_in[0];
      for (int iX = 0; iX<=nX_Func_in-1; iX++){
      Func_out[2*iX+1] = 0.5 * Func_in[iX] 
		       + 0.5 * Func_in[iX+1];
      }

      for (int iX = 0; iX<=nX_Func_in-1; iX++){
	  Func_out[2*iX]=Func_in[iX];
      }
  }
  //-----------------------------------------------------------------------------------------------
  
  //-----------------------------------------------------------------------------------------------
  void Restrict_RO_Single (vector<double>& Func_in,  int& nX_Func_in,
			   vector<double>& Func_out, int& nX_Func_out) {
      fill(Func_out.begin(),Func_out.end(),0.0);
      for (int iX = 1; iX<=nX_Func_out-1; iX++) {
	  Func_out[iX] = 0.25 * Func_in[iX*2-1]
		       + 0.5  * Func_in[(iX*2)]
		       + 0.25 * Func_in[(iX*2)+1]; 
      }

      Func_out[0] = 0.25 * Func_in[nX_Func_in-1]
		  + 0.5  * Func_in[0] 
		  + 0.25 * Func_in[+1]; 
   
	       
      Func_out[nX_Func_out] = 0.25 * Func_in[nX_Func_in-1]
			    + 0.5  * Func_in[0] 
			    + 0.25 * Func_in[1];
  
   }
  //-----------------------------------------------------------------------------------------------
 
  //-----------------------------------------------------------------------------------------------
  void Residual_Computation_Single(vector<double>& PHI_input,vector<double>& RO_input,
				   int & nX_input, double& dX_input,vector<double>& Resid_out) {
      vector<double> Delta_PHI(nX_input+1);
      fill(Delta_PHI.begin(), Delta_PHI.end(), 0.0);
      fill(Resid_out.begin(), Resid_out.end(), 0.0);


      //----- calculation of Delta_PHI -------------------
      for (int iX = 1; iX<=nX_input-2; iX++) {
	  Delta_PHI[iX] = PHI_input[iX-1]
			- 2.0*PHI_input[iX] 
			+ PHI_input[iX+1];
      }
      Delta_PHI[0] =  PHI_input[nX_input-1] 
		   - 2.0 * PHI_input[0] 
		   + PHI_input[1];
		 
      Delta_PHI[nX_input] = PHI_input[nX_input-1] 
			  - 2.0 * PHI_input[nX_input] 
			  + PHI_input[1];
			
      Delta_PHI[nX_input-1] = PHI_input[nX_input-2] 
			    - 2.0 * PHI_input[nX_input-1]
			    + PHI_input[nX_input];
  
      for (int iX = 0; iX<=nX_input; iX++){
	  Delta_PHI[iX] = Delta_PHI[iX]/(pow(dX_input,2));
      }
      
      //---------------------------------------------------
      for (int iX = 0; iX<=nX_input; iX++) {
	  Resid_out[iX] = RO_input[iX] + Delta_PHI[iX];
      }

  }//END Residual_Computation_Single
  //***********************************************************************************************
 
	//***********************************************************************************************
	void E_x_Calculator_2points() {
		fill(E_X.begin(), E_X.end(), 0.0);
		P->Ghost_Grid_Points_Vector(PHI_X,P->nX_cpu,1);

		E_X[0]              = -P->minus_1           +  PHI_X [1];   
		for (int iX = 1; iX<P->nX_cpu-1; iX++){
			E_X[iX] = -PHI_X [ iX-1 ] +  PHI_X [ iX+1 ];
		}
		E_X[P->nX_cpu-1]    = -PHI_X[P->nX_cpu-2]   +  P->ghost_nX_cpu;  
		for (int iX = 0; iX<P->nX_cpu; iX++){
			E_X[iX] = -E_X[iX]/ (2.0*P->dX);
		}
	} //End E_x_Calculator
	//***********************************************************************************************
  
  //***********************************************************************************************
  void E_x_Calculator() {
// 	  f_x = (1*f[i-2]-8*f[i-1]+0*f[i+0]+8*f[i+1]-1*f[i+2])/(12*1.0*h**1)
    fill(E_X.begin(), E_X.end(), 0.0);
    P->Ghost_Grid_Points_Vector(PHI_X,P->nX_cpu,1);

    E_X[0] = -  (-PHI_X [2] + 8.0* PHI_X [1] - 8.0* P->minus_1 + P->minus_2) / ( 12.0*P->dX );//CPU Boundary Correction
    E_X[1] = -  (-PHI_X [3] + 8.0* PHI_X [2] - 8.0* PHI_X [0] + P->minus_1) / ( 12.0*P->dX );//CPU Boundary Correction
    
    E_X[P->nX_cpu-1] = -  (-P->nX_cpu_plus_1 + 8.0* P->ghost_nX_cpu 
			  - 8.0* PHI_X[P->nX_cpu-2] + PHI_X[P->nX_cpu-3])
			/ ( 12.0*P->dX );//CPU Boundary Correction
			
    E_X[P->nX_cpu-2] = -  (-P->ghost_nX_cpu + 8.0* PHI_X[P->nX_cpu-1] 
			  - 8.0* PHI_X[P->nX_cpu-3] + PHI_X[P->nX_cpu-4])
			/ ( 12.0*P->dX );  //CPU Boundary Correction    
      
    for (int iX = 2; iX<=P->nX_cpu-3; iX++) //with four neighbouring PhasePoints1
     {
      E_X[iX] = -  (-PHI_X [ (iX+1)+1 ] 
		     + 8.0* PHI_X [ iX+1 ]
		     - 8.0* PHI_X [ iX-1 ]
		     + PHI_X [  (iX-1)-1 ]
							) / ( 12.0*P->dX );
      
                      // Five Point Stencil

     }
  } //End E_x_Calculator
  //***********************************************************************************************
 
	//*********************************************************************
	void Zalesak_limiter_fourth_order(vector<double> & array_in) {
		/* this is coming from Kazeminezhad paper PRE  67, 026704, 2003. 
		 * which is cited 
		 * Journal of Computational Physics Volume 31, Issue 3, June 1979, Pages 335-362
		 * as the main source for this equation. Paper by Zalesak, Steven T.P. has this 
		 * equation in its appendix and proposed higher order versions of it. 
		 * 2th order  enhanced [i] =   1/2   (f[i+1]+ f[i]) 
		 * 4th order  enhanced [i] =   7/12  (f[i+1]+ f[i]) -   1/12  (f[i+2]+f[i-1]) 
		 * 6th order  enhanced [i] =  37/60  (f[i+1]+ f[i]) -   2/15  (f[i+2]+f[i-1]) + 1/60   (f[i+3]+f[i-2])
		 * 8th order  enhanced [i] = 533/840 (f[i+1]+ f[i]) - 139/840 (f[i+2]+f[i-1]) + 29/840 (f[i+3]+f[i-2]) - 1/280 (f[i+4]+f[i-3])
		 */
		vector<double> array_enhanced(P->nX_cpu);
		double a_coef,b_coef;
		fill(array_enhanced.begin(), array_enhanced.end(), 0.0);
		P->Ghost_Grid_Points_Vector(array_in,P->nX_cpu,1);
		a_coef = 7.0/12.0; b_coef = -1.0/12.0;
		
		
		array_enhanced[0]           = a_coef * (array_in[1]           + array_in[0]          ) + b_coef * ((array_in[2]         + P->minus_1        ));
		array_enhanced[P->nX_cpu-1] = a_coef * (P->ghost_nX_cpu    + array_in[P->nX_cpu-1]) + b_coef * ((P->nX_cpu_plus_1 + array_in[P->nX_cpu-2]));
		array_enhanced[P->nX_cpu-2] = a_coef * (array_in[P->nX_cpu-1] + array_in[P->nX_cpu-2]) + b_coef * ((P->ghost_nX_cpu  + array_in[P->nX_cpu-3]));

		for (int iX = 1; iX<P->nX_cpu-2; ++iX) {
			array_enhanced[iX]  = a_coef * (array_in[iX+1]        + array_in[iX]         ) + b_coef * ((array_in[iX+2]      + array_in[iX-1]        ));

		}
		for (int iX = 0; iX<P->nX_cpu; ++iX) {
			array_in[iX] = array_enhanced[iX];
		}
	} 
	//*********************************************************************

  //***********************************************************************************************
  void Ex_Ghost_Talk()  {
    
//     double E_X_snd, E_X_rcv;
    
    
    for (int i=0; i<=P->nX_cpu-1; ++i){
      E_X_ghost[i+1]=E_X[i]; //since there is no way to start the array from -1, then we have to push forward one element so we cap assign Ex(-1) to E_X_ghost(0)
    }
    
    P->Ghost_Grid_Points_Vector(E_X,P->nX_cpu,1);
//     E_X_snd=E_X(P->nX_cpu-1);
//     MPI_Sendrecv(&E_X_snd, 1, MPI_DOUBLE, Right_cpu[rank], 233,
// 		 &E_X_rcv, 1, MPI_DOUBLE, Left_cpu[rank],  233,
// 		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    E_X_ghost[0]=P->minus_1;// E_X_rcv; //this is -1 of Ex
    
//     E_X_snd=E_X(1);
//     MPI_Sendrecv(&E_X_snd, 1, MPI_DOUBLE, Left_cpu[rank],  233,
// 		 &E_X_rcv, 1, MPI_DOUBLE, Right_cpu[rank], 233,
// 		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     E_X_ghost[P->nX_cpu+2]=P->nX_cpu_plus_1;// E_X_rcv; //this is nX_cpu+1
    
    
//     E_X_snd=E_X(0);
//     MPI_Sendrecv(&E_X_snd, 1, MPI_DOUBLE, Left_cpu[rank],  233,
// 		 &E_X_rcv, 1, MPI_DOUBLE, Right_cpu[rank], 233,
// 		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    E_X_ghost[P->nX_cpu+1]=P->ghost_nX_cpu;  // E_X_rcv; //this is nX_cpu
    
    
  } 
  //***********************************************************************************************
 
  //***********************************************************************************************
  void find_Ep( double& X,double& Ex_p)  {
      
      int lX;
      double Ex_1,Ex_2,Ex_3,Ex_4;
      double X_0,X_1,X_2,X_3,X_4;
      lX = int ( (X-P->X_min) / P->dX );
      if (lX<0 || lX>=P->nX_cpu) cout <<"error: out of indexing border in X direction, rank"<<rank
	<<" X="<<X<<endl;
      
      X_0=X-P->X_min;      
      X_2=lX*P->dX;
      X_1=X_2-P->dX;  
      X_3=X_2+P->dX;  
      X_4=X_2+2*P->dX;

      /*Ex_1=E_X((lX)-1);
      Ex_2=E_X(lX);
      Ex_3=E_X((lX)+1);
      Ex_4=E_X(((lX)+1)+1);*/
      
      Ex_1=E_X_ghost[(lX+1)-1]; //+1 has to be added since the E_X_ghost is shifted ahead compared to E_X
      Ex_2=E_X_ghost[lX+1];
      Ex_3=E_X_ghost[(lX+1)+1];
      Ex_4=E_X_ghost[((lX+1)+1)+1];
      
      
      
      Ex_p = (X_0-X_2)*(X_0-X_3)*(X_0-X_4)*Ex_1 /(-6)         
            +(X_0-X_1)*(X_0-X_3)*(X_0-X_4)*Ex_2 /(+2)      
            +(X_0-X_2)*(X_0-X_1)*(X_0-X_4)*Ex_3 /(-2)      
            +(X_0-X_2)*(X_0-X_3)*(X_0-X_1)*Ex_4 /(+6);
      Ex_p = Ex_p/(P->dX*P->dX*P->dX);
  } //End find_Ep
  //***********************************************************************************************
  
  //***********************************************************************************************
  void Ene_ele_Calculator() {
    Energy_ele_cpu ();
    if (size>1){
        P->Sum_up_Double(Ene_ele_cpu,Ene_ele_Tot,1);
    }else{
        Ene_ele_Tot = Ene_ele_cpu;
    }
  } 
  //***********************************************************************************************
  
  //***********************************************************************************************
  void Energy_ele_cpu() {

    Ene_ele_cpu = 0.0;
         
       P->Ghost_Grid_Points_Vector(E_X,P->nX_cpu,1);

    
    Ene_ele_cpu = Ene_ele_cpu + ((
				   pow(P->minus_1,2)
			      + 4* pow(E_X[0]    ,2)
			      +    pow(E_X[1]    ,2)
				    ) /3.0) * P->dX; //CPU Boundary Correction
				    
    Ene_ele_cpu = Ene_ele_cpu + ((
				   pow(E_X[P->nX_cpu-2],2)
			      + 4* pow(E_X[P->nX_cpu-1],2)
			      +    pow(P->ghost_nX_cpu ,2)
				    ) /3.0) * P->dX; //CPU Boundary Correction    
    for (int iX = 1; iX<=P->nX_cpu-2; iX++) {
      Ene_ele_cpu = Ene_ele_cpu + (( 
				   pow(E_X[iX-1],2)
			      + 4* pow(E_X[iX]  ,2)
			      +    pow(E_X[iX+1],2)
				    ) /3.0) * P->dX;
     }
      Ene_ele_cpu = 0.5 * Ene_ele_cpu;

     
  } 
  //***********************************************************************************************
     
  //***********************************************************************************************
  void Write_Temporal_Ex_Field_Techplot() {
      file_Time_Ex << P->R_Time <<"   "<< E_X[0]<< endl;
      file_Time_Ex.flush();      
  } 
  //***********************************************************************************************
  
  //***********************************************************************************************
  void Write_Temporal_Ene_ele_Techplot() {
      file_Time_Ene_ele << P->R_Time <<"   "<< Ene_ele_Tot<< endl;
      file_Time_Ene_ele.flush();
      
  } 
  //***********************************************************************************************
   
 
  
  //***********************************************************************************************
  void  Write_Dset_RO_Ex_PHI_HDF5(){    
    

    //-------------------------------------------------------------------------
   
    // 2) Create the dataspace for the dataset  -------------------------------
    P->Space_ID_RO = H5Screate_simple(P->Rank_X, P->Dim_X , NULL); 
    P->Memory_ID_RO = H5Screate_simple(P->Rank_X, P->Dim_X_Chunk, NULL);
    //-------------------------------------------------------------------------
   
    // 3) Create chunked dataset ----------------------------------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(P->PList_ID, P->Rank_X, P->Dim_X_Chunk);

    P->DSet_ID_RO = H5Dcreate(P->Group_ID_S, "RO", H5T_NATIVE_DOUBLE, P->Space_ID_RO,
			  H5P_DEFAULT, P->PList_ID, H5P_DEFAULT);
    H5Pclose(P->PList_ID);
    //------------------------------------------------------------------------
    H5Sclose(P->Space_ID_RO);
  
    // 4) Select hyperslab in the file----------------------------------------
    P->Space_ID_RO = H5Dget_space(P->DSet_ID_RO);
    P->status_hdf5 = H5Sselect_hyperslab(P->Space_ID_RO, H5S_SELECT_SET, P->OffSet_X, P->Stride_X, P->Count_X, P->Block_X);
    //-----------------------------------------------------------------------

    // 5) Create property list for collective dataset write------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(P->PList_ID, H5FD_MPIO_COLLECTIVE);
    //-----------------------------------------------------------------------
    
    P->status_hdf5 = H5Dwrite(P->DSet_ID_RO, H5T_NATIVE_DOUBLE, P->Memory_ID_RO, P->Space_ID_RO,
 		      P->PList_ID, RO_X.data());

   
    //H5Fflush(File_ID_S, H5F_SCOPE_GLOBAL ); 
    
    // Close/release resources    
    H5Dclose(P->DSet_ID_RO);
    H5Pclose(P->PList_ID);
    
    H5Sclose(P->Space_ID_RO);
    H5Sclose(P->Memory_ID_RO);
    
    
    
    
    //-------------------------------------------------------------------------
   
    // 2) Create the dataspace for the dataset  -------------------------------
    P->Space_ID_E = H5Screate_simple(P->Rank_X, P->Dim_X , NULL); 
    P->Memory_ID_E = H5Screate_simple(P->Rank_X, P->Dim_X_Chunk, NULL);
    //-------------------------------------------------------------------------
   
    // 3) Create chunked dataset ----------------------------------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(P->PList_ID, P->Rank_X, P->Dim_X_Chunk);

    P->DSet_ID_E = H5Dcreate(P->Group_ID_S, "Ex", H5T_NATIVE_DOUBLE, P->Space_ID_E,
			  H5P_DEFAULT, P->PList_ID, H5P_DEFAULT);
    H5Pclose(P->PList_ID);
    //------------------------------------------------------------------------
    H5Sclose(P->Space_ID_E);
  
    // 4) Select hyperslab in the file----------------------------------------
    P->Space_ID_E = H5Dget_space(P->DSet_ID_E);
    P->status_hdf5 = H5Sselect_hyperslab(P->Space_ID_E, H5S_SELECT_SET, P->OffSet_X, P->Stride_X, P->Count_X, P->Block_X);
    //-----------------------------------------------------------------------

    // 5) Create property list for collective dataset write------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(P->PList_ID, H5FD_MPIO_COLLECTIVE);
    //-----------------------------------------------------------------------
    
    P->status_hdf5 = H5Dwrite(P->DSet_ID_E, H5T_NATIVE_DOUBLE, P->Memory_ID_E, P->Space_ID_E,
		      P->PList_ID, E_X.data());
   
    //H5Fflush(File_ID_S, H5F_SCOPE_GLOBAL ); 
    
    // Close/release resources    
    H5Dclose(P->DSet_ID_E);
    H5Pclose(P->PList_ID);
    
    H5Sclose(P->Space_ID_E);
    H5Sclose(P->Memory_ID_E);
    
    
    
    //-------------------------------------------------------------------------
   
    // 2) Create the dataspace for the dataset  -------------------------------
    P->Space_ID_PHI = H5Screate_simple(P->Rank_X, P->Dim_X , NULL); 
    P->Memory_ID_PHI = H5Screate_simple(P->Rank_X, P->Dim_X_Chunk, NULL);
    //-------------------------------------------------------------------------
   
    // 3) Create chunked dataset ----------------------------------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(P->PList_ID, P->Rank_X, P->Dim_X_Chunk);

    P->DSet_ID_PHI = H5Dcreate(P->Group_ID_S, "PHI", H5T_NATIVE_DOUBLE, P->Space_ID_PHI,
			  H5P_DEFAULT, P->PList_ID, H5P_DEFAULT);
    H5Pclose(P->PList_ID);
    //------------------------------------------------------------------------
    H5Sclose(P->Space_ID_PHI);
  
    // 4) Select hyperslab in the file----------------------------------------
    P->Space_ID_PHI = H5Dget_space(P->DSet_ID_PHI);
    P->status_hdf5 = H5Sselect_hyperslab(P->Space_ID_PHI, H5S_SELECT_SET, P->OffSet_X, P->Stride_X, P->Count_X, P->Block_X);
    //-----------------------------------------------------------------------

    // 5) Create property list for collective dataset write------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(P->PList_ID, H5FD_MPIO_COLLECTIVE);
    //-----------------------------------------------------------------------
    
    P->status_hdf5 = H5Dwrite(P->DSet_ID_PHI, H5T_NATIVE_DOUBLE, P->Memory_ID_PHI, P->Space_ID_PHI,
		      P->PList_ID, PHI_X.data());
   
    //H5Fflush(File_ID_S, H5F_SCOPE_GLOBAL ); 
    
    // Close/release resources    
    H5Dclose(P->DSet_ID_PHI);
    H5Pclose(P->PList_ID);
    
    H5Sclose(P->Space_ID_PHI);
    H5Sclose(P->Memory_ID_PHI);
   
    }
  //***********************************************************************************************

    
  };
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#endif //GridClass_Included
