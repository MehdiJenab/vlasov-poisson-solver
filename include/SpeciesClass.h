#ifndef SpeciesClass_Included
#define SpeciesClass_Included




//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
class DistributionFunctionClass{
public:
	  ParametersClass * P;
	//---- MPI variables ------
  	int size, rank,NODE;
	int * Left_cpu;
	int * Right_cpu;
	

	double Phi_imposed, Phi_imposed_CoMoving_Correction, Velocity_imposed, Beta_imposed, Phi_maximum_imposed, Alpha_Schamel_Main_Soliton,v_soliton_Other_Soliton,phi_max_Other_Soliton,Alpha_Schamel_Other_Soliton;
	
	bool Schamel_DF_Flag = false;

	vector<double> Potential_imposed_array_cpu,Potential_imposed_array_total, v_soliton_imposed_array_cpu, Beta_soliton_imposed_array_cpu, Phi_maximum_soliton_imposed_array_cpu, Alpha_Schamel_Main_Soliton_array_cpu,v_soliton_Other_Soliton_array_cpu,phi_max_Other_Soliton_array_cpu,Alpha_Schamel_Other_Soliton_array_cpu;
	
	vector<double> Energy_Kinentic_Jenab_DF, Beta_Jenab_DF;
	double gamma1 = 0.0, gamma2 = 0.0;
	double X0_CoMoving_correction; int iX_CoMoving_correction;
	
	
	
	//=========================================================================
	DistributionFunctionClass (ParametersClass * P_in, MpiClass * MPI){
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
		
		v_soliton_imposed_array_cpu.resize (P->nX_cpu); 
		Beta_soliton_imposed_array_cpu.resize (P->nX_cpu);
		Phi_maximum_soliton_imposed_array_cpu.resize (P->nX_cpu);
		Alpha_Schamel_Main_Soliton_array_cpu.resize (P->nX_cpu);
		v_soliton_Other_Soliton_array_cpu.resize (P->nX_cpu);
		phi_max_Other_Soliton_array_cpu.resize (P->nX_cpu);
		Alpha_Schamel_Other_Soliton_array_cpu.resize (P->nX_cpu);
		
		Potential_imposed_array_cpu.resize (P->nX_cpu);
		Potential_imposed_array_total.resize (P->nX);
		
		

		
		initialization_DF_Kappa(gamma1,gamma2);
		
		
		
		fill(Potential_imposed_array_cpu.begin(), Potential_imposed_array_cpu.end(), 0.0); 
		fill(Potential_imposed_array_total.begin(), Potential_imposed_array_total.end(), 0.0);
		
		fill(v_soliton_imposed_array_cpu.begin(), v_soliton_imposed_array_cpu.end(), 0.0);
		fill(Beta_soliton_imposed_array_cpu.begin(), Beta_soliton_imposed_array_cpu.end(), 0.0);
		fill(Phi_maximum_soliton_imposed_array_cpu.begin(), Phi_maximum_soliton_imposed_array_cpu.end(), 0.0);
		fill(Alpha_Schamel_Main_Soliton_array_cpu.begin(), Alpha_Schamel_Main_Soliton_array_cpu.end(), 0.0);
		fill(v_soliton_Other_Soliton_array_cpu.begin(), v_soliton_Other_Soliton_array_cpu.end(), 0.0);
		fill(phi_max_Other_Soliton_array_cpu.begin(), phi_max_Other_Soliton_array_cpu.end(), 0.0);
		fill(Alpha_Schamel_Other_Soliton_array_cpu.begin(), Alpha_Schamel_Other_Soliton_array_cpu.end(), 0.0);
		
		if (P->DF_NotShifted_Shifted == 1) {
			Read_potential_imposed (Potential_imposed_array_cpu,Potential_imposed_array_total,v_soliton_imposed_array_cpu,Beta_soliton_imposed_array_cpu, Phi_maximum_soliton_imposed_array_cpu, Alpha_Schamel_Main_Soliton_array_cpu,v_soliton_Other_Soliton_array_cpu,phi_max_Other_Soliton_array_cpu,Alpha_Schamel_Other_Soliton_array_cpu );
		}
		
		
		#ifdef multiple_beta_for_trapped_DF
		Read_Beta_array_for_Jenab_DF(Energy_Kinentic_Jenab_DF, Beta_Jenab_DF);
		//if (rank == 0){cout<<P->Name<<"*************  Energy="<<Energy_Kinentic_Jenab_DF[1] <<"  DF_Jenab="<< Beta_Jenab_DF[2]<<endl<<flush;}
		#endif
		#ifdef single_beta_for_trapped_DF
		Beta_Jenab_DF.resize (2);
		Energy_Kinentic_Jenab_DF.resize (2);
		for (int i_Beta = 0; i_Beta <= 1; i_Beta++){
			Beta_Jenab_DF[i_Beta] = 0.0;
			Energy_Kinentic_Jenab_DF[i_Beta]=0.0; //dummy arrays
		}
		#endif
			
			
	};
	//=========================================================================
	
	
	  #ifdef multiple_beta_for_trapped_DF
	//***********************************************************************************************
	void Read_Beta_array_for_Jenab_DF(vector<double>& Energy_Kinentic_Jenab_DF, vector<double>& Beta_Jenab_DF) {
		/*take note: arrays are starting from 0, 
		and if the pointers goes beyond the domain of the array, no error will be shown
		the type of the element has to agree with the type of the function */
		size_t flags = 0;
		json_error_t error;
		int i_element;
		json_t* Pair_Energy_Beta_json, *Beta_json, *Energy_Kinentic_json;
// 		vector<double> Beta_Jenab_DF, Energy_Kinentic_Jenab_DF;
		int len_of_array;
		
		
		if (rank==0) cout << "-----------------------------------------------------"<<endl;
		
		
		//--- loading the JSON file ------------------------------------------
		json_t* file_json = json_load_file("./Jenab_Distribution_Function.json", flags, &error);
		if(!file_json){
			cerr << "in line " << error.line << ": " << error.text << endl;
			exit(20);
		}
		//--------------------------------------------------------------------
		
		//--- returning the data object ------------------------------------  
		json_t* data_json = json_object_get(file_json, "data");
		if(!data_json){
			cerr << "in line " << error.line << ": " << error.text << endl;
			exit(1);
		}
		//------ reading the length of the array ---------------
		json_t* element_json      = json_array_get(data_json, 0);
		if (!element_json) cout<<"can't catch element len_of_array object of " <<endl;
		element_json = json_object_get(element_json, "len_of_array");
		if (!element_json) cout<<"can't catch len_of_array of "  <<endl;
		len_of_array = json_integer_value(element_json);
		if (rank==0) cout << "len_of_array of " <<" = " << len_of_array << endl;
		//--------------------------------------------------------------------
		Beta_Jenab_DF.resize (len_of_array);  Energy_Kinentic_Jenab_DF.resize (len_of_array);
		
		
		//--- returning the DF object --------------------------------  
		json_t* DF_json = json_array_get(data_json, 1);
		if (!DF_json) cout<<"can't catch DF in {}" <<endl;
		DF_json = json_object_get(DF_json, "DF");
		if (!DF_json) cout<<"can't catch DF object" <<endl;
		//--------------------------------------------------------------------

		i_element = 0;
		//---------------------
		for (int i_Beta = 0; i_Beta <= len_of_array-1; i_Beta++){
			Pair_Energy_Beta_json    = json_array_get(DF_json, i_element);
			if (!Pair_Energy_Beta_json) cout<<"can't catch element "<< i_element<<" of Jenab_DF" <<flush<<endl;
			Energy_Kinentic_json = json_array_get(Pair_Energy_Beta_json, 0);
			if (!Energy_Kinentic_json) cout<<"can't catch Energy_Kinentic_json" <<flush<<endl;
			double Energy_Kinentic_value = json_real_value(Energy_Kinentic_json);
// 			if (rank==0) cout << "Energy_Kinentic_json " <<" = " << Energy_Kinentic_value << endl;
			Beta_json = json_array_get(Pair_Energy_Beta_json, 1);
			if (!Beta_json) cout<<"can't catch Beta_json" <<flush<<endl;
			double Beta_value = json_real_value(Beta_json);
// 			if (rank==0) cout << "Beta value" <<" = " << Beta_value << endl;
			Energy_Kinentic_Jenab_DF [i_Beta] = Energy_Kinentic_value;
			Beta_Jenab_DF[i_Beta] = Beta_value;
			++i_element;
		}
		if (rank==0) cout << "-----------------------------------------------------"<<endl;
	}
	//***********************************************************************************************
	#endif
  
  //***********************************************************************************************
  //-----------------------------------------------------------------------------------------------
  void Read_potential_imposed(vector<double>& Potential_imposed_array_cpu, vector<double>& Potential_imposed_array_total, vector<double>& v_soliton_imposed_array_cpu, vector<double>& Beta_soliton_imposed_array_cpu, vector<double>& Phi_maximum_soliton_imposed_array_cpu, 
      vector<double>& Alpha_Schamel_Main_Soliton_array_cpu,vector<double>& v_soliton_Other_Soliton_array_cpu,
      vector<double>& phi_max_Other_Soliton_array_cpu,vector<double>& Alpha_Schamel_Other_Soliton_array_cpu) {
	#ifdef debug_initial
		P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="$$$ Read_potential_imposed, Started $$$");
	#endif
      /*take note: arrays are starting from 0, 
      and if the pointers goes beyond the domain of the array, no error will be shown
      the type of the element has to agree with the type of the function */
    if (rank == 0 ){
      size_t flags = 0;
      json_error_t error;
      int i_element, X_integer;
      json_t* Pair_X_Potential_json, *Potential_json, *X_json, *velocity_soliton, *Beta_soliton, *Phi_maximum_soliton;
      json_t* Alpha_Schamel_Main_Soliton_json, *v_soliton_Other_Soliton_json, *phi_max_Other_Soliton_json, *Alpha_Schamel_Other_Soliton_json;
      json_t* file_json;
      
      vector<double> Potential_imposed_total, v_soliton_imposed_total, Beta_imposed_total, Phi_maximum_imposed_total; 
      
      vector<double> Alpha_Schamel_Main_Soliton_total, v_soliton_Other_Soliton_total, phi_max_Other_Soliton_total, Alpha_Schamel_Other_Soliton_total;
      
      Potential_imposed_total.resize (P->nX);  v_soliton_imposed_total.resize (P->nX);  Beta_imposed_total.resize (P->nX); Phi_maximum_imposed_total.resize (P->nX); 
      
      Alpha_Schamel_Main_Soliton_total.resize (P->nX); v_soliton_Other_Soliton_total.resize (P->nX);
      phi_max_Other_Soliton_total.resize (P->nX); Alpha_Schamel_Other_Soliton_total.resize (P->nX);
      
      
      fill(Potential_imposed_total.begin(), Potential_imposed_total.end(), 1.0);
      fill(v_soliton_imposed_total.begin(), v_soliton_imposed_total.end(), 0.0);
      fill(Beta_imposed_total.begin(), Beta_imposed_total.end(), 0.0);
      fill(Phi_maximum_imposed_total.begin(), Phi_maximum_imposed_total.end(), 0.0);
      fill(Alpha_Schamel_Main_Soliton_total.begin(), Alpha_Schamel_Main_Soliton_total.end(), 0.0);
      fill(v_soliton_Other_Soliton_total.begin(), v_soliton_Other_Soliton_total.end(), 0.0);
      fill(phi_max_Other_Soliton_total.begin(), phi_max_Other_Soliton_total.end(), 0.0);
      fill(Alpha_Schamel_Other_Soliton_total.begin(), Alpha_Schamel_Other_Soliton_total.end(), 0.0);

      if (rank==0) cout << "-----------------------------------------------------"<<endl;
    
      //--- loading the JSON file ------------------------------------------
      file_json = json_load_file("./Potential_imposed.json", flags, &error);

      if(!file_json){
	cout << "in line " << error.line << ": " << error.text <<flush<< endl;
	exit(10);      
      }
      //--------------------------------------------------------------------
    
      //--- returning the potential object ------------------------------------  
      json_t* potential_array_json = json_object_get(file_json, "PHI");
      if(!potential_array_json){
	cerr << "in line " << error.line << ": " << error.text <<flush<<endl;
	exit(20);
      }
      size_t imposed_array_size =  json_array_size(potential_array_json);
      int Potential_array_size_int = static_cast<int>(imposed_array_size);
//        cout<<"size="<<Potential_array_size_int<<endl;
      i_element=0;
      for (int i_X = 0; i_X <= Potential_array_size_int-1; i_X++){
		Pair_X_Potential_json    = json_array_get(potential_array_json, i_element);
		if (!Pair_X_Potential_json) cout<<"can't catch element "<< i_element<<" of PHI" <<flush<<endl;

		X_json = json_array_get(Pair_X_Potential_json, 0);
		if (!X_json) cout<<"can't catch X of PHI" <<flush<<endl;
		double X_value = json_real_value(X_json);
//   		if (rank==0) cout << "X of Grid" <<" = " << X_value << endl;
      
		Potential_json = json_array_get(Pair_X_Potential_json, 1);
		if (!Potential_json) cout<<"can't catch Value of PHI" <<flush<<endl;
		double Potential_imposed_value = json_real_value(Potential_json);
//    		if (rank==0) cout << "PHI of Grid" <<" = " << Potential_imposed_value << endl;
		X_integer = int((0.0+X_value+0.01*P->dX)/P->dX);
//  		cout<<"{"<<X_value<<" "<< X_integer<<"  "<<flush;
		Potential_imposed_total[X_integer] = Potential_imposed_value;
        
                velocity_soliton = json_array_get(Pair_X_Potential_json, 2);
		if (!velocity_soliton) cout<<"can't catch Value of v_soliton" <<flush<<endl;
		double v_soliton_imposed_value = json_real_value(velocity_soliton);
		v_soliton_imposed_total[X_integer] = v_soliton_imposed_value;
        
                Beta_soliton = json_array_get(Pair_X_Potential_json, 3);
		if (!Beta_soliton) cout<<"can't catch Value of beta_soliton" <<flush<<endl;
		double Beta_soliton_imposed_value = json_real_value(Beta_soliton);
		Beta_imposed_total[X_integer] = Beta_soliton_imposed_value;
        
                Phi_maximum_soliton = json_array_get(Pair_X_Potential_json, 4);
 		if (!Phi_maximum_soliton) cout<<"can't catch Value of Phi_maximum_soliton" <<flush<<endl;
 		double Phi_maximum_soliton_imposed_value = json_real_value(Phi_maximum_soliton);
 		Phi_maximum_imposed_total[X_integer] = Phi_maximum_soliton_imposed_value;
        
                Alpha_Schamel_Main_Soliton_json = json_array_get(Pair_X_Potential_json, 5);
 		if (!Alpha_Schamel_Main_Soliton_json) cout<<"can't catch Value of Alpha_Schamel_Main_Soliton" <<flush<<endl;
 		double Alpha_Schamel_Main_Soliton_value = json_real_value(Alpha_Schamel_Main_Soliton_json);
 		Alpha_Schamel_Main_Soliton_total[X_integer] = Alpha_Schamel_Main_Soliton_value;        


                
                v_soliton_Other_Soliton_json = json_array_get(Pair_X_Potential_json, 6);
 		if (!v_soliton_Other_Soliton_json) cout<<"can't catch Value of v_soliton_Other_Soliton" <<flush<<endl;
 		double v_soliton_Other_Soliton_value = json_real_value(v_soliton_Other_Soliton_json);
 		v_soliton_Other_Soliton_total[X_integer] = v_soliton_Other_Soliton_value;        
		
                phi_max_Other_Soliton_json = json_array_get(Pair_X_Potential_json, 7);
 		if (!phi_max_Other_Soliton_json) cout<<"can't catch Value of phi_max_Other_Soliton" <<flush<<endl;
 		double phi_max_Other_Soliton_value = json_real_value(phi_max_Other_Soliton_json);
 		phi_max_Other_Soliton_total[X_integer] = phi_max_Other_Soliton_value;                 

                Alpha_Schamel_Other_Soliton_json = json_array_get(Pair_X_Potential_json, 8);
 		if (!Alpha_Schamel_Other_Soliton_json) cout<<"can't catch Value of Alpha_Schamel_Other_Soliton" <<flush<<endl;
 		double Alpha_Schamel_Other_Soliton_value = json_real_value(Alpha_Schamel_Other_Soliton_json);
 		Alpha_Schamel_Other_Soliton_total[X_integer] = Alpha_Schamel_Other_Soliton_value; 

                
                ++i_element;
	
      }
      		for (int iX = 0; iX <P->nX; iX++){
			Potential_imposed_array_total [iX] = Potential_imposed_total[iX];
		}
// 	  for (int i = 0; i <P->nX; i++){
//  	    if (Potential_imposed_total[i]==0.0){ cout<<"  {"<<i<<"  "<<Potential_imposed_total[i]<<"}";}
// 	  }
      

// 	if (rank == 0 ){
		for (int i_X = 0; i_X <= P->nX_cpu-1; i_X++){
                        Potential_imposed_array_cpu [i_X]   = Potential_imposed_total[i_X];
                        v_soliton_imposed_array_cpu [i_X]   = v_soliton_imposed_total[i_X];
                        Beta_soliton_imposed_array_cpu [i_X]= Beta_imposed_total[i_X];
                        Phi_maximum_soliton_imposed_array_cpu [i_X]= Phi_maximum_imposed_total[i_X];
                        Alpha_Schamel_Main_Soliton_array_cpu [i_X]= Alpha_Schamel_Main_Soliton_total[i_X];
                        v_soliton_Other_Soliton_array_cpu [i_X]= v_soliton_Other_Soliton_total[i_X];
                        phi_max_Other_Soliton_array_cpu [i_X]= phi_max_Other_Soliton_total[i_X];
                        Alpha_Schamel_Other_Soliton_array_cpu [i_X]= Alpha_Schamel_Other_Soliton_total[i_X];          
		}
		for (int i_rank = 1; i_rank <= size-1; i_rank++){
			int starting_nx = i_rank * P->nX_cpu;
			MPI_Send(&Potential_imposed_total[starting_nx] , P->nX_cpu, MPI_DOUBLE,i_rank,610, MPI_COMM_WORLD);
                        MPI_Send(&v_soliton_imposed_total[starting_nx] , P->nX_cpu, MPI_DOUBLE,i_rank,611, MPI_COMM_WORLD);
                        MPI_Send(&Beta_imposed_total[starting_nx] , P->nX_cpu, MPI_DOUBLE,i_rank,612, MPI_COMM_WORLD);
                        MPI_Send(&Phi_maximum_imposed_total[starting_nx] , P->nX_cpu, MPI_DOUBLE,i_rank,613, MPI_COMM_WORLD);
                        
                        MPI_Send(&Alpha_Schamel_Main_Soliton_total[starting_nx] , P->nX_cpu, MPI_DOUBLE,i_rank,614, MPI_COMM_WORLD);
                        MPI_Send(&v_soliton_Other_Soliton_total[starting_nx] , P->nX_cpu, MPI_DOUBLE,i_rank,615, MPI_COMM_WORLD);
                        MPI_Send(&phi_max_Other_Soliton_total[starting_nx] , P->nX_cpu, MPI_DOUBLE,i_rank,616, MPI_COMM_WORLD);
                        MPI_Send(&Alpha_Schamel_Other_Soliton_total[starting_nx] , P->nX_cpu, MPI_DOUBLE,i_rank,617, MPI_COMM_WORLD);
            
		}
		

		
	}else{
		MPI_Recv(&Potential_imposed_array_cpu[0], P->nX_cpu, MPI_DOUBLE,0 ,610, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&v_soliton_imposed_array_cpu[0], P->nX_cpu, MPI_DOUBLE,0 ,611, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&Beta_soliton_imposed_array_cpu[0], P->nX_cpu, MPI_DOUBLE,0 ,612, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&Phi_maximum_soliton_imposed_array_cpu[0], P->nX_cpu, MPI_DOUBLE,0 ,613, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&Alpha_Schamel_Main_Soliton_array_cpu[0], P->nX_cpu, MPI_DOUBLE,0 ,614, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&v_soliton_Other_Soliton_array_cpu[0], P->nX_cpu, MPI_DOUBLE,0 ,615, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&phi_max_Other_Soliton_array_cpu[0], P->nX_cpu, MPI_DOUBLE,0 ,616, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&Alpha_Schamel_Other_Soliton_array_cpu[0], P->nX_cpu, MPI_DOUBLE,0 ,617, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

		MPI_Bcast(&Potential_imposed_array_total[0], P->nX,  MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		
		
		
      	#ifdef debug_initial
		P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="$$$ Read_potential_imposed, finished $$$");
	#endif
  }
  //-----------------------------------------------------------------------------------------------
    
    
  //-----------------------------------------------------------------------------------------------
  void Maxwellian_Reflected_distribution_phi_imposed(double & X0, double& Vx0, double& Fp0, double& MPP0, double& Phi_imposed ) {        
		//---------------------------------------------------------------------------------
		double U0_Schamel, A_constant;
		U0_Schamel = 0.4; 
		//---------------------------------------------------------------------------------

		//---------------------------------------------------------------------------------
		double Energy_particle = ((P->norm_factor/2.0) * pow( (Vx0 -U0_Schamel) , 2)) //kinetic energy of a partilce in the soliton frame
					  + (P->Charge * Phi_imposed/P->Temp); //electric potential energy 
					  
		A_constant = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) ; 
		//---------------------------------------------------------------------------------

		
		//---------------------------------------------------------------------------------
		if       (Vx0 > 0.0 ){
			Fp0 = A_constant * exp(-pow( ( (sqrt(P->norm_factor/2.0) * U0_Schamel ) + sqrt(Energy_particle) ),2) );
			MPP0 = 1.0;
		
		  
		}else if (Vx0 < 0.0 ){
			Fp0 = A_constant * exp(-pow( ( (sqrt(P->norm_factor/2.0) * U0_Schamel ) - sqrt(Energy_particle) ),2) ); 
			MPP0 = 2.0;
		}
		//---------------------------------------------------------------------------------    
  }
  
  //-----------------------------------------------------------------------------------------------
  void initialization_DF_Kappa(double& gamma1, double& gamma2){
    double kappa1 = P->Kappa_DF;
	double kappa2 = P->Kappa_DF - 0.5;
	
	double Kmax = 100.0;
	int Kn = 100000;
	double d_kappa = Kmax/float(Kn);
    double k1, k2;
	for (int i = 0; i <= Kn-1; i++){
	    k1 = i*d_kappa;
	    k2 = (i+1) * d_kappa;
	    gamma1 = gamma1 + ( ( ( pow(k1,(kappa1-1)) * exp(-k1))+( pow(k2,(kappa1-1)) * exp(-k2)) )/2 ) * d_kappa;
	    gamma2 = gamma2 + ( ( ( pow(k1,(kappa2-1)) * exp(-k1))+( pow(k2,(kappa2-1)) * exp(-k2)) )/2 ) * d_kappa;
	}
	
  }
	//-----------------------------------------------------------------------------------------------
	
  double distribution_function_energy(double& Energy_kinetic,double& Temperature,  double& gamma1, double& gamma2) {
      double Fp0;
      if (P->DF_Vx_MK == 0){
            Fp0 =  Maxwellian_distribution_energy(Energy_kinetic, Temperature);
        }else if (P->DF_Vx_MK == 1) {
            Fp0 =   Kappa_distribution_energy(Energy_kinetic, Temperature, gamma1, gamma2);
        }    
        
        return Fp0;
  }
	
  //-----------------------------------------------------------------------------------------------
  double Maxwellian_distribution_energy(double& Energy_kinetic,double& Temperature){
		return exp(-Energy_kinetic/Temperature);
  }
      
  //-----------------------------------------------------------------------------------------------
  
  //-----------------------------------------------------------------------------------------------
  double Kappa_distribution_energy(double& Energy_kinetic,double& Temperature,  double& gamma1, double& gamma2) {
      return  sqrt(1.0/(P->Kappa_DF-1.5)) * 
	     (gamma1/gamma2)  * 
	     pow( (1.0 + (  (Energy_kinetic/Temperature) *  (1.0/(P->Kappa_DF-1.5)) )  ),(-P->Kappa_DF));  
  }
  //-----------------------------------------------------------------------------------------------
  
  
  //-----------------------------------------------------------------------------------------------
  double sign_function(double& input){
		if (input > 0) return 1.0;
		if (input < 0) return -1.0;
		return 0.0;
  }
  double Vx_phi_function(double& Energy_potential){
        return sqrt(abs(2.0 * Energy_potential / P->Mass));	
  }
  double Energy_potential_function(double& phi_input){
        return P->Charge * phi_input;
  }
  double Energy_kinetic_function(double& v_input){
		return(P->Mass/2.0) * pow(v_input,2);
  } 
 
 
  void Energy_kinetic_shifted_function(double& Vx_input, double& Phi_input, double& v_soliton, double& Energy_kinetic_shifted, double& Energy_kinetic_shifted_SF){
		double Energy_potential          = P->Charge * Phi_input;
		double Vx_input_SF               = Vx_input-v_soliton;
		double Energy_kinetic_SF         = Energy_kinetic_function(Vx_input_SF);
		       Energy_kinetic_shifted_SF = (Energy_kinetic_SF)  + (Energy_potential); 
		double Vx_input_shifted_SF       = sign_function(Vx_input_SF) * sqrt(2.0*Energy_kinetic_shifted_SF /P->Mass);
		double Vx_input_shifted          = Vx_input_shifted_SF  + v_soliton;
               Energy_kinetic_shifted    = Energy_kinetic_function(Vx_input_shifted);
  }
  
  double Fp_base_reflected_function(double& v_soliton, double& Energy_potential_total,double& gamma1, double& gamma2, double& A_constant ){
				double Energy_kinetic_soliton_reflected , Energy_kinetic_soliton_reflected_SF, Phi_zero ;
                double Energy_kinetic_soliton_reflected_upper , Energy_kinetic_soliton_reflected_upper_SF, Energy_kinetic_soliton_reflected_base ;
                double Vx_phi_total = Vx_phi_function(Energy_potential_total) + v_soliton;
                double Vx_phi_total_minus = -Vx_phi_function(Energy_potential_total) + v_soliton;
                double Fp0_base_reflected, Fp0; 
                
                Energy_kinetic_soliton_reflected = 0.0, Energy_kinetic_soliton_reflected_SF = 0.0;
                Energy_kinetic_soliton_reflected_upper = 0.0, Energy_kinetic_soliton_reflected_upper_SF = 0.0;
                Phi_zero = 0.0;
                Energy_kinetic_shifted_function(Vx_phi_total, Phi_zero, v_soliton,
                        Energy_kinetic_soliton_reflected, Energy_kinetic_soliton_reflected_SF);
                
                Energy_kinetic_shifted_function(Vx_phi_total_minus, Phi_zero, v_soliton,
                        Energy_kinetic_soliton_reflected_upper, Energy_kinetic_soliton_reflected_upper_SF);
                
				Fp0 = A_constant * distribution_function_energy(Energy_kinetic_soliton_reflected,P->Temp,gamma1,gamma2);
				Fp0 = Fp0 + A_constant * distribution_function_energy(Energy_kinetic_soliton_reflected_upper,P->Temp,gamma1,gamma2);
				Fp0_base_reflected = Fp0 / 2.0;
                Energy_kinetic_soliton_reflected_base = (Energy_kinetic_soliton_reflected+Energy_kinetic_soliton_reflected_upper)/2.0;
                return Fp0_base_reflected,Energy_kinetic_soliton_reflected_base ;
  }
  
  double Fp_reflected_side(double& Vx0, double& Fp0_base_reflected, double& Vx_phi_reflected, double& phi_imposed, double& Phi_maximum_domain, double& v_soliton, double& Alpha_Schamel, double& Vx_phi_reflected_side, double& A_constant, double& gamma1, double& gamma2){
        double v_1_negative, Energy_kinetic_shifted_reflected_side_1, Energy_kinetic_shifted_reflected_SF_side_1;
        double f_1_negative, v_1_positive, Energy_kinetic_shifted_reflected_side_2, Energy_kinetic_shifted_reflected_SF_side_2;
        double v_2_negative, f_1_positive, v_2_positive, Fp0, Vx0_SF, phi_imposed_reflected;
        double f_2_negative, f_2_positive, Alpha_mulitplication;
        
    v_1_negative = v_soliton - Vx_phi_reflected;
    phi_imposed_reflected = phi_imposed-Phi_maximum_domain;
    Energy_kinetic_shifted_function(v_1_negative,phi_imposed_reflected, v_soliton, Energy_kinetic_shifted_reflected_side_1, Energy_kinetic_shifted_reflected_SF_side_1 );
    Alpha_mulitplication = Alpha_Schamel * Energy_kinetic_shifted_reflected_SF_side_1;
	f_1_negative = Fp0_base_reflected * Maxwellian_distribution_energy(Alpha_mulitplication,P->Temp);
	Vx0_SF = Vx0-v_soliton;
	v_2_negative = v_soliton - Vx_phi_reflected_side;
	
	v_1_positive = v_soliton + Vx_phi_reflected;
    Energy_kinetic_shifted_function(v_1_positive,phi_imposed_reflected, v_soliton, Energy_kinetic_shifted_reflected_side_1, Energy_kinetic_shifted_reflected_SF_side_1 );
    
    Alpha_mulitplication = Alpha_Schamel * Energy_kinetic_shifted_reflected_SF_side_1;
	f_1_positive = Fp0_base_reflected * Maxwellian_distribution_energy(Alpha_mulitplication,P->Temp);
	
	Energy_kinetic_shifted_function(v_2_negative, phi_imposed, v_soliton, Energy_kinetic_shifted_reflected_side_2, Energy_kinetic_shifted_reflected_SF_side_2 );
	
	f_2_negative = A_constant * distribution_function_energy(Energy_kinetic_shifted_reflected_side_2,P->Temp,gamma1,gamma2);
	
	v_2_positive = v_soliton + Vx_phi_reflected_side;
	Energy_kinetic_shifted_function(v_2_positive, phi_imposed, v_soliton, Energy_kinetic_shifted_reflected_side_2, Energy_kinetic_shifted_reflected_SF_side_2 );
	
	
	f_2_positive = A_constant * distribution_function_energy(Energy_kinetic_shifted_reflected_side_2,P->Temp,gamma1,gamma2);
    
    if (Vx0_SF<0.0){
			Fp0 = ( (f_2_negative-f_1_negative) / (v_2_negative-v_1_negative) ) * (Vx0 - v_1_negative ) + f_1_negative;
    }
    else{
			Fp0 = ( (f_2_positive-f_1_positive) / (v_2_positive-v_1_positive) ) * (Vx0 - v_1_positive ) + f_1_positive;
    }
    return Fp0;
  }
  
  /*
  double assign_MPP0 (double & Energy_kinetic_shifted,double & Energy_potential_maximum,double & Energy_kinetic_shifted_SF,
                      double & Phi_maximum_imposed,double & Velocity_imposed,double & v_soliton_Other_Soliton,double & phi_max_Other_Soliton){
      double MPP0;  
      if (Energy_potential_maximum<0.0){
          if (Energy_kinetic_shifted_SF<-0.01){ //Free
              MPP0 = -1.0;
          }else if (Energy_kinetic_shifted_SF>0.01){ //Trapped
              MPP0 = -10.0;
          }else{ //Boundary between them
              MPP0 = -5.0;
          }
          
      }else{
           if (Energy_kinetic_shifted_SF<-Energy_potential_maximum-0.01){ //Free
              MPP0 = 7.0;
          }else if (Energy_kinetic_shifted_SF>-Energy_potential_maximum+0.01){ //Reflected?
              MPP0 = 1.0;
          }else{ //Boundary between them
              MPP0 = 5.0;
          }
      }

       MPP0 =    Energy_kinetic_shifted;  
    return MPP0;
  }
  
  double assign_MPP0_old (double & maximum_value, double & value, double & value_SF){
      double MPP0, delta_value;
      double base_mpp ;
      double number_of_mpp = 15.0;      
      double maximum_mpp = 10.0;
      delta_value = abs(maximum_value)/number_of_mpp;
      
      
      if (maximum_value>0.0){ // reflected particles
          base_mpp = 2.0;             
      }else{ // trapped particles
          base_mpp = 0.0; maximum_mpp = 8.0;
      }
      
          
      MPP0 = (int( value_SF/delta_value  ) * (( maximum_mpp-base_mpp)/number_of_mpp)  )+ base_mpp;
      base_mpp = -100;
      if (maximum_value<0.0){            
            delta_value = abs(delta_value);
//             if (value< (maximum_value+delta_value)){value =value;}
            
            if (value_SF<0.0){if (int( abs(value_SF/delta_value) )%2 != 0){ MPP0 = base_mpp;}} // trapped electrons
            else if (0.0<value_SF and value_SF<delta_value){MPP0 = base_mpp ;} //border between free and trapped electrons 
            else if (delta_value<value) { if (int( sqrt(value/delta_value) )%3 != 0){ MPP0 = base_mpp;}} // free electrons

      }
      
      if (maximum_value>0.0){
            if (0.0<value and value<delta_value){value = value;} // inner circle of reflected ions
            else if (delta_value<value_SF and value_SF<=maximum_value) {if (int( (value_SF/delta_value) )%2 != 0){ MPP0 = base_mpp;}} // reflected ions
            else if (maximum_value<value_SF and value_SF<=maximum_value+(delta_value)) { MPP0 = base_mpp;} // border between reflected and free ions
            else if (value > (maximum_value+delta_value) ) { if (int( sqrt(value/delta_value) )%3 != 0){ MPP0 = base_mpp;}}  // free ions
      }
      
      
      
      
      
//       if (int( (value/maximum_value) *10.0 )%2 == 1){ MPP0 = 10.0;      }
      return MPP0;
    }
    */
  
  
  
	#ifdef multiple_beta_for_trapped_DF
	//==============================================================================
	double searching_distribution_function_of_trapped_population(double & Energy_kinetic_shifted_SF, vector<double> & Energy_Kinentic_Jenab_DF, vector<double> &  Beta_Jenab_DF, double &Temperature){
		double beta_minimum,beta, DF;
		int i_search;
		bool flag_search;
		flag_search = true;
		i_search = 0;
		DF = -1.0;
		do {
			if (i_search >= Beta_Jenab_DF.size()){
				flag_search = false;
			}
			if (abs(Energy_kinetic_shifted_SF) > abs(Energy_Kinentic_Jenab_DF[i_search-1]) 
				and abs(Energy_kinetic_shifted_SF) < abs(Energy_Kinentic_Jenab_DF[i_search])){
				beta = Beta_Jenab_DF[i_search]* Energy_kinetic_shifted_SF;
				DF = Maxwellian_distribution_energy( beta ,Temperature);
				flag_search = false;
			}
			i_search = i_search + 1;
		}while(flag_search == true );
		
		if (DF == -1.0){
			beta_minimum = Beta_Jenab_DF[Beta_Jenab_DF.size()-1]* Energy_kinetic_shifted_SF;
			DF = Maxwellian_distribution_energy(beta_minimum ,Temperature);
		}
		
// 		if (abs(Energy_kinetic_shifted_SF) > max_element(Energy_Kinentic_Jenab_DF) ){
// 				beta_minimum = Beta_Jenab_DF[Beta_Jenab_DF.size()-1]* Energy_kinetic_shifted_SF;
// 				DF = Maxwellian_distribution_energy(beta_minimum ,Temperature);
// 		}
// 		if (abs(Energy_kinetic_shifted_SF) < min_element(Energy_Kinentic_Jenab_DF)){
// 				beta = 0.0 * Energy_kinetic_shifted_SF;
// 				DF = Maxwellian_distribution_energy(beta ,Temperature);
// 				}
		return DF;
		
		}
	#endif	
	//==============================================================================
  //-----------------------------------------------------------------------------------------------
  void Schamel_distribution_phi_imposed(double & X0, double& Vx0, double& Fp0, double& EPP0, double& MPP0, 
                                        double& Phi_imposed, double& Phi_imposed_CoMoving_Correction, double& Beta_Schamel, double& Phi_maximum_imposed, 
                                        double& Velocity_imposed, double&  Alpha_Schamel_Main_Soliton, 
                                        double&  v_soliton_Other_Soliton, double&  phi_max_Other_Soliton, 
                                        double&  Alpha_Schamel_Other_Soliton, double& gamma1, double& gamma2,
				       vector<double> & Energy_Kinentic_Jenab_DF, vector<double> & Beta_Jenab_DF) { 

        //-----------------------------------------------------------------------------------------
        double v_soliton, A_constant;
        double Vx0_SF, Energy_kinetic_shifted_SF, Energy_kinetic_shifted, Energy_kinetic_soliton, Energy_potential,Vx_phi, Energy_potential_maximum;
	double Energy_kinetic_shifted_CoMoving_Correction, Energy_kinetic_shifted_SF_CoMoving_Correction;
        
        double phi_imposed_reflected, Energy_potential_reflected, Vx_phi_reflected;
        double Energy_kinetic_shifted_reflected_SF, Energy_kinetic_shifted_reflected;
        double phi_imposed_reflected_2, Fp0_base_reflected, Energy_potential_maximum_side, Energy_kinetic_soliton_reflected_base;
        
        double Phi_maximum_domain_side,phi_imposed_reflected_side,  Vx_phi_reflected_side;
        //-----------------------------------------------------------------------------------------
        
        
        //-----------------------------------------------------------------------------------------
        double Vx_phi_reflected_side_other_soliton,phi_max_side_Other_Soliton, Energy_potential_max_side_Other_Soliton;
        double Energy_kinetic_shifted_other_soliton, Energy_kinetic_shifted_other_soliton_SF, v_other_soliton_SF, v_other_soliton_in_SF;
        double Energy_kinetic_shifted_reflected_other_soliton, Energy_kinetic_shifted_reflected_other_soliton_SF;
        double  phi_imposed_reflected_2_Other_Soliton;
        double Fp0_base_reflected_Other_Soliton, Energy_kinetic_soliton_reflected_base_OS, Energy_potential_maximum_OS;
        double minus_phi_imposed;
        //-----------------------------------------------------------------------------------------
        
        
        //-----------------------------------------------------------------------------------------
        A_constant = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor);

        phi_imposed_reflected_2 = Phi_imposed-Phi_maximum_imposed;
        phi_imposed_reflected_2_Other_Soliton = 0.0-phi_max_Other_Soliton;
        
        Fp0 = 0.0;
	double Fp0_trapped;
        
        
        v_soliton = Velocity_imposed;

        Energy_potential_maximum = Energy_potential_function(Phi_maximum_imposed) ;
        Energy_potential_maximum_OS = Energy_potential_function(phi_max_Other_Soliton) ;
        
        
        Energy_potential         = Energy_potential_function(Phi_imposed); 
        Vx_phi                   = Vx_phi_function  (Energy_potential);
        
        //-----------------------------------------------------------------------------------------
        Phi_maximum_domain_side = Phi_maximum_imposed * 1.0;//1.1;
        phi_imposed_reflected_side = Phi_maximum_domain_side - Phi_imposed;
        if (abs(Phi_imposed) >  abs(Phi_maximum_domain_side)){
                phi_imposed_reflected_side = 0.0;
        }
        Energy_potential_maximum_side = Energy_potential_function( phi_imposed_reflected_side);
        Vx_phi_reflected_side  = Vx_phi_function(Energy_potential_maximum_side);
        //----------------------------------------------------------------------------------------
        
        //-----------------------------------------------------------------------------------------
        Vx0_SF = Vx0-v_soliton;
        Fp0_base_reflected, Energy_kinetic_soliton_reflected_base = Fp_base_reflected_function(v_soliton,Energy_potential_maximum,gamma1,gamma2,A_constant);
        Fp0_base_reflected_Other_Soliton, Energy_kinetic_soliton_reflected_base_OS = Fp_base_reflected_function(v_soliton_Other_Soliton,Energy_potential_maximum_OS,gamma1,gamma2,A_constant);
   
        Energy_kinetic_shifted_function(Vx0,Phi_imposed, v_soliton, Energy_kinetic_shifted, Energy_kinetic_shifted_SF);
	
	Energy_kinetic_shifted_function(Vx0,Phi_imposed_CoMoving_Correction, v_soliton, Energy_kinetic_shifted_CoMoving_Correction, Energy_kinetic_shifted_SF_CoMoving_Correction);

        
        phi_imposed_reflected = Phi_maximum_imposed - Phi_imposed;
        
        if (abs(Phi_imposed)>abs(Phi_maximum_imposed)){
            phi_imposed_reflected = 0.0;
        }
        Energy_potential_reflected = Energy_potential_function(phi_imposed_reflected); 
        Vx_phi_reflected = Vx_phi_function(Energy_potential_reflected); 
    
        Energy_kinetic_soliton = Energy_kinetic_function(v_soliton);
        
        
        Energy_kinetic_shifted_function(Vx0,phi_imposed_reflected_2, v_soliton, Energy_kinetic_shifted_reflected, Energy_kinetic_shifted_reflected_SF);        
        Energy_kinetic_shifted_function(Vx0,phi_imposed_reflected_2_Other_Soliton, v_soliton_Other_Soliton, Energy_kinetic_shifted_reflected_other_soliton, Energy_kinetic_shifted_reflected_other_soliton_SF);
        
        //-----------------------------------------------------------------------------------------
        phi_max_side_Other_Soliton = phi_max_Other_Soliton * 1.0;//1.1;
        Energy_potential_max_side_Other_Soliton = Energy_potential_function(phi_max_side_Other_Soliton );
	Vx_phi_reflected_side_other_soliton  = Vx_phi_function(Energy_potential_max_side_Other_Soliton);
        
        minus_phi_imposed = -Phi_imposed;
        Energy_kinetic_shifted_function(v_soliton_Other_Soliton,minus_phi_imposed, v_soliton, Energy_kinetic_shifted_other_soliton, Energy_kinetic_shifted_other_soliton_SF);
        v_other_soliton_in_SF = v_soliton_Other_Soliton-v_soliton;
        v_other_soliton_SF =  sign_function(v_other_soliton_in_SF ) * sqrt(2.0*Energy_kinetic_shifted_other_soliton_SF  /P->Mass);
        //-----------------------------------------------------------------------------------------
        
        
        EPP0 = 0.0;
        MPP0 = 0.0;
        
        double Maxwellian_Parameter_reflected_other_soliton_SF = Energy_kinetic_shifted_reflected_other_soliton_SF * Alpha_Schamel_Other_Soliton;
        double Maxwellian_Parameter_reflected_SF = Energy_kinetic_shifted_reflected_SF * Alpha_Schamel_Main_Soliton;
        double Maxwellian_Parameter_shifted_SF = Energy_kinetic_shifted_SF_CoMoving_Correction * Beta_Schamel;//Energy_kinetic_shifted_SF * Beta_Schamel;
        //=========================================================================================
        if (Energy_potential_maximum > 0.0){            // Sr) Species partitioning, reflected species
		#ifdef flag_reflected_distribution_full_mixture_mehtod
		if ( abs(Vx0_SF) <  Vx_phi_reflected_side 
                && Vx0_SF>= v_other_soliton_SF - abs(Vx_phi_reflected_side_other_soliton) 
                && Vx0_SF<= v_other_soliton_SF + abs(Vx_phi_reflected_side_other_soliton) ){ 
                                        Fp0 = (Fp0_base_reflected_Other_Soliton + Fp0_base_reflected)/2.0;//overlapping
                                        EPP0 = ( (Energy_kinetic_soliton_reflected_base_OS + Maxwellian_Parameter_reflected_other_soliton_SF) / P->Temp 
                                             + (Energy_kinetic_soliton_reflected_base + Maxwellian_Parameter_reflected_SF) / P->Temp ) / 2.0; 
                                        MPP0 = 40.0;
            } else if (abs(Vx0_SF) <  Vx_phi_reflected_side){ //main soliton
                                        Fp0 = Fp0_base_reflected * Maxwellian_distribution_energy( Maxwellian_Parameter_reflected_SF,P->Temp);
                                        EPP0 = (Energy_kinetic_soliton_reflected_base+ Maxwellian_Parameter_reflected_SF)/P->Temp;
                                        if (X0>200.0){MPP0=20.0;}else{MPP0=30.0;}
            } else if (Vx0_SF>= v_other_soliton_SF - abs(Vx_phi_reflected_side_other_soliton) 
                    && Vx0_SF<= v_other_soliton_SF + abs(Vx_phi_reflected_side_other_soliton)){ // other soliton
                                        Fp0 = Fp0_base_reflected_Other_Soliton * Maxwellian_distribution_energy (Maxwellian_Parameter_reflected_other_soliton_SF,P->Temp);
                                        EPP0 = (Energy_kinetic_soliton_reflected_base_OS + Maxwellian_Parameter_reflected_other_soliton_SF) / P->Temp;
                                        if (X0>200.0){MPP0=20.0;}else{MPP0=30.0;}
            } else {
                                        Fp0 = A_constant * distribution_function_energy(Energy_kinetic_shifted,P->Temp,gamma1,gamma2);
                                        EPP0 = Energy_kinetic_shifted / P->Temp;
                                        MPP0 = 10.0;
            }
		#endif
		#ifdef flag_reflected_distribution_simple_cut
					Fp0 = A_constant * distribution_function_energy(Energy_kinetic_shifted,P->Temp,gamma1,gamma2);
                                        EPP0 = Energy_kinetic_shifted / P->Temp;
                                        MPP0 = 10.0;
		#endif
                        /*if ( abs(Vx0_SF) >=  Vx_phi_reflected_side ){ //Sr-Vf) V-direction partitioning, bended population
                            if (Vx0_SF>= v_other_soliton_SF - abs(Vx_phi_reflected_side_other_soliton) 
                                && Vx0_SF<= v_other_soliton_SF + abs(Vx_phi_reflected_side_other_soliton) ){ //Sr-Vf-Vfr) V-direction partitioning, bended population of Free-Reflected particles for other soliton
                                     Fp0 = Fp0_base_reflected
                                      * Maxwellian_distribution_energy
                                      (Maxwellian_Parameter_reflected_other_soliton_SF,P->Temp)        ;
                                      EPP0 = (Energy_kinetic_soliton_reflected_base + Maxwellian_Parameter_reflected_other_soliton_SF) / P->Temp;
                                      MPP0 = 30.0;//30*number_of_zero
                            } 
                            else {
                                Fp0 = A_constant * distribution_function_energy(Energy_kinetic_shifted,P->Temp,gamma1,gamma2);
                                EPP0 = Energy_kinetic_shifted / P->Temp;
                                MPP0 = 20.0;//*number_of_zeros;
                            }
                        }
                        else { //Sr-Vs) V-direction partitioning, reflected population without side correction
                            Fp0 = Fp0_base_reflected * Maxwellian_distribution_energy( Maxwellian_Parameter_reflected_SF,P->Temp);
                            EPP0 = (Energy_kinetic_soliton_reflected_base+ Maxwellian_Parameter_reflected_SF)/P->Temp;
                            MPP0 = 10.0;//*number_of_zeros;
                         }*/
// jenab ... add the side correction part!
        }//=========================================================================================
        else { // St) Species partitioning, trapped species
                    if ( abs(Vx0_SF) >=  Vx_phi ){ // St-Vf) V-direction partitioning, free population
                                    Fp0 = A_constant * distribution_function_energy(Energy_kinetic_shifted,P->Temp,gamma1,gamma2);
                                    EPP0 = Energy_kinetic_shifted/P->Temp;
                                    MPP0 = 10.0;//*number_of_zeros);
                    }
                    else{ // St-Vt) V-direction partitioning, trapped population
                        Fp0 = A_constant 
                            * distribution_function_energy(Energy_kinetic_soliton,P->Temp,gamma1,gamma2) ;
			#ifdef single_beta_for_trapped_DF
				Fp0_trapped = Maxwellian_distribution_energy(Maxwellian_Parameter_shifted_SF,P->Temp);
			#endif
			#ifdef multiple_beta_for_trapped_DF
				Fp0_trapped = searching_distribution_function_of_trapped_population (Energy_kinetic_shifted_SF_CoMoving_Correction, Energy_Kinentic_Jenab_DF, Beta_Jenab_DF, P->Temp);
			#endif
			Fp0 = Fp0 * Fp0_trapped;
			if (Beta_Schamel ==0.0){
				EPP0 = (Energy_kinetic_soliton - Energy_kinetic_shifted_SF_CoMoving_Correction)/P->Temp;
			}else{
				EPP0 = (Energy_kinetic_soliton + Maxwellian_Parameter_shifted_SF)/P->Temp;
			}
			if (X0>P->iX_switching_EPP){ EPP0 = EPP0 + 2000.0; }else{ EPP0 = EPP0 + 3000.0;}
// 			  if (X0>1000.0){ MPP0 =  30.0; }else{ MPP0  = 20.0;}
//                             if      (X0<750.0)            { MPP0  = 20.0; }
// 			    else if (X0>750.0 && X0<1750) { MPP0 =  30.0; }
// 			    else if (X0>1750.0 && X0<2750){ MPP0 =  40.0; }
// 			    else if (X0>2750.0)           { MPP0 =  50.0; }
			    
			    //---- tracker for tracing speed of the solitons ----- 
// 			    if ( abs(Energy_kinetic_shifted_SF) > 0.99*abs(Energy_potential_maximum)  ){
// 				    if (X0>1000.0){ MPP0 =  2000.0; }else{ MPP0  = 1000.0;}
// 					if      (X0<750.0)            { MPP0  = 1000.0; }
// 					else if (X0>750.0 && X0<1750) { MPP0 =  2000.0; }
// 					else if (X0>1750.0 && X0<2750){ MPP0 =  3000.0; }
// 					else if (X0>2750.0)           { MPP0 =  4000.0; }
// 			    }
			    
                             
                    }
        }//=========================================================================================

        
        
        //=========================================================================================
//         MPP0 = assign_MPP0_old (Energy_potential_maximum, Energy_kinetic_shifted, Energy_kinetic_shifted_SF);
//         MPP0 = assign_MPP0 (Energy_kinetic_shifted, Energy_potential_maximum,Energy_kinetic_shifted_SF, 
//                             Phi_maximum_imposed, Velocity_imposed, v_soliton_Other_Soliton, phi_max_Other_Soliton);
//         MPP0 = Energy_kinetic_shifted_SF;
        //=========================================================================================
	}
  //-----------------------------------------------------------------------------------------------
  
  
	/*
	//***********************************************************************************************
	void DF_INITIALIZATION () {
	double X0,Vx0,Fp0;	
	double random_number, random_number_F;
	unsigned int seed_number;
	bool Schamel_DF_Flag = false;
	vector<double> * DF_one_Vx;
	DF_one_Vx  = new vector<double>(P->nVx);
	//     DF_one_Vx->push_back((*iter).Data_PP->at(P->x));
	fill(DF_one_Vx->begin(), DF_one_Vx->end(), 0.0);
	
	double Velocity_Border_Soliton = 0.0;
	
	double Kmax, d_kappa, gamma1, gamma2, kappa1, kappa2, k1,k2;
	int Kn;
	//       int time_name = 9900;
	
	// feeding and reading parameters -------------------------------
	vector<int> Num_PP_X(P->nX_cpu);
	vector<int> Integer_PP_X_cpu(P->nX_cpu);
	
	vector<double>  PP_X(1),  PP_Vx(1), PP_Fp(1) ;
	PP_X[0] = 0.0; PP_Vx[0]=0.0; PP_Fp[0]=0.0;
	
	//---------------------------------------------------------------

	//       Read_Num_PhasePoints (Num_PP_X, time_name);
	MPI_Barrier(MPI_COMM_WORLD);
	
	//       Read_X_Vx_Fp (PP_X, PP_Vx, PP_Fp, Integer_PP_X_cpu, Num_PP_X, time_name);
	MPI_Barrier(MPI_COMM_WORLD);
	
		
	seed_number = (unsigned)time(NULL)*(rank+1); // in order to have different sequece of pesudo-random number for each cpu, rank must be inside the formula for seeding
	srand(seed_number); //careful not to use srand() twice in the code unless you are sure that they are apart from each other in time at least for one second
	
		//---------------------------------------------------------------------------------

	//==================================================================================
	vector<double> Num_Density_imposed_array_cpu;
	Num_Density_imposed_array_cpu.resize (P->nX_cpu);      
	fill(Num_Density_imposed_array_cpu.begin(), Num_Density_imposed_array_cpu.end(), 1.0);      
	Read_Phi_imposed (Num_Density_imposed_array_cpu);      
	//==================================================================================
	
		
		kappa1 = P->Kappa_DF;
		kappa2 = P->Kappa_DF - 0.5;
		
		Kmax = 100.0;
		Kn = 100000;
		d_kappa = Kmax/Kn;
	
		for (int i = 0; i <= Kn-1; i++){
		k1 = i*d_kappa;
		k2 = (i+1) * d_kappa;
		gamma1 = gamma1 + ( ( ( pow(k1,(kappa1-1)) * exp(-k1))+( pow(k2,(kappa1-1)) * exp(-k2)) )/2 ) * d_kappa;
		gamma2 = gamma2 + ( ( ( pow(k1,(kappa2-1)) * exp(-k1))+( pow(k2,(kappa2-1)) * exp(-k2)) )/2 ) * d_kappa;
		}
		//---------------------------------------------------------------------------------
	
	fill(DF_g_1D->begin(), DF_g_1D->end(), 10.0);
	fill(MPP_g_1D->begin(), MPP_g_1D->end(), 0.0);
	//---------------------------------------------------------------------------------
	for (int iX = 0; iX<=P->nX_cpu-1; iX++)  {
	// 	if (rank==2){cout<<", X="<<P->X_min+ iX*P->dX<<",  "<<flush;}
		// preparation for feeding/constructing phase points
		int trapped_flag = 0;
		//double X_total  = P->X_min  + double(iX)  * P->dX ;
		bool feeding_flag = false;
	// 	double X_feeding_max, X_feeding_min;
	// 	X_feeding_max = P->X_max; X_feeding_min = 0; // full restart condition 
	// 	X_feeding_max = 2000.0; X_feeding_min = 1000.0;
		
		// to enable feeding process, next three lines should be activated.
	// 	if (X_total<=X_feeding_max && X_total>=X_feeding_min){
	// 	  feeding_flag = true;
	// 	}
		
		if (feeding_flag == true){
		Injection (PP_X, PP_Vx, PP_Fp, iX, Num_PP_X, Integer_PP_X_cpu);	  
		}else{
		random_number_F =  ((double)(rand() % 100)/100.0);
		if (rand() % 2==0){
		random_number_F =  random_number_F*(-1);//((double) rand() / (double)(RAND_MAX)) ;// / (double)(RAND_MAX)) ;
		} 
		
		//==================================================================================
		double Num_Density_imposed_element = Num_Density_imposed_array_cpu[iX];
		double U0_Schamel = 8.8;//9.7526;//10.0; 
		Velocity_Border_Soliton = U0_Schamel;
		double MPP0 = 0.0;
		
		
	// 	//------------------------------------
	//  	double A_IGP = 2.5;
	//  	double delta_IGP = 25.0;
	//  	double X_IGP = 51.2;
	//  	X0  = P->X_min  + double(iX)  * P->dX;
	//  	Velocity_Border_Soliton =  U0_Schamel + A_IGP * exp( -pow ( (X0  - X_IGP)/delta_IGP ,2) ) ;
	// 	Creat_One_Line_DF(DF_one_Vx, Velocity_Border_Soliton,U0_Schamel);
	// 	//------------------------------------
		
		
	//   		if (abs(Num_Density_imposed_element-1.0)>0.01 && iX>0){
	// 		  Velocity_Border_Soliton = Velocity_Border_Soliton - 2.0*P->dVx;
	// 		}else{
	// 		  Velocity_Border_Soliton = U0_Schamel;
	// 		}
		if (Num_Density_imposed_element > 1.0){
			Initial_Step_finding_Soliton_Border_more_than_one(DF_one_Vx,iX, Num_Density_imposed_element, Velocity_Border_Soliton,U0_Schamel);
		}
		if (Num_Density_imposed_element < 1.0){
			Initial_Step_finding_Soliton_Border_less_than_one(DF_one_Vx,iX, Num_Density_imposed_element, Velocity_Border_Soliton,U0_Schamel);
		}
		
		for (int iVx = 0; iVx<=P->nVx-1; iVx++) {
			DF_g_1D->at(iX * (P->nVx) + iVx) = DF_one_Vx->at(iVx);
		}
	// 	if (rank == 0){cout<<P->Name<<"*************  X="<<iX*P->dX+P->X_min<<" out-->"<<"V_border="<<Velocity_Border_Soliton-U0_Schamel <<"  Num_Den="<< Num_Density_imposed_element<<endl<<flush;}
		//==================================================================================
		
	// 	cout<<rank<<" finished the initial step iX="<<iX<<endl<<flush;
		for (int iVx = 0; iVx<=P->nVx-1; iVx++) {
		for (int i = 1; i<=P->pX; i++) 	{
		for (int j = 1; j<=P->pVx; j++) 	  {

			random_number = ((double) rand() / (double)(RAND_MAX)) ;
			Vx0 = P->Vx_min + double(iVx) * P->dVx +random_number * P->dVx; //j*(dVx_p_s/4.0)//
			
			random_number = ((double) rand() / (double)(RAND_MAX)) ;
			X0  = P->X_min  + double(iX)  * P->dX  + random_number * P->dX;//i*(dX_p_s/4.0)//
			if(P->DF_Vx_MK==0){
			Maxwellian_Jet_distribution(Vx0, Fp0);// Vx1_J, Temp1_J, Den1_J)
			}else if(P->DF_Vx_MK==1){
			Kappa_Distribution(Vx0, Fp0, gamma1,gamma2);
			}else if (P->DF_Vx_MK==2) {
			//Maxwellian_Schamel_distribution(X0,Vx0,Fp0,PHI_Sagdeev,trapped_flag);
			Maxwellian_Schamel_distribution_Num_Den_Based(X0,Vx0,Fp0,Velocity_Border_Soliton,U0_Schamel,trapped_flag);
			Schamel_DF_Flag = true;
			}
			
			if (Schamel_DF_Flag == false){
			if(P->DF_X_NMRE==1){
				DF_add_modes(X0,Fp0);
			}else if(P->DF_X_NMRE==2){
				DF_add_randomness(Fp0,random_number_F);
			}else if(P->DF_X_NMRE==3){
				DF_add_Envelope_Soliton(X0,Fp0);
	// 			DF_add_3Modes(X0,Fp0);
			}
			}
			
			//IF (Fp0/= 0.0){ // adding the phase PhasePoints  
			if (Fp0 >= 0.0){
			PhasePointClass PhasePoint (X0,Vx0,Fp0,MPP0);
			//PhasePoint= new PhasePointClass(X0,Vx0,Fp0);
			PhasePoints->push_back(PhasePoint);
			// delete PhasePoint;
			}
		}
		}
		}
		} // end of if for feeding or constructing
		} // end of for on iX
		//---------------------------------------------------------------------------------
		delete DF_one_Vx;
	//        exit(0);
	}
	//-----------------------------------------------------------------------------------------------
	void Initial_Step_finding_Soliton_Border_more_than_one(vector<double>* DF_one_Vx, int& iX, double& Num_Density_imposed_element, double& Velocity_Border_Soliton, double& U0_Schamel){
			
			fill(DF_one_Vx->begin(), DF_one_Vx->end(), 0.0);
			double Num_Den_Sub = 0.0;
			double Num_Den_Sub_0 = 0.0;
			//==== rought tunning ==========================================================================
			double delta_Velocity = 50.0*P->dVx;
			int NumberIteration_OuterLevel = 0;
			int NumberIteration_InnerLevel = 0;
			
			//------------------------------------------

			//------------------------------------------

			//------------------------------------------
	// 		if (rank == 7){cout<<"In -->"<<P->Name<<", X="<<iX*P->dX+P->X_min<<"V_border="<<Velocity_Border_Soliton-U0_Schamel<<"   "<< Num_Density_imposed_element<<endl<<flush;}
			do {  
				do {
	//      			if (rank == 7){cout<<P->Name<<" X="<<iX*P->dX+P->X_min<<" rough-->"<<"V_border="<<Velocity_Border_Soliton-U0_Schamel <<"  Num_Den="<<Num_Den_Sub <<"   "<< Num_Density_imposed_element<<endl<<flush;}
					Velocity_Border_Soliton = Velocity_Border_Soliton + delta_Velocity;
					Creat_One_Line_DF(DF_one_Vx, Velocity_Border_Soliton,U0_Schamel);
					Integrate_for_Num_Density(DF_one_Vx,Num_Den_Sub);
					if (Num_Den_Sub >  Num_Density_imposed_element){
						Velocity_Border_Soliton = Velocity_Border_Soliton - delta_Velocity;
						if (rank == 7){cout<<"$"<<flush;}
					} else{
					Num_Den_Sub_0 = Num_Den_Sub;
					if (rank == 7){cout<<"("<<NumberIteration_OuterLevel<<
						"%"<<((Num_Den_Sub_0-Num_Density_imposed_element)/Num_Den_Sub_0)*100.0<<")"<<flush;}
					
					}
					NumberIteration_InnerLevel += 1;
					if (NumberIteration_InnerLevel==1000){if (rank == 7){cout<<"***"<<flush;}}
				}while ( Num_Den_Sub <  Num_Density_imposed_element && NumberIteration_InnerLevel<1000);
				NumberIteration_OuterLevel += 1;
				delta_Velocity = delta_Velocity * 0.95;//P->dVx/(NumberIteration_OuterLevel);
				if (NumberIteration_OuterLevel==1000){if (rank == 7){cout<<"  Not Converged"<<flush;}}
			}while( (abs(Num_Den_Sub_0-Num_Density_imposed_element)/Num_Density_imposed_element) > (0.001/100.0) && NumberIteration_OuterLevel<1000 );
			
	// 		
	// 		do {
	// //       			if (rank == 7){cout<<P->Name<<"  X="<<iX*P->dX+P->X_min<<" Rough-->"<<"V_border="<<Velocity_Border_Soliton-U0_Schamel<<"  Num_Den="<<Num_Den_Sub <<"   "<< Num_Density_imposed_element<<endl<<flush;}
	//  			if (rank == 7){cout<<"*";}
	// 			Velocity_Border_Soliton = Velocity_Border_Soliton + delta_Velocity;
	//   			Creat_One_Line_DF(DF_one_Vx, Velocity_Border_Soliton,U0_Schamel);
	//   			Integrate_for_Num_Density(DF_one_Vx,Num_Den_Sub);
	//    			if (Num_Den_Sub >  Num_Density_imposed_element){Velocity_Border_Soliton = Velocity_Border_Soliton - delta_Velocity;if (rank == 7){cout<<"$";}}
	// 		}while ( Num_Den_Sub <  Num_Density_imposed_element);
	// 		//==============================================================================================
	// 
	// 		//==== fine tunning ============================================================================
	// 		delta_Velocity = 0.1*P->dVx;
	// 		do {   
	//  			if (rank == 7){cout<<P->Name<<"  X="<<iX*P->dX+P->X_min<<" Fine-->"<<"V_border="<<Velocity_Border_Soliton-U0_Schamel<<"  Num_Den="<<Num_Den_Sub <<"   "<< Num_Density_imposed_element<<endl<<flush;}
	// // 			if (rank == 7){cout<<"=";}
	// 			Velocity_Border_Soliton = Velocity_Border_Soliton + delta_Velocity;
	// 			Creat_One_Line_DF(DF_one_Vx, Velocity_Border_Soliton,U0_Schamel);
	// 			Integrate_for_Num_Density(DF_one_Vx,Num_Den_Sub);
	//    			if (Num_Den_Sub >  Num_Density_imposed_element){Velocity_Border_Soliton = Velocity_Border_Soliton - delta_Velocity;if (rank == 7){cout<<"$";}}
	// 		}while ( Num_Den_Sub <  Num_Density_imposed_element);
	// 		//==============================================================================================  
	// 
	// 		//==== fine tunning ============================================================================
	// 		delta_Velocity = 0.01*P->dVx;
	// 		do {
	//       			if (rank == 7){cout<<P->Name<<"  X="<<iX*P->dX+P->X_min<<" Super-Fine-->"<<"V_border="<<Velocity_Border_Soliton-U0_Schamel<<"  Num_Den="<<Num_Den_Sub <<"   "<< Num_Density_imposed_element<<endl<<flush;}
	// // 			if (rank == 7){cout<<"-";}
	// 			Velocity_Border_Soliton = Velocity_Border_Soliton + delta_Velocity;
	// 			Creat_One_Line_DF(DF_one_Vx, Velocity_Border_Soliton,U0_Schamel);
	// 			Integrate_for_Num_Density(DF_one_Vx,Num_Den_Sub);
	//    			if (Num_Den_Sub > Num_Density_imposed_element){Velocity_Border_Soliton = Velocity_Border_Soliton - delta_Velocity;if (rank == 7){cout<<"$";}}
	// 		}while ( Num_Den_Sub < Num_Density_imposed_element);
	// 		//============================================================================================== 
	// 
			if (rank == 7){cout<<endl<<"Out -->"<<P->Name<<", X="<<iX*P->dX+P->X_min<<"V_border="<<Velocity_Border_Soliton-U0_Schamel<<"  Num Den="<<Num_Den_Sub_0<<",  Imposed="<< Num_Density_imposed_element<<
			",  %"<<((Num_Den_Sub_0-Num_Density_imposed_element)/Num_Den_Sub_0)*100.0<<endl<<flush;}
	}
	//-----------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------
	void Initial_Step_finding_Soliton_Border_less_than_one(vector<double>* DF_one_Vx, int& iX, double& Num_Density_imposed_element, double& Velocity_Border_Soliton, double& U0_Schamel){
			
			fill(DF_one_Vx->begin(), DF_one_Vx->end(), 0.0);
			double Num_Den_Sub = 0.0;
			double Num_Den_Sub_0 = 0.0;
			//==== rought tunning ==========================================================================
			double delta_Velocity = 5.0*P->dVx;
			int NumberIteration_OuterLevel = 0;
			int NumberIteration_InnerLevel = 0;
		
			//------------------------------------------

			//------------------------------------------

			//------------------------------------------
			if (rank == 7){cout<<P->Name<<" >>>>>>>  X="<<iX*P->dX+P->X_min<<" in-->"<<"V_border="<<Velocity_Border_Soliton-U0_Schamel<<"   "<< Num_Density_imposed_element<<endl<<flush;}
			do {  
				do {
	//      			if (rank == 7){cout<<P->Name<<" X="<<iX*P->dX+P->X_min<<" rough-->"<<"V_border="<<Velocity_Border_Soliton-U0_Schamel <<"  Num_Den="<<Num_Den_Sub <<"   "<< Num_Density_imposed_element<<endl<<flush;}
					if (rank == 7){cout<<"-"<<NumberIteration_OuterLevel<<flush;}
					Num_Den_Sub_0 = Num_Den_Sub;
					Velocity_Border_Soliton = Velocity_Border_Soliton + delta_Velocity;
					Creat_One_Line_DF(DF_one_Vx, Velocity_Border_Soliton,U0_Schamel);
					Integrate_for_Num_Density(DF_one_Vx,Num_Den_Sub);
					if (Num_Den_Sub <  Num_Density_imposed_element){
						Velocity_Border_Soliton = Velocity_Border_Soliton - delta_Velocity;
						if (rank == 7){cout<<"$"<<flush;}
					}
					NumberIteration_InnerLevel += 1;
				}while ( Num_Den_Sub >  Num_Density_imposed_element && NumberIteration_InnerLevel<1000);
				NumberIteration_OuterLevel += 1;
				delta_Velocity = delta_Velocity/2.0;//P->dVx/(NumberIteration_OuterLevel);
			}while( (abs(Num_Den_Sub_0-Num_Density_imposed_element)/Num_Density_imposed_element) > (0.1/100.0) && NumberIteration_OuterLevel<100 );
			//==============================================================================================


			if (rank == 7){cout<<endl<<P->Name<<" *********  X="<<iX*P->dX+P->X_min<<" finish-->"<<"V_border="<<Velocity_Border_Soliton-U0_Schamel <<"  Num_Den="<<Num_Den_Sub <<"   "<< Num_Density_imposed_element<<endl<<flush;}

	}
	//-----------------------------------------------------------------------------------------------
 
 
 
	//-----------------------------------------------------------------------------------------------
	void Creat_One_Line_DF(vector<double>* DF_one_Vx, double& Velocity_Border_Soliton, double& U0_Schamel) {
			double Vx0, Fp0, X0;
			int iVx_displaced;
	// 		double random_number;
			vector<double> CNT(P->nVx);
			fill(CNT.begin(), CNT.end(), 0.0);
			fill(DF_one_Vx->begin(), DF_one_Vx->end(), 0.0);
			X0 = 0.0;
			if (P->Charge == 1.0){
				for (int iVx = 0; iVx<=P->nVx-1; iVx++) {
					for (int j = 1; j<=100; j++) 	  {
	// 					random_number = ((double) rand() / (double)(RAND_MAX)) ;
						Vx0 = P->Vx_min + double(iVx) * P->dVx +(double(j)/100.0)* P->dVx;
						Maxwellian_Reflected_distribution_Num_Den_Based(X0,Vx0,Fp0,Velocity_Border_Soliton,U0_Schamel);
						iVx_displaced = int( (Vx0 - P->Vx_min)/ P->dVx );
					
						if (Fp0 > 0){
							DF_one_Vx->at(iVx_displaced) += Fp0; 
							CNT[iVx_displaced] += 1.0;
						}
					}
				}
				for (int iVx = 0; iVx<=P->nVx-1; iVx++) {
					if (CNT[iVx]>0.0){DF_one_Vx->at(iVx) = DF_one_Vx->at(iVx)/CNT[iVx];}
				
				}
			}//ions
			if (P->Charge == -1.0){
				for (int iVx = 0; iVx<=P->nVx-1; iVx++) { // untrapped and pushed particles
					Vx0 = P->Vx_min + double(iVx) * P->dVx;
					Maxwellian_Trapped_distribution_Num_Den_Based(X0,Vx0,Fp0,Velocity_Border_Soliton,U0_Schamel);
					iVx_displaced = int( (Vx0 - P->Vx_min)/ P->dVx );
					if (Fp0 > 0){
					DF_one_Vx->at(iVx_displaced) = Fp0;
	// 				  if (iVx_displaced>=P->nVx){cout<<" out of -- "<<P->nVx<<" < "<<iVx_displaced<<" for Vx0="<<Vx0<<", iVx="<<iVx<<", Vx_min="<<P->Vx_max - P->dVx<<endl;}
					}
				}
				for (int iVx = int( (U0_Schamel - abs(Velocity_Border_Soliton-U0_Schamel) -P->Vx_min )/P->dVx) ; 
					iVx<=int( (U0_Schamel + abs(Velocity_Border_Soliton-U0_Schamel) -P->Vx_min )/P->dVx); iVx++) { //trapped particles
					if (iVx>0 && iVx<P->nVx){	
						Vx0 = P->Vx_min + double(iVx) * P->dVx;	
						Maxwellian_Jet_distribution(U0_Schamel, Fp0);
						DF_one_Vx->at(iVx) =Fp0;
					}else{
					cout<<"@";
					}

				}
			} //electrons
	}
	//-----------------------------------------------------------------------------------------------
  
	//-----------------------------------------------------------------------------------------------
	void Integrate_for_Num_Density(vector<double>* DF_one_Vx, double & Num_Den_Sub) {   
		Num_Den_Sub = 0.0;
		for (int iVx = 0; iVx<=P->nVx-1; iVx++) {
			if (iVx%2 == 0.0){ 
				Num_Den_Sub = Num_Den_Sub + 2.0 * DF_one_Vx->at(iVx);  
			}else{ 
				Num_Den_Sub = Num_Den_Sub + 4.0 * DF_one_Vx->at(iVx);   
			}
		}   
		Num_Den_Sub = (1.0/3.0) * P->dVx * Num_Den_Sub;
	}
	//-----------------------------------------------------------------------------------------------
  
	//-----------------------------------------------------------------------------------------------
	void Read_Phi_imposed(vector<double>& Num_Density_imposed_array_cpu) {
		#ifdef debug_initial
			P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="$$$ Read_Phi_imposed, Started $$$");
		#endif
	//       take note: arrays are starting from 0, 
	//       and if the pointers goes beyond the domain of the array, no error will be shown
	//       the type of the element has to agree with the type of the function 
	size_t flags = 0;
	json_error_t error;
	int i_element, X_integer;
	json_t* Pair_X_Den_json, *Num_Den_json, *X_json;
	json_t* file_json;
	
	vector<double> Num_Density_imposed_total;
	Num_Density_imposed_total.resize (P->nX);      
	fill(Num_Density_imposed_total.begin(), Num_Density_imposed_total.end(), 1.0);

	if (rank==0) cout << "-----------------------------------------------------"<<endl;
	
	//--- loading the JSON file ------------------------------------------
	if (P->Charge == 1.0){
			file_json = json_load_file("./Num_Den_ions_imposed.json", flags, &error);
	}else{
			file_json = json_load_file("./Num_Den_electrons_imposed.json", flags, &error);
	}
	if(!file_json){
		cout << "in line " << error.line << ": " << error.text <<flush<< endl;
		exit(10);      
	}
	//--------------------------------------------------------------------
	
	//--- returning the Num_Den object ------------------------------------  
	json_t* Num_Den_array_json = json_object_get(file_json, "Num_Den");
	if(!Num_Den_array_json){
		cerr << "in line " << error.line << ": " << error.text <<flush<<endl;
		exit(20);
	}
	size_t Num_Den_size =  json_array_size(Num_Den_array_json);
	int Num_Den_size_int = static_cast<int>(Num_Den_size);
	//        cout<<"size="<<Num_Den_size_int<<endl;
	i_element=0;
	for (int i_X = 0; i_X <= Num_Den_size_int-1; i_X++){
			Pair_X_Den_json    = json_array_get(Num_Den_array_json, i_element);
			if (!Pair_X_Den_json) cout<<"can't catch element "<< i_element<<" of Num_Den" <<flush<<endl;

			X_json = json_array_get(Pair_X_Den_json, 0);
			if (!X_json) cout<<"can't catch X of Num_Den" <<flush<<endl;
			double X_value = json_real_value(X_json);
	//   		if (rank==0) cout << "X of Grid" <<" = " << X_value << endl;
	
			Num_Den_json = json_array_get(Pair_X_Den_json, 1);
			if (!Num_Den_json) cout<<"can't catch Value of Num_Den" <<flush<<endl;
			double Num_Den_imposed_value = json_real_value(Num_Den_json);
	//    		if (rank==0) cout << "Num_Den of Grid" <<" = " << Num_Den_imposed_value << endl;
			X_integer = int((0.0+X_value+0.01*P->dX)/P->dX);
	//  		cout<<"{"<<X_value<<" "<< X_integer<<"  "<<flush;
			Num_Density_imposed_total[X_integer] = Num_Den_imposed_value;
			++i_element;

		
	}
	// 	  for (int i = 0; i <P->nX; i++){
	//  	    if (Num_Density_imposed_total[i]==0.0){ cout<<"  {"<<i<<"  "<<Num_Density_imposed_total[i]<<"}";}
	// 	  }
	

		if (rank == 0 ){
			for (int i_X = 0; i_X <= P->nX_cpu-1; i_X++){
			Num_Density_imposed_array_cpu [i_X]= Num_Density_imposed_total[i_X];
			}
			for (int i_rank = 1; i_rank <= size-1; i_rank++){
				int starting_nx = i_rank * P->nX_cpu;
				MPI_Send(&Num_Density_imposed_total[starting_nx] , P->nX_cpu, MPI_DOUBLE,i_rank,610, MPI_COMM_WORLD);
			}
		}else{
			MPI_Recv(&Num_Density_imposed_array_cpu[0], P->nX_cpu, MPI_DOUBLE,0 ,610, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

			
			
			
		#ifdef debug_initial
			P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="$$$ Read_Phi_imposed, finished $$$");
		#endif
	}
	//-----------------------------------------------------------------------------------------------
		
	//-----------------------------------------------------------------------------------------------
	void Injection(vector<double> & PP_X, vector<double> & PP_Vx, vector<double> & PP_Fp, int & iX, vector<int> & Num_PP_X,  vector<int> & Integer_PP_X_cpu) { 
	int integer_pp = Integer_PP_X_cpu[iX];
	double X0,Vx0,Fp0;
	double MPP0 = 0.0;
	int i_count =0;
	for (int i_pp=0; i_pp<=Num_PP_X[iX]-1; i_pp++ ){
	i_count++;
	X0  = PP_X[integer_pp+i_pp];
	Vx0 = PP_Vx[integer_pp+i_pp];
	Fp0 = PP_Fp[integer_pp+i_pp];
	PhasePointClass PhasePoint (X0,Vx0,Fp0,MPP0);
	if (X0>=P->X_min && X0<=P->X_min+(P->nX_cpu*P->dX)){
	// 	cout <<"rank="<<rank<<"   "<<P->X_min <<"   "<<P->X_min+(P->nX_cpu*P->dX)<<endl;
		PhasePoints->push_back(PhasePoint);
	}
	}
	//     if (rank==1){cout <<"rank="<<rank<<"  "<<iX<<"   "<<i_count<<"   "<<PhasePoints->getFill()<<endl;}
	}
	//-----------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------
	void Maxwellian_Schamel_distribution_Num_Den_Based (double & X0, double& Vx0, double& Fp0, double& Velocity_Border_Soliton,double& U0_Schamel,int & trapped_flag) {
		double random_number;
		double MPP0 = 0.0;
		if (P->Charge == 1.0){
			Maxwellian_Reflected_distribution_Num_Den_Based(X0,Vx0,Fp0,Velocity_Border_Soliton,U0_Schamel);} //ions
	
		if (P->Charge == -1.0){
			Maxwellian_Trapped_distribution_Num_Den_Based(X0,Vx0,Fp0,Velocity_Border_Soliton,U0_Schamel);
			//=============================================================================
			if (trapped_flag == 0){			
				for (int iVx = 0; iVx<=int( (Velocity_Border_Soliton - U0_Schamel) /P->dVx); iVx++) {
					for (int i = 1; i<=P->pX; i++) 	{
						for (int j = 1; j<=P->pVx; j++) {
							// upper side of the trapped population
							random_number = ((double) rand() / (double)(RAND_MAX)) ;
							Vx0 = U0_Schamel + double(iVx) * P->dVx +random_number * P->dVx; 
			
							random_number = ((double) rand() / (double)(RAND_MAX)) ;
							X0  =  double(int(X0/P->dX))  * P->dX  + random_number * P->dX;
							Maxwellian_Jet_distribution(U0_Schamel, Fp0);
							MPP0 = 1.0;
							//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
							if (Fp0>0){
								PhasePointClass PhasePoint_upper (X0,Vx0,Fp0,MPP0);
								PhasePoints->push_back(PhasePoint_upper);
							}
							//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
							
							// lower side of the trapped population
							random_number = ((double) rand() / (double)(RAND_MAX)) ;
							Vx0 = U0_Schamel - double(iVx) * P->dVx +random_number * P->dVx; 
			
							random_number = ((double) rand() / (double)(RAND_MAX)) ;
							X0  = double(int(X0/P->dX))  * P->dX  + random_number * P->dX;
							Maxwellian_Jet_distribution(U0_Schamel, Fp0);
							MPP0 = 2.0;
							//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
							if (Fp0>0){
								PhasePointClass PhasePoint_lower (X0,Vx0,Fp0,MPP0);
								PhasePoints->push_back(PhasePoint_lower);
							}
							//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
							
						}
					}
				}
			}//electrons
		}
		trapped_flag = 1;
	}
	//-----------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------
	void Maxwellian_Reflected_distribution_Num_Den_Based(double& X0, double& Vx0, double& Fp0, double& Velocity_Border_Soliton, double& U0_Schamel ) { 
		double Vx0_soliton_frame = Vx0-U0_Schamel; 
		double Vx_power_2 = pow(Vx0_soliton_frame,2) - pow(Velocity_Border_Soliton-U0_Schamel,2) ;
		Fp0 = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) 
			* exp( (  -0.5 * (pow(Vx0, 2)) ) * P->norm_factor  );

	
		if (Vx0_soliton_frame >= 0.0 ){ 
			if (Vx_power_2 > 0){Vx0_soliton_frame  = sqrt(Vx_power_2);}
			if (Vx_power_2 < 0){Vx0_soliton_frame  = -sqrt(-Vx_power_2);Fp0 = -1.0; } //the extra phase points which I think should stay  
		} else if (Vx0_soliton_frame < 0.0 ){
			if (Vx_power_2 > 0){Vx0_soliton_frame  = -sqrt(Vx_power_2);}
			if (Vx_power_2 < 0){Vx0_soliton_frame  = sqrt(-Vx_power_2); Fp0 = -1.0; }//(0.5*delta_IGP/Vx_phi) //(400.0 / sqrt(2.0 *   P->Charge * PHI_imposed / P->Mass ))
		} 
			
		Vx0 = Vx0_soliton_frame + U0_Schamel;
		if (Vx0 > P->Vx_max - P->dVx ){Fp0 = -1.0;}//Vx0  = P->Vx_max - P->dVx; }
		if (Vx0 < P->Vx_min + P->dVx ){Fp0 = -1.0;}//Vx0  = P->Vx_min + P->dVx; }
	
	
	
	
	
	// 	Fp0 = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) 
	// 			* exp( (  -0.5 * (pow(Vx0, 2)) ) * P->norm_factor  );
	// 			
	// 	
	// 	if (Vx0 < U0_Schamel ){
	// 		Vx0 = Vx0 + (abs(Velocity_Border_Soliton-U0_Schamel)* (- (exp(-pow(((Vx0-U0_Schamel)/1.0),10))-1.0)));
	//   		if (Vx0 > U0_Schamel){
	//    			Fp0 = -1.0;
	//   		}
	// 	}else if (Vx0 > U0_Schamel ){
	// 		Vx0 = Vx0 - (abs(Velocity_Border_Soliton-U0_Schamel)* (- (exp(-pow(((Vx0-U0_Schamel)/1.0),10))-1.0)));
	//  		if (Vx0 < U0_Schamel){
	//    			Fp0 = -1.0;
	//  		}
	// 	}
	// 	if (Vx0 > P->Vx_max - P->dVx ){Fp0 = -1.0; }
	// 	if (Vx0 < P->Vx_min + P->dVx ){Fp0 = -1.0; } 
		
	}
	//-----------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------
	void Maxwellian_Trapped_distribution_Num_Den_Based(double& X0, double& Vx0, double& Fp0, double& Velocity_Border_Soliton, double& U0_Schamel ) { 
		double Vx0_soliton_frame = Vx0-U0_Schamel; 
		double Vx_power_2 = pow(Vx0_soliton_frame,2) + pow(Velocity_Border_Soliton-U0_Schamel,2) ;
		Fp0 = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) 
			* exp( (  -0.5 * (pow(Vx0, 2)) ) * P->norm_factor  );

	
		if (Vx0_soliton_frame > 0.0 ){       
			Vx0_soliton_frame  = sqrt(Vx_power_2);
		} else if (Vx0_soliton_frame < 0.0  ){
			Vx0_soliton_frame  = -sqrt(Vx_power_2);
		}
			
		Vx0 = Vx0_soliton_frame + U0_Schamel;
		if (Vx0 > P->Vx_max - P->dVx ){Fp0 = -1.0;}//Vx0  = P->Vx_max - P->dVx; }
		if (Vx0 < P->Vx_min + P->dVx ){Fp0 = -1.0;}//Vx0  = P->Vx_min + P->dVx; }
		
	
	// 	Fp0 = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) 
	// 			* exp( (  -0.5 * (pow(Vx0, 2)) ) * P->norm_factor  );
	// 	if (Vx0 < U0_Schamel ){
	// 		Vx0 = Vx0 - abs(Velocity_Border_Soliton-U0_Schamel);
	// 	}
	// 	if (Vx0 > U0_Schamel ){
	// 		Vx0 = Vx0 + abs(Velocity_Border_Soliton-U0_Schamel);
	// 	}
	// 	if (Vx0 > (P->Vx_max - P->dVx) ){Fp0 = -1.0; }
	//	if (Vx0 < (P->Vx_min + P->dVx) ){Fp0 = -1.0; } 
	}
	//-----------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------
	void Maxwellian_Schamel_distribution(double & X0, double& Vx0, double& Fp0, vector<double>& PHI_Sagdeev,int & trapped_flag) { 
	double Beta_Schamel,U0_Schamel, Phi_imposed;
	double MPP0 = 0.0;
	
	double PHI_Sagdeev_max,PHI_imposed_width;
	int i_X_PHI_Sagdeev_max;
	vector<double>::iterator result;
	
	PHI_imposed_width = 1.0;
	result = max_element(PHI_Sagdeev.begin(), PHI_Sagdeev.end());
	i_X_PHI_Sagdeev_max = distance(PHI_Sagdeev.begin(), result);
	PHI_Sagdeev_max = PHI_Sagdeev[i_X_PHI_Sagdeev_max];
	
	
	//     double Charge_IGP = P->Charge;
	//     if (Charge_IGP > 0.0) {Charge_IGP = -Charge_IGP;}
	
	
	//             double A_soliton,X_soliton;
	//      if (X0<abs( P->X_min_Tot - P->X_max_Tot)/2.0){
	//         A_soliton = 0.04;X_soliton = 8.0;//abs( P->X_min_Tot - P->X_max_Tot)/2.0;
	//         Beta_Schamel = -6.7974, U0_Schamel = 0.4;
	//         PHI_Soliton_Mandel(Phi_imposed,X0,A_soliton,X_soliton,Beta_Schamel,U0_Schamel);
	
	//     Beta_Schamel = -1.0, U0_Schamel = 2.675; //if Vx0>0, otherwise U0 = -2.675
	//     PHI_Manfredi_test(Phi_imposed,X0,P->Alpha,P->Kappa);
	
	//      X_IGP = abs( P->X_min_Tot - P->X_max_Tot)/2; 
	double X_IGP = 512.0 ,delta_IGP , A_IGP, PHI_imposed;
	PHI_imposed = PHI_Sagdeev[ int(X0 / P->dX)];
	//     
	U0_Schamel = 8.8; 
	//     U0_Schamel = U0_Schamel /sqrt(P->norm_factor) ; 
	
	if (P->Charge == 1.0){
	Beta_Schamel = 0.0;A_IGP = 0.1;delta_IGP = 50.0;
	//       PHI_imposed = A_IGP * exp( -pow ( (X0  - X_IGP)/delta_IGP ,2) ) ;
	Maxwellian_Reflected_distribution(X0,Vx0,Fp0,PHI_imposed,PHI_imposed_width,U0_Schamel);} //ions
	
	if (P->Charge == -1.0){
	Beta_Schamel = 0.0;A_IGP = 0.1;delta_IGP = 50.0;
	//       PHI_imposed = A_IGP * exp( -pow ( (X0  - X_IGP)/delta_IGP ,2) ) ;
	Maxwellian_Trapped_distribution (X0,Vx0,Fp0,PHI_imposed,U0_Schamel);
	//=============================================================================
	if (trapped_flag == 0){
	// 	if (X0<X_IGP + 4.0 * delta_IGP && X0>X_IGP - 4.0 * delta_IGP ){
		for (int i = 1; i<=P->pX; i++) {
		double random_number = ((double) rand() / (double)(RAND_MAX)) ;
		double X0_trapped  = (float(int(X0/P->dX)) + random_number) * P->dX;
		PHI_imposed = PHI_Sagdeev[ int(X0_trapped / P->dX)];
	// 	      PHI_imposed = A_IGP * exp( -pow ( (X0_trapped  - X_IGP)/delta_IGP ,2) ) ;
		
	// 	      int pp_number = int(((PHI_imposed/A_IGP)+1.0) *100.0);
		int pp_number = int((PHI_imposed/PHI_Sagdeev_max+0.25) *100.0);
		double Vx_phi = sqrt(abs(2.0 *   P->Charge * PHI_imposed / P->Mass  ));
		if (pp_number != 0){
			for (int i_pp = 0; i_pp<= pp_number   ; i_pp++) { 
			Vx0 = U0_Schamel + i_pp * (Vx_phi/float(pp_number)); 
			//Maxwellian_Trapped_distribution (X0,Vx0,Fp0,A_IGP,X_IGP,delta_IGP,U0_Schamel);
			Maxwellian_Jet_distribution(U0_Schamel, Fp0);
			//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			PhasePointClass PhasePoint (X0_trapped,Vx0,Fp0,MPP0);
			PhasePoints->push_back(PhasePoint);
			//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			
			}
			for (int i_pp = 0; i_pp<= pp_number   ; i_pp++) { 
			Vx0 = U0_Schamel - i_pp * (Vx_phi/float(pp_number)); 
			//Maxwellian_Trapped_distribution (X0,Vx0,Fp0,A_IGP,X_IGP,delta_IGP,U0_Schamel);
			Maxwellian_Jet_distribution(U0_Schamel, Fp0);
			//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			PhasePointClass PhasePoint (X0_trapped,Vx0,Fp0,MPP0);
			PhasePoints->push_back(PhasePoint);
			//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			
			}
		}
		}
	trapped_flag = 1;
	}
	//=============================================================================
	}//electrons
	
	if (P->Charge == -3.0){//add hole to the DF   
		PHI_IGP(Phi_imposed,X0,A_IGP,X_IGP,delta_IGP);
		Phi_imposed = Phi_imposed * P->Temp/ (P->Charge ) ;
		// 	if (P->Mass == 1.0) {Beta_Schamel = -0.5;}
		// 	if (P->Mass == 0.01) {Beta_Schamel = -10.0;}
		double Energy_particle = ((P->norm_factor/2.0) * pow( (Vx0 -U0_Schamel) , 2)) - (P->Charge * Phi_imposed/P->Temp); //jenab this minus here doesn't have any physical backup!!! 

	
	
	
	
		Fp0 = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) ; //norm_factor = Mass/Temp;
	
		if		(Vx0 > (U0_Schamel + sqrt(abs(2.0*P->Charge * Phi_imposed/P->Mass))) ){//temp should be added
		Fp0 = Fp0 * exp(-pow( ( (sqrt(P->norm_factor/2.0) * U0_Schamel ) + sqrt(Energy_particle) ),2) ); //Maxwellian outside toward infinity
		}else if 	(Vx0 < (U0_Schamel - sqrt(abs(2.0*P->Charge * Phi_imposed/P->Mass))) ){
		Fp0 = Fp0 * exp(-pow( ( (sqrt(P->norm_factor/2.0) * U0_Schamel ) - sqrt(Energy_particle) ),2) ); // maxwellian around the zero velocity (inside DF)
		}else {//trapped particles
		Fp0 = Fp0 * exp(- ( ( (P->norm_factor/2.0) * ( pow(U0_Schamel,2) )) + Beta_Schamel * Energy_particle  )  );
		}
		
		
	//     }else{    
	// 	  A_soliton = 0.04;X_soliton = 40.0;//abs( P->X_min_Tot - P->X_max_Tot)*4.0/5.0;
	//           Beta_Schamel = -6.7974, U0_Schamel = -0.4;
	//           PHI_Soliton_Mandel(Phi_imposed,X0,A_soliton,X_soliton,Beta_Schamel,U0_Schamel);
	//       
	// 	  double Energy_particle = ((P->norm_factor/2.0) * pow( (Vx0 -U0_Schamel) , 2)) + (P->Charge * Phi_imposed/P->Temp);
	//   
	// 	  Fp0 = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) ; //norm_factor = Mass/Temp;
	//        
	// 	  if		(Vx0 > (U0_Schamel + sqrt(abs(2.0*P->Charge * Phi_imposed/P->Mass))) ){//temp should be added
	// 	      Fp0 = Fp0 * exp(-pow( ( (sqrt(P->norm_factor/2.0) * U0_Schamel ) + sqrt(Energy_particle) ),2) ); //Maxwellian outside toward infinity
	// 	  }else if 	(Vx0 < (U0_Schamel - sqrt(abs(2.0*P->Charge * Phi_imposed/P->Mass))) ){
	// 	      Fp0 = Fp0 * exp(-pow( ( (sqrt(P->norm_factor/2.0) * U0_Schamel ) - sqrt(Energy_particle) ),2) ); // maxwellian around the zero velocity (inside DF)
	// 	  }else {//trapped particles
	// 	      Fp0 = Fp0 * exp(- ( ( (P->norm_factor/2.0) * ( pow(U0_Schamel,2) )) + Beta_Schamel * Energy_particle  )  );
	// 	  }
	//     }
	
		}
	}
	//-----------------------------------------------------------------------------------------------
	*/
	//-----------------------------------------------------------------------------------------------
	void PHI_IGP (double& Phi_imposed,double& X0, double& A_IGP, double& X_IGP, double& delta_IGP) {  
	Phi_imposed = A_IGP * exp( -pow ( (X0  - X_IGP)/delta_IGP ,2) ) ;
	}
	//-----------------------------------------------------------------------------------------------  
	
	//-----------------------------------------------------------------------------------------------
	void PHI_Soliton_Mandel (double& Phi_imposed,double& X0,double& A_soliton,double& X_soliton, double Beta, double U0) {  
	double b = ( 1.0/sqrt(M_PI) ) * ( 1.0 - Beta - pow(U0,2) ) * exp(-0.5*pow(U0,2)); //jenab do we need to have norm_factor in this equation
	Phi_imposed = -P->Charge * A_soliton * 1.0/ pow(cosh(sqrt( (b/15.0) * sqrt(A_soliton))* (X0-X_soliton)),4);    
	}
	//-----------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------
	void PHI_Manfredi_test (double& Phi_imposed,double& X0,double& Alpha,double& Kappa) { 
	Phi_imposed = Alpha * sin(Kappa*X0)+1.0; //adding +1 make sure that phi is always positive.
	}
	//-----------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------
	void Maxwellian_Jet_distribution(double& Vx0, double& Fp0) { 

	Fp0 = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) 
			* exp(-0.5 * (pow(Vx0, 2)) * P->norm_factor); 
	Fp0 += (1.0-P->Den_Ratio_b)* sqrt(1.0/(2.0*M_PI)) * sqrt(4.0) 
			* exp(-0.5 * (pow(Vx0-4.5, 2)) * 4.0); 
			//double& Vx1_J,double& Temp1_J,double& Den1_J)       
	//     double Vx_J, Temp_J, Den_Ratio_J, norm_factor_J
	}
	//-----------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------
	void Maxwellian_Reflected_distribution(double& X0, double& Vx0, double& Fp0, double& PHI_imposed,double& PHI_imposed_width, double& U0_Schamel ) { 
	
	double Vx_phi = abs(2.0 *   P->Charge * PHI_imposed / P->Mass  );
	double Vx0_soliton_frame = Vx0-U0_Schamel; 
	double Vx_power_2 = pow(Vx0_soliton_frame,2) - Vx_phi ;
	Fp0 = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) 
			* exp( (  -0.5 * (pow(Vx0, 2)) ) * P->norm_factor  );
	// 	           * exp( (  -0.5 * (pow(Vx0, 2)) + 0.8125*P->Charge * PHI_imposed ) * P->norm_factor  ); 
	//     if (Vx0 > 0){Vx0 = Vx0 - sqrt(2.0*PHI_imposed) ;}
	//     if (Vx0 < 0){Vx0 = Vx0 + sqrt(2.0*PHI_imposed) ;}
	// 	  Vx0 = Vx0 * sqrt( 1 - (  (2.0 *   P->Charge * PHI_imposed) / (P->Mass * exp(-(pow(Vx0, 2))))  ));
	//        Vx0 = Vx0 *   (  1.0 - P->norm_factor*P->Charge *( (PHI_imposed/A_IGP)  *  exp(-(pow(Vx0, 2)))   ) );
	//       	   Vx0 = Vx0 *   (  1.0 - (P->Charge * PHI_imposed/P->Mass)  *  exp(-(pow(Vx0, 2)))   ) ;* exp(-(pow(Vx0 -U0_Schamel, 2)))
	
	if (Vx0_soliton_frame > 0.0 ){       
	if (Vx_power_2 > 0){Vx0_soliton_frame  = sqrt(Vx_power_2);}
	if (Vx_power_2 < 0){Vx0_soliton_frame  = -sqrt(-Vx_power_2); X0 = X0 + Vx0_soliton_frame*(PHI_imposed_width/sqrt(Vx_phi));} //the extra phase points which I think should stay
	} else
	if (Vx0_soliton_frame < 0.0  ){
	if (Vx_power_2 > 0){Vx0_soliton_frame  = -sqrt(Vx_power_2);}
	if (Vx_power_2 < 0){Vx0_soliton_frame  = sqrt(-Vx_power_2);X0 = X0 + Vx0_soliton_frame*(PHI_imposed_width/sqrt(Vx_phi));}//(0.5*delta_IGP/Vx_phi) //(400.0 / sqrt(2.0 *   P->Charge * PHI_imposed / P->Mass ))
	}
	Vx0 = Vx0_soliton_frame + U0_Schamel;
	if (Vx0 > P->Vx_max - P->dVx ){Vx0  = P->Vx_max - P->dVx; }
	if (Vx0 < P->Vx_min + P->dVx ){Vx0  = P->Vx_min + P->dVx; } 
	
	}
	//-----------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------
	void Maxwellian_Trapped_distribution(double& X0, double& Vx0, double& Fp0, double& PHI_imposed, double& U0_Schamel ) { 
	//     double PHI_imposed = A_IGP * exp( -pow ( (X0  - X_IGP)/delta_IGP ,2) ) ;  
	double Vx_phi = abs(2.0 *   P->Charge * PHI_imposed / P->Mass  );
	double Vx0_soliton_frame = Vx0-U0_Schamel;
	double Vx_power_2 = pow(Vx0_soliton_frame,2) + Vx_phi ;
	Fp0 = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) 
			* exp( (  -0.5 * (pow(Vx0, 2)) ) * P->norm_factor  );

	if (Vx0_soliton_frame > 0 ){       
	Vx0_soliton_frame  = sqrt(Vx_power_2);
	} else
	if (Vx0_soliton_frame < 0  ){
	Vx0_soliton_frame  = -sqrt(Vx_power_2);
	}
	Vx0 = Vx0_soliton_frame + U0_Schamel;
	if (Vx0 > P->Vx_max - P->dVx ){Vx0  = P->Vx_max - P->dVx; }
	if (Vx0 < P->Vx_min + P->dVx ){Vx0  = P->Vx_min + P->dVx; } 

	
	}
	//-----------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------
	void Kappa_Distribution(double& Vx0, double& Fp0, double& gamma1, double& gamma2) {
	Fp0 =  sqrt(1.0/(2.0*M_PI)) 	 *
		sqrt(1.0/(P->Kappa_DF-1.5)) * 
		sqrt(P->norm_factor) 	 * 
		(gamma1/gamma2) 		 * 
		pow( (1.0 + ( pow(Vx0,2) * 0.5 * P->norm_factor *  (1.0/(P->Kappa_DF-1.5)) )  ),(-P->Kappa_DF));  
	}
	//-----------------------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------------------
	void DF_add_randomness(double& Fp0,double& random_number_F) { 
	Fp0 = Fp0 * (1.0+ P->Alpha *(5.0)* random_number_F);
			//(P->Alpha * ( 1.0 + (5.0)* random_number_F )*cos(P->Kappa*X0))  
	}
	//-----------------------------------------------------------------------------------------------
	
		//-----------------------------------------------------------------------------------------------
		void DF_add_modes(double& X0_regular_grid,double& Vx0_regular_grid,  double& Fp0, double& EPP0) { 
			double Purtub_1,Purtub_0, Purtub;
			Purtub_0 = P->Alpha*cos(P->Kappa* (X0_regular_grid));
			Purtub_1 = P->Alpha*cos(P->Kappa* (X0_regular_grid+P->dX));
			Purtub = (Purtub_0 + Purtub_1)/2.0;
			
			EPP0 = 0.5 * (pow(Vx0_regular_grid, 2)) * P->norm_factor;
			if (X0_regular_grid> ((P->X_max_Tot - P->X_min_Tot)/2.0 - 5.0*P->dX) && X0_regular_grid< ((P->X_max_Tot - P->X_min_Tot)/2.0 + 5.0*P->dX))
			{ EPP0 = EPP0 + 2000.0; }
			Fp0 = Fp0 * (1.0+Purtub); 
		}
		//-----------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------
	void DF_add_3Modes(double& X0, double& Fp0) { 
	Fp0 = Fp0 * (1.0 + ( 
				(P->Alpha/3.0) * ( cos(1.0 * P->Kappa*X0) + cos(2.0 * P->Kappa*X0 + 0.25) + cos(3.0 * P->Kappa*X0 + 1.05) )
				)
			); 
	}
	//----------------------------------------------------------------------------------------------- 
	
	//-----------------------------------------------------------------------------------------------
	void DF_add_Envelope_Soliton(double& X0, double& Fp0) { 
	Fp0 = Fp0 * (1.0 + ( 
				(0.1) *exp(-pow(((X0-250.0)/20.8),2)) * ( cos(X0) )
				)
			); 
	}
	//-----------------------------------------------------------------------------------------------    

	void Phase_Point_DF_generator(double& X0, double& Vx0, double& Fp0, double & X0_regular_grid, double & Vx0_regular_grid, int & iX, double & EPP0, double & MPP0 ){
		Velocity_imposed = v_soliton_imposed_array_cpu[iX];		
		X0_CoMoving_correction = X0 + (Vx0-Velocity_imposed)* (0.0*P->d_Time); // this correction is to compensate the tilting effect in propgation of solitons
		iX_CoMoving_correction = int (( X0_CoMoving_correction  - P->X_min_Tot) / P->dX  );
		Phi_imposed = Potential_imposed_array_cpu[iX];
		if (iX_CoMoving_correction < 0 || iX_CoMoving_correction > P->nX ) {
			Phi_imposed_CoMoving_Correction = 0.0;
		}else{
			Phi_imposed_CoMoving_Correction = Potential_imposed_array_total[iX_CoMoving_correction];
		}
		
		Beta_imposed = Beta_soliton_imposed_array_cpu[iX];
		Phi_maximum_imposed = Phi_maximum_soliton_imposed_array_cpu[iX];
		
		Alpha_Schamel_Main_Soliton = Alpha_Schamel_Main_Soliton_array_cpu[iX];
		v_soliton_Other_Soliton = v_soliton_Other_Soliton_array_cpu[iX];
		phi_max_Other_Soliton = phi_max_Other_Soliton_array_cpu[iX];
		Alpha_Schamel_Other_Soliton = Alpha_Schamel_Other_Soliton_array_cpu[iX];
		
					
					//----------------------------------------------------------
		if (P->DF_NotShifted_Shifted == 0){
			if(P->DF_Vx_MK==0){
				Maxwellian_Jet_distribution(Vx0, Fp0);// Vx1_J, Temp1_J, Den1_J) // don't use Vx0_regular_grid in combination with open boundary condition, it creates asymmetry in upper and lower parts of Maxwellian, even in case of free streaming it can be clearly observed. 2018-12-11
			}else if(P->DF_Vx_MK==1){
				Kappa_Distribution(Vx0, Fp0, gamma1,gamma2);
			}
		}else if (P->DF_NotShifted_Shifted == 1) {
			Schamel_distribution_phi_imposed (X0,Vx0,Fp0,EPP0, MPP0,Phi_imposed,Phi_imposed_CoMoving_Correction, Beta_imposed, Phi_maximum_imposed, Velocity_imposed, Alpha_Schamel_Main_Soliton,v_soliton_Other_Soliton,phi_max_Other_Soliton,Alpha_Schamel_Other_Soliton, gamma1,gamma2, Energy_Kinentic_Jenab_DF, Beta_Jenab_DF);
			Schamel_DF_Flag = true;
		}
					//----------------------------------------------------------
				
					
		//----------------------------------------------------------
		if (Schamel_DF_Flag == false){
			if(P->DF_X_NMRE==1){
				DF_add_modes(X0_regular_grid,Vx0_regular_grid,Fp0,EPP0); // I am not sure if this works fine for Manfredi test 2018-12-11
			}else if(P->DF_X_NMRE==2){
				//DF_add_randomness(Fp0,random_number_F);
			}else if(P->DF_X_NMRE==3){
				DF_add_Envelope_Soliton(X0,Fp0);
				//DF_add_3Modes(X0,Fp0);
			}
		}
		//----------------------------------------------------------
	}
	//-------------------------------------------------------------------------
	
};


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
class SpeciesClass { 

public:
  ParametersClass * P;
	//---- MPI variables ------
	MpiClass * MPI;
  	int size, rank,NODE;
	int * Left_cpu;
	int * Right_cpu;
  
  vector<double> * DF_g_1D;
  #ifdef print_EPP
	vector<double> * EPP_g_1D;
  #endif
  #ifdef print_MPP
	vector<double> * MPP_g_1D;    
  #endif
  
  vector<double> * Num_Den_g; 
  vector<double> * Ene_kin_x_cpu;
  

  
  Storage<PhasePointClass> * PhasePoints; //creating the phase space class
  
  ofstream file_Num_Den, file_Time_S_F2, file_Time_S_FlnF, file_Time_Ene_kin; 
  
  //6) Entropy and Energy Parameters-------------------------------
  double S_F2_Tot, S_F2_cpu, S_F2_Tot_ini;
  double S_FlnF_Tot, S_FlnF_cpu, S_FlnF_Tot_ini;
  double Ene_kin_Tot, Ene_kin_cpu;
  
   vector<double> DF_g_RCV, DF_g_SND;
   #ifdef print_EPP
	vector<double> EPP_g_RCV, EPP_g_SND;
   #endif
   #ifdef print_MPP
	vector<double> MPP_g_RCV, MPP_g_SND;
   #endif


  //===============================================================================================
  SpeciesClass(ParametersClass * P_in, MpiClass * MPI_in){
    P = P_in;
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
  };
  //===============================================================================================
     
  //***********************************************************************************************
  void Initialization_SpeciesClass() {
    
    DF_g_RCV.resize(P->nVx+1);
    DF_g_SND.resize(P->nVx+1);
    
    #ifdef print_EPP
	EPP_g_RCV.resize(P->nVx+1);
	EPP_g_SND.resize(P->nVx+1);
    #endif
    #ifdef print_MPP
	MPP_g_RCV.resize(P->nVx+1);
	MPP_g_SND.resize(P->nVx+1);
    #endif
    
    DF_g_1D = new vector<double> ((P->nX_cpu+1)*P->nVx);
    
    #ifdef print_EPP
	EPP_g_1D = new vector<double> ((P->nX_cpu+1)*P->nVx);
    #endif
    #ifdef print_MPP
	MPP_g_1D = new vector<double> ((P->nX_cpu+1)*P->nVx);
    #endif
    
    Num_Den_g = new vector<double> (P->nX_cpu);
    Ene_kin_x_cpu = new vector<double> (P->nX_cpu);

    PhasePoints = new Storage<PhasePointClass>;
    
    stringstream pathname;
      
    #ifdef print_Temporal_Techplot
      pathname.str("");
      pathname << "./temporal/S_F2_"<<P->Name<<".dat";
      file_Time_S_F2.open(pathname.str().c_str());

      pathname.str("");
      pathname << "./temporal/S_FlnF_"<<P->Name<<".dat";
      file_Time_S_FlnF.open(pathname.str().c_str());

      pathname.str("");
      pathname << "./temporal/Ene_Kin_"<<P->Name<<".dat";
      file_Time_Ene_kin.open(pathname.str().c_str());
    #endif     

  }
  //***********************************************************************************************  
  
  

	void DF_INITIALIZATION_2() {
		
		//==================================================================================	
		double X0,Vx0,Fp0;
		double EPP0 =0.0;
		double MPP0 =0.0;
		double random_number;
		unsigned int seed_number;
		double X0_regular_grid, Vx0_regular_grid;
		
		// Data_PP = [Fp, X, X0, Vx, Vx0, EPP, MPP]
		vector<double> Data_store_PP;
		Data_store_PP.resize(P->number_of_data);
		
		#ifdef debug_initial
		P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="--- --- SpeciesClass: DF initialization");
		#endif
		DistributionFunctionClass DF_Object(P,MPI) ;
		#ifdef debug_initial
		P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="--- --- SpeciesClass: finished DF_Object");
		#endif
		seed_number = (unsigned)time(NULL)*(rank+1); // in order to have different sequece of pesudo-random number for each cpu, rank must be inside the formula for seeding
		srand(seed_number); //careful not to use srand() twice in the code unless you are sure that they are apart from each other in time at least for one second
	


        //=====================================================================
		double Fp_maximum = 0.0;
		for (int iX = 0; iX<=P->nX_cpu-1; iX++)  {
			for (int iVx = 0; iVx<=P->nVx-1; iVx++) {
				Vx0 = (P->Vx_min ) + double(iVx) * P->dVx;
				X0  = P->X_min  + double(iX)  * P->dX;
				DF_Object.Phase_Point_DF_generator(X0, Vx0,Fp0,X0,Vx0,iX,EPP0,MPP0);
				if (Fp0>Fp_maximum) Fp_maximum = Fp0;
			}
		}
		vector<double> Fp_maximum_array(size);
		MPI_Gather(&Fp_maximum, 1, MPI_DOUBLE,&Fp_maximum_array.at(0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (rank==0){
			for (int i_rank=0; i_rank<size; ++i_rank) {
				if(Fp_maximum>Fp_maximum_array[i_rank]) {
					Fp_maximum = Fp_maximum_array[i_rank];
				}
			}
		}
		MPI_Bcast(&Fp_maximum, 1,  MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //=====================================================================
		if (rank==0) {cout<<" >>> Fp_maximum="<<Fp_maximum<<" <<<"<<endl;}
						
		double log10_Fp_max;
		log10_Fp_max= log10(Fp_maximum);
		int number_of_PhasePoints;
		//=====================================================================
		
		for (int iX = 0; iX<=P->nX_cpu-1; iX++)  {
// 			if (rank==0){cout<<"-----"<<endl<<std::flush;}
			for (int iVx = 0; iVx<=P->nVx-1; iVx++) { 
				X0_regular_grid  = P->X_min  + double(iX)  * P->dX;
				Vx0_regular_grid =  (P->Vx_min ) + double(iVx) * P->dVx  ;// using this with open boundary condition might cause asymmetry in upper and lower part of DF
				DF_Object.Phase_Point_DF_generator(X0_regular_grid, Vx0_regular_grid,Fp0,X0_regular_grid,Vx0_regular_grid,iX,EPP0,MPP0);
				number_of_PhasePoints = int (25.0 / (abs(log10(Fp0))-abs(log10_Fp_max)+1.0));				
				if (number_of_PhasePoints<4) {number_of_PhasePoints = 4;}
// 				if (rank==0){cout<<number_of_PhasePoints<<" "<<std::flush;}
// 				for (int i = 0; i<number_of_PhasePoints; ++i) 	{
 				for (int i = 1; i<=P->pX; i++) 	{
 					for (int j = 1; j<=P->pVx; j++) 	  {
						random_number = ((double) rand() / (double)(RAND_MAX+1.0)) ;//create random number between [0,1)
						Vx0 =Vx0_regular_grid + random_number * P->dVx; //j*(dVx_p_s/4.0)// 
						 
						random_number = ((double) rand() / (double)(RAND_MAX+1.0)) ;
						X0  = X0_regular_grid  + random_number * P->dX;//i*(dX_p_s/4.0)//
 						
						DF_Object.Phase_Point_DF_generator(X0, Vx0,Fp0,X0_regular_grid,Vx0_regular_grid,iX,EPP0,MPP0);
									
						//IF (Fp0/= 0.0){ // adding the phase PhasePoints  
						if (Fp0 >= 0.0){
							//if (Fp0 >= 0.000001){ // jenab, serious change in the algorithm	
							
							/* Euler-Leapfrog data structure
							 * Data_store_PP[P->fp]=Fp0; 
							 * Data_store_PP[P->x]=X0;
							 * Data_store_PP[P->x0]=0.0;
							 * Data_store_PP[P->vx]= Vx0;
							 * Data_store_PP[P->vx0]=0.0; */
							
							Data_store_PP[P->fp]=Fp0; 
							Data_store_PP[P->x]=X0;
							Data_store_PP[P->vx]=Vx0;
							Data_store_PP[P->vx0]= 0.0;
							Data_store_PP[P->vx1]=0.0;
							#ifdef print_EPP
							Data_store_PP[P->epp]=EPP0;
							#endif
							
							#ifdef print_MPP
							Data_store_PP[P->mpp]=MPP0;	
							#endif
							PhasePointClass PhasePoint (Data_store_PP, P->number_of_data);
							
	// 						PhasePointClass PhasePoint (X0,Vx0,Fp0,EPP0,MPP0);
							//PhasePoint= new PhasePointClass(X0,Vx0,Fp0);
							PhasePoints->push_back(PhasePoint);
	// 						delete PhasePoint;
						}
 					}
				}
			}
		}
	}
	//***********************************************************************************************

  
	
	//***********************************************************************************************
	void Read_Num_PhasePoints(vector<int> & Num_PP_X, int & time_name){
	hid_t       File_ID_Read, Group_ID_Read, dataset_ID_Read;          //handles 
	hid_t       memspace_Read;
	hid_t       dataspace_Read; //datatype_Read
	int rank_read;
	
	hsize_t Dim_DFp_Chunk[1];
	hsize_t Count_DFp[1], OffSet_DFp[1], Stride_DFp[1], Block_DFp[1],OffSet_DFp_dump[1];
	//     hsize_t dimsm_Read[1];              // memory space dimensions 
	

	
	
	Dim_DFp_Chunk [0] = P->nX_cpu; 
	OffSet_DFp    [0] = rank*P->nX_cpu;

	OffSet_DFp_dump[0] = 0;
	fill(Num_PP_X.begin(), Num_PP_X.end(), 0);
	
	
	
	stringstream Group_name;
	Group_name << "Timestep_"<<time_name;  
	
	//-----------------------------------------------------------------
	File_ID_Read = H5Fopen("Contour_feeding.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	Group_ID_Read = H5Gopen2( File_ID_Read, Group_name.str().c_str(), H5P_DEFAULT ); 
	//-----------------------------------------------------------------
	
	stringstream DSet_name;
	DSet_name << "Num_"<<P->Name; 
	dataset_ID_Read = H5Dopen(Group_ID_Read, DSet_name.str().c_str(),H5P_DEFAULT);

	Block_DFp [0] = Dim_DFp_Chunk [0];
	Count_DFp[0]= 1;
	Stride_DFp[0]=1;
	//     cout << fixed <<" rank = "<<rank<<"   from "<<OffSet_DFp [0]<<"   end at "<<OffSet_DFp [0]+Dim_DFp_Chunk [0]-1<<endl;
	
	
	

	
	dataspace_Read = H5Dget_space(dataset_ID_Read);
	rank_read      = H5Sget_simple_extent_ndims(dataspace_Read);
	
	P->status_hdf5 = H5Sselect_hyperslab(dataspace_Read, H5S_SELECT_SET, OffSet_DFp, Stride_DFp, Count_DFp, Block_DFp);
	if (P->status_hdf5<0){cout <<" --- Reading-Feeding error, not able to select dataspace_Read in rank = "<<rank<<endl;}
	//     dimsm_Read[0]=64;//dims_out[0];
	//     memspace_Read = H5Screate_simple(rank_read,dimsm_Read,NULL);
	memspace_Read = H5Screate_simple(rank_read,Dim_DFp_Chunk,NULL);

	
	P->status_hdf5 = H5Sselect_hyperslab(memspace_Read, H5S_SELECT_SET, OffSet_DFp_dump, Stride_DFp, Count_DFp, Block_DFp);
	if (P->status_hdf5<0){cout <<" --- Reading-Feeding error, not able to select memspace_Read in rank = "<<rank<<endl;}
	
	
	P->status_hdf5 = H5Dread(dataset_ID_Read, H5T_NATIVE_INT, memspace_Read, dataspace_Read, H5P_DEFAULT, Num_PP_X.data());  
	if (P->status_hdf5<0){cout <<" --- Reading-Feeding error, not able to read in rank = "<<rank<<endl;}


	
	H5Dclose(dataset_ID_Read);
	
	//-----------------------------------------------------------------
	P->status_hdf5 = H5Gclose(Group_ID_Read);
	H5Fclose(File_ID_Read);
	//-----------------------------------------------------------------

	
	}
	//***********************************************************************************************


  //***********************************************************************************************
  void Read_X_Vx_Fp(vector<double> & PP_X, vector<double> & PP_Vx, vector<double> & PP_Fp, vector<int> & Integer_PP_X_cpu, vector<int> & Num_PP_X, int & time_name){
  
    hid_t       File_ID_Read, Group_ID_Read, dataset_ID_Read;         // handles 
    hid_t       memspace_Read;
    hid_t       dataspace_Read; //datatype_Read

    int rank_read;
  
    hsize_t Dim_DFp_Chunk[1];
    hsize_t Count_DFp[1], OffSet_DFp[1], Stride_DFp[1], Block_DFp[1],OffSet_DFp_dump[1];
    hsize_t dimsm_Read[1];              // memory space dimensions 
  
  
    vector<int> Num_PP_cpu_array(size);
    vector<int> cpu_PP_border_array(size);
    
    int Num_PP_cpu = 0;
    fill(Integer_PP_X_cpu.begin(), Integer_PP_X_cpu.end(), 0);
    for (int i = 0; i<=P->nX_cpu-1; ++i){
        Num_PP_cpu += Num_PP_X[i];
	
    }
    int temp=0;
    for (int i = 1; i<=P->nX_cpu-1; ++i){
        temp += Num_PP_X[i-1];
	Integer_PP_X_cpu [i]= temp;
	
    }

    //cout<<endl;
    //cout <<"rank="<<rank<<"  "<<Num_PP_cpu;
  
  
  
    MPI_Gather(&Num_PP_cpu,          1, MPI_INT,
	    &Num_PP_cpu_array.at(0), 1, MPI_INT,
	    0, MPI_COMM_WORLD);
    if (rank==0) {
      int temp= 0;
      fill(cpu_PP_border_array.begin(), cpu_PP_border_array.end(), 0);
      for (int i=1; (i<=NODE); ++i){
 	      temp += Num_PP_cpu_array.at(i-1);
	      cpu_PP_border_array[i] = temp;
      }
//       for (int i=0; (i<=NODE); ++i){ cout <<"rank="<<rank<<"  "<<cpu_PP_border_array[i]<<endl; }
      
      
    }
 
    MPI_Bcast(&cpu_PP_border_array[0], size,  MPI_INT, 0, MPI_COMM_WORLD);

//      for (int i = 0; i<=P->nX_cpu-1; ++i){ if (rank==1) {cout<<"   iX=" <<i<<"   Num per iX="<<Num_PP_X[i]<<"   integer="<<Integer_PP_X_cpu [i]<<endl;}}
    
    PP_X.resize(Num_PP_cpu);
    fill(PP_X.begin(), PP_X.end(), 0.0);
   
    PP_Vx.resize(Num_PP_cpu);
    fill(PP_Vx.begin(), PP_Vx.end(), 0.0);
      
    PP_Fp.resize(Num_PP_cpu);
    fill(PP_Fp.begin(), PP_Fp.end(), 0.0);

    Dim_DFp_Chunk [0] = Num_PP_cpu;
    OffSet_DFp[0]= cpu_PP_border_array[rank];

    dimsm_Read[0]=Num_PP_cpu;
    OffSet_DFp_dump[0]=0;
  
    stringstream Group_name;
    Group_name << "Timestep_"<<time_name; 
  
    
    //-----------------------------------------------------------------
    File_ID_Read = H5Fopen("Contour_feeding.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    Group_ID_Read = H5Gopen2( File_ID_Read, Group_name.str().c_str(), H5P_DEFAULT ); 
  
    Block_DFp [0] = Dim_DFp_Chunk [0];
    Count_DFp[0]= 1;
    Stride_DFp[0]=1;
    
    
    //-----------------------------------------------------------------------------------------------------------------  
    stringstream DSet_name;
    DSet_name << "X_"<<P->Name;
  
    dataset_ID_Read = H5Dopen(Group_ID_Read, DSet_name.str().c_str(),H5P_DEFAULT);
    dataspace_Read = H5Dget_space(dataset_ID_Read);
    rank_read      = H5Sget_simple_extent_ndims(dataspace_Read);
    
    P->status_hdf5 = H5Sselect_hyperslab(dataspace_Read, H5S_SELECT_SET, OffSet_DFp, Stride_DFp, Count_DFp, Block_DFp);
    if (P->status_hdf5<0){cout <<"------- error, not able to select dataspace_Read in rank = "<<rank<<endl;}
    
    memspace_Read = H5Screate_simple(rank_read,dimsm_Read,NULL);
    P->status_hdf5 = H5Sselect_hyperslab(memspace_Read, H5S_SELECT_SET, OffSet_DFp_dump, Stride_DFp, Count_DFp, Block_DFp);
    if (P->status_hdf5<0){cout <<"------- error, not able to select memspace_Read in rank = "<<rank<<endl;}
    
    P->status_hdf5 = H5Dread(dataset_ID_Read, H5T_NATIVE_DOUBLE, memspace_Read, dataspace_Read, H5P_DEFAULT, PP_X.data());  
    if (P->status_hdf5<0){cout <<"------- error, not able to read in rank = "<<rank<<endl;}
    
    H5Dclose(dataset_ID_Read);
    //-----------------------------------------------------------------------------------------------------------------

    //-----------------------------------------------------------------------------------------------------------------  
    DSet_name.str("");
    DSet_name << "Vx_"<<P->Name;

    dataset_ID_Read = H5Dopen(Group_ID_Read, DSet_name.str().c_str(),H5P_DEFAULT);
    dataspace_Read = H5Dget_space(dataset_ID_Read);
    rank_read      = H5Sget_simple_extent_ndims(dataspace_Read);
    
    P->status_hdf5 = H5Sselect_hyperslab(dataspace_Read, H5S_SELECT_SET, OffSet_DFp, Stride_DFp, Count_DFp, Block_DFp);
    if (P->status_hdf5<0){cout <<"------- error, not able to select dataspace_Read in rank = "<<rank<<endl;}
    
    memspace_Read = H5Screate_simple(rank_read,dimsm_Read,NULL);
    P->status_hdf5 = H5Sselect_hyperslab(memspace_Read, H5S_SELECT_SET, OffSet_DFp_dump, Stride_DFp, Count_DFp, Block_DFp);
    if (P->status_hdf5<0){cout <<"------- error, not able to select memspace_Read in rank = "<<rank<<endl;}
    
    P->status_hdf5 = H5Dread(dataset_ID_Read, H5T_NATIVE_DOUBLE, memspace_Read, dataspace_Read, H5P_DEFAULT, PP_Vx.data());  
    if (P->status_hdf5<0){cout <<"------- error, not able to read in rank = "<<rank<<endl;}
 
    H5Dclose(dataset_ID_Read);
    //-----------------------------------------------------------------------------------------------------------------    
  
    
    
    //-------------------------------------------------------------------------------------------------------------------  
    DSet_name.str("");
    DSet_name << "Fp_"<<P->Name;

    dataset_ID_Read = H5Dopen(Group_ID_Read, DSet_name.str().c_str(),H5P_DEFAULT);
    dataspace_Read  = H5Dget_space(dataset_ID_Read);
    rank_read       = H5Sget_simple_extent_ndims(dataspace_Read);
    
    P->status_hdf5 = H5Sselect_hyperslab(dataspace_Read, H5S_SELECT_SET, OffSet_DFp, Stride_DFp, Count_DFp, Block_DFp);
    if (P->status_hdf5<0){cout <<"------- error, not able to select dataspace_Read in rank = "<<rank<<endl;}
    
    memspace_Read = H5Screate_simple(rank_read,dimsm_Read,NULL);
    P->status_hdf5 = H5Sselect_hyperslab(memspace_Read, H5S_SELECT_SET, OffSet_DFp_dump, Stride_DFp, Count_DFp, Block_DFp);
    if (P->status_hdf5<0){cout <<"------- error, not able to select memspace_Read in rank = "<<rank<<endl;}
    
    P->status_hdf5 = H5Dread(dataset_ID_Read, H5T_NATIVE_DOUBLE, memspace_Read, dataspace_Read, H5P_DEFAULT, PP_Fp.data());  
    if (P->status_hdf5<0){cout <<"------- error, not able to read in rank = "<<rank<<endl;}
  
    H5Dclose(dataset_ID_Read);
    //----------------------------------------------------------------------------------------------------------------- 
    
    
    
    //-----------------------------------------------------------------
    P->status_hdf5 = H5Gclose(Group_ID_Read);
    H5Fclose(File_ID_Read);
    //-----------------------------------------------------------------     
}
//***********************************************************************************************



      
  //***********************************************************************************************
  void Interpolate_Integrate(vector<double>& Den_g) {
        
      
//         Interpolation_Bilinear ();
	bilinear_weighting_Kazeminezhad();
        DF_Boundary_correction_for_iX_0 ();
	DF_Boundary_correction_for_iX_nX ();
    
        
         Num_Density_Calculator();
//         Num_Density_Calculator_zeor_order();
        for (int iX=0; iX<=P->nX_cpu-1; iX++){
                Den_g[iX] = P->Charge*P->Density*Num_Den_g->at(iX);
        }
  }
  //***********************************************************************************************

  //***********************************************************************************************
  void T_nD_1D (int & i_1D,int & i_X,int & i_Vx){
      i_1D = i_X * (P->nVx) + i_Vx;    
  }
  //***********************************************************************************************

  //***********************************************************************************************
  void T_1D_nD (int & i_1D,int & i_X,int & i_Vx){
      i_X = int (i_1D/(P->nVx));
      i_Vx = i_1D%(P->nVx);
  }
  //***********************************************************************************************
  
  
	//*********************************************************************
	void bilinear_weighting_Kazeminezhad(){
		int lX,lVx, i_Fg_ii;
		double initial_value	= 0.0;
		int size_shadow_x	= 2;
		int size_shadow_Vx	= 2;
		std::vector<int> i_X_array, i_Vx_array;
		i_X_array.resize (size_shadow_x,   0);
		i_Vx_array.resize(size_shadow_Vx,  0);
                
		std::vector<double>	Distance_X_array, Distance_Vx_array;
		Distance_X_array.resize (size_shadow_x,   initial_value);
		Distance_Vx_array.resize(size_shadow_Vx,  initial_value);
                

                        
		std::vector<std::vector<double>> Weight_PP_matrix;
		Weight_PP_matrix.resize(size_shadow_x, std::vector<double>(size_shadow_Vx, initial_value));
		
		vector<double> wg_PP_1D (P->nX_cpu*P->nVx);	
		
		fill(wg_PP_1D.begin(), wg_PP_1D.end(), 0.0);
		fill(DF_g_SND.begin(), DF_g_SND.end(), 0.0);
		fill(DF_g_1D->begin(), DF_g_1D->end(), 0.0);
		
		#ifdef print_EPP
		vector<double> CNT_EPP (P->nX_cpu*P->nVx);
		fill(CNT_EPP.begin(),   CNT_EPP.end(),   0.0);
		fill(EPP_g_SND.begin(), EPP_g_SND.end(), 0.0);
		fill(EPP_g_1D->begin(), EPP_g_1D->end(), 0.0);
		#endif
	
		#ifdef print_MPP
		vector<double> CNT_MPP (P->nX_cpu*P->nVx);
		fill(CNT_MPP.begin(),   CNT_MPP.end(),   0.0);
		fill(MPP_g_SND.begin(), MPP_g_SND.end(), 0.0);
		fill(MPP_g_1D->begin(), MPP_g_1D->end(), 0.0);
		#endif
		
		//-------------------------------------------------------------
		for(Storage<PhasePointClass>::Iterator particle_iter = PhasePoints->begin(); particle_iter != PhasePoints->end(); ++particle_iter) {
			lX			= int ( ((*particle_iter).Data_PP->at(P->x)  - P->X_min ) / P->dX  );
			lVx 			= int ( ((*particle_iter).Data_PP->at(P->vx_interpol) - P->Vx_min) / P->dVx );
			
			Distance_X_array [0]	= (abs ( ((*particle_iter).Data_PP->at(P->x)  - P->X_min ) - (double(lX)  * P->dX  ) ) )/P->dX;
			Distance_X_array [1]	= 1.0 - Distance_X_array[0]; 
			Distance_Vx_array[0]	= (abs ( ((*particle_iter).Data_PP->at(P->vx_interpol) - P->Vx_min) - (double(lVx) * P->dVx ) ) )/P->dVx;
			Distance_Vx_array[1]	= 1.0 - Distance_Vx_array[0]; 
			

                        
			for (int iX=0; iX<size_shadow_x; ++iX){
				for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
					Weight_PP_matrix[iX][iVx] = Distance_X_array[iX] * Distance_Vx_array[iVx];
				}
			}
			/*
			double Wtotal = 0.0;
			for (int iX=0; iX<size_shadow_x; ++iX){
				for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
					Wtotal += Weight_PP_matrix[iX][iVx];
				}
			} // wtotal = 1
			//if (rank==0) cout<<Wtotal-1.0<<"  "<<std::flush;
			for (int iX=0; iX<size_shadow_x; ++iX){
				for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
					Weight_PP_matrix[iX][iVx] = Weight_PP_matrix[iX][iVx]/Wtotal;
				}
			}
			*/
			i_X_array[0] = lX;
			int rX = lX;
			for (int iX=1; iX<size_shadow_x; ++iX){
				i_X_array[iX] = rX;
				++rX;
			}
			
			i_Vx_array[0] = lVx;
			int rVx = P->VxUp[lVx];
			for (int iVx=1; iVx<size_shadow_Vx; ++iVx){
				i_Vx_array[iVx] = rVx;
				rVx = P->VxUp[rVx];
			}			
			for (int iX=0; iX<size_shadow_x; ++iX){
				for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
					T_nD_1D (i_Fg_ii,i_X_array[iX],i_Vx_array[iVx]);
					DF_g_1D->at(i_Fg_ii)	+= Weight_PP_matrix[iX][iVx] * (*particle_iter).Data_PP->at(P->fp);//(*particle_iter).Data_PP->at(P->fp);
					wg_PP_1D[i_Fg_ii]	+= Weight_PP_matrix[iX][iVx];//1.0;
					
					#ifdef print_EPP
					EPP_g_1D->at(i_Fg_ii)	+=      (*particle_iter).Data_PP->at(P->epp);
					CNT_EPP [i_Fg_ii]	+= 1.0;
					#endif
					
					#ifdef print_MPP
					MPP_g_1D->at(i_Fg_ii)	+=      (*particle_iter).Data_PP->at(P->mpp);
					CNT_MPP [i_Fg_ii]	+= 1.0;
					#endif
				}
			}
		}
		//-------------------------------------------------------------
		

		//-------------------------------------------------------------
		for (int iX = 0; iX<P->nX_cpu; ++iX) {
			for (int iVx = 0; iVx<P->nVx; ++iVx) {
				T_nD_1D (i_Fg_ii,iX,iVx);			
				//----------------------------------------------------------------------
				if (wg_PP_1D[i_Fg_ii]!=0.0) {
					DF_g_1D->at(i_Fg_ii)	 = DF_g_1D->at(i_Fg_ii)/wg_PP_1D[i_Fg_ii];
				} else {
					//cout<<"zero grid point of"<<P->Name<<": iX="<<iX<<"  iVx="<<iVx<<endl;
					DF_g_1D->at(i_Fg_ii)	 = 0.0;    
				}
				//----------------------------------------------------------------------
				
				//----------------------------------------------------------------------
				#ifdef print_EPP
				if (CNT_EPP[i_Fg_ii]!=0.0) {
					EPP_g_1D->at(i_Fg_ii)	 = EPP_g_1D->at(i_Fg_ii)/CNT_EPP[i_Fg_ii];
				} else {
					EPP_g_1D->at(i_Fg_ii)	 = 0.0;
				}
				#endif
				#ifdef print_MPP
				if (CNT_MPP[i_Fg_ii]!=0.0) {
					MPP_g_1D->at(i_Fg_ii) 	 = MPP_g_1D->at(i_Fg_ii)/CNT_MPP[i_Fg_ii];
				} else {
					MPP_g_1D->at(i_Fg_ii)	 = 0.0; 
				}
				#endif
				//----------------------------------------------------------------------
			}
		}
		//-------------------------------------------------------------
		
		
		//-------------------------------------------------------------
		for (int iVx = 0; iVx<P->nVx; ++iVx) {
			T_nD_1D (i_Fg_ii,P->nX_cpu,iVx);
			DF_g_SND[iVx] = DF_g_1D->at(i_Fg_ii);
			#ifdef print_EPP
			EPP_g_SND[iVx] = EPP_g_1D->at(i_Fg_ii);
			#endif
			#ifdef print_MPP
			MPP_g_SND[iVx] = MPP_g_1D->at(i_Fg_ii);
			#endif
		}
		//-------------------------------------------------------------
		
	}
	//*********************************************************************
	
  //***********************************************************************************************
  void Interpolation_Bilinear (){
  //Inverse distance weighting
	double Wtotal;
	int size_shadow_Vx = P->size_shadow_Vx;
	int size_shadow_x = 2;
	double initial_value = 0.0;
// 		left,up (0,1) -------- right, up (1,1)
// 			|			|
// 			|			|
// 			|	PP*		|
// 			|			|
// 			|			|
// 		left,down (0,0) -------- right, down (1,0)
							// 		      (0,3) ----------------- (1,3)
							// 			|			|
							// 			|			|
							// 			|			|
							// 			|			|
							// 			|			|
							// 		      (0,2) ----------------- (1,2)	
							// 	Velocity	|			|
							// 		Dir	|			|
							// 			|	PP*		|
							// 			|			|
							// 			|			|
							// 		      (0,1) ----------------- (1,1)
							// 			|			|
							// 			|			|
							// 			|			|
							// 			|			|
							// 			|			|
							// 		      (0,0) ----------------- (1,0)	
	
	std::vector<std::vector<double>> Weight_PP_matrix;
	Weight_PP_matrix.resize(size_shadow_x, std::vector<double>(size_shadow_Vx, initial_value));
	
	std::vector<double> Distance_X_array, Distance_Vx_array;
	Distance_X_array.resize (size_shadow_x,   initial_value);
	Distance_Vx_array.resize(size_shadow_Vx,  initial_value);
	
	std::vector<int> i_X_array, i_Vx_array;
	i_X_array.resize (size_shadow_x,   0);
	i_Vx_array.resize(size_shadow_Vx,  0);
	
	std::vector<int> value_iVx, coeff_iVx;
	value_iVx.resize(size_shadow_Vx,  0);
	coeff_iVx.resize(size_shadow_Vx,  0);
	
	std::vector<std::vector<int>> i_Fg_matrix;
	i_Fg_matrix.resize(size_shadow_x, std::vector<int>(size_shadow_Vx, 0));
	int lX,lVx;
	int i_Fg_ii, half_shadow_size;
	double Vx_in_cell, Distance_Vx_closest, largest_distance, half_shadow, Distance_Vx_opposite ;
    
    
	vector<double> CNT ((P->nX_cpu+1)*P->nVx);
	
	
	fill(CNT.begin(), CNT.end(), 0.0);
	

	fill(DF_g_SND.begin(), DF_g_SND.end(), 0.0);
	fill(DF_g_1D->begin(), DF_g_1D->end(), 0.0);

	#ifdef print_EPP
		vector<double> CNT_EPP ((P->nX_cpu+1)*P->nVx);
		fill(CNT_EPP.begin(), CNT_EPP.end(), 0.0);
		fill(EPP_g_SND.begin(), EPP_g_SND.end(), 0.0);
		fill(EPP_g_1D->begin(), EPP_g_1D->end(), 0.0);
	#endif
	
	#ifdef print_MPP
		vector<double> CNT_MPP ((P->nX_cpu+1)*P->nVx);
		fill(CNT_MPP.begin(), CNT_MPP.end(), 0.0);
		fill(MPP_g_SND.begin(), MPP_g_SND.end(), 0.0);
		fill(MPP_g_1D->begin(), MPP_g_1D->end(), 0.0);
	#endif
	// size = 1  2  3  4  5  6  7 ...
	// half = 0  0  1  1  2  2  3 ...		
					//size=5// iVx       0    1  2  3   4 	//size=6// iVx       0    1  2  3   4   5
					//half=2// value_0  -2   -1  0  1   2 	//half=2// value_0  -2   -1  0  1   2   3
						// value     2    1  0  1   2 		// value     2    1  0  1   2   3
						// coeff     +    +  +  -   - 	 	// coeff     +    +  +  -   -   - 
					// Distance_Vx      2+d  1+d d 1-d 2-d			    2+d  1+d d 1-d 2-d 3-d
	half_shadow_size = int((float(size_shadow_Vx)/2.0)-0.25); // -0.25 is to make sure it works for both odd and even numbers
	for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
		value_iVx[iVx]	= iVx - half_shadow_size; // value_0
	}
	for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
		coeff_iVx[iVx] = 1;
		if (value_iVx[iVx]>0){
			coeff_iVx[iVx] = -1;			
		}		
	}
	for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
		value_iVx[iVx] = abs(value_iVx[iVx]);
	}
	
	
	largest_distance = sqrt(pow(float(size_shadow_Vx),2)+ pow (float(size_shadow_x),2));  
			//sqrt( pow( (int(float(size_shadow_Vx/2.0))), 2) + pow (1.0,2));
	for(Storage<PhasePointClass>::Iterator particle_iter = PhasePoints->begin(); particle_iter != PhasePoints->end(); ++particle_iter) { //interpolation main loop
		lX  = int (( (*particle_iter).Data_PP->at(P->x)  - P->X_min  ) / P->dX  ); //lx can be used to find out of boundry error for cpus
		Vx_in_cell = ((*particle_iter).Data_PP->at(P->vx) - P->Vx_min);
		lVx = int ( Vx_in_cell / P->dVx );
		
// 		if (lX<0 || lX>P->nX_cpu ) { 
// 			cout << "out of boundary in X =" << (*particle_iter).Data_PP->at(P->x) << "  lX= " << lX << "-- for ="<< P->Name <<" ---rank="<<rank<< endl;
// 		}
// 		if (lVx<0 || lVx>P->nVx ){ 
// 			cout << "out of boundary in Vx= " << (*particle_iter).Data_PP->at(P->vx) << " lVx="<< lVx << "-- for ="<< P->Name<<" ---rank="<<rank<<" X =" << (*particle_iter).Data_PP->at(P->x) << "  lX= " << lX <<endl;
// 		}
		Distance_X_array[0]=  ( abs( ((*particle_iter).Data_PP->at(P->x) -  P->X_min  ) - (double(lX)  * P->dX ) ) ) / P->dX;
		Distance_X_array[1]= 1.0-Distance_X_array[0]; 
		
		Distance_Vx_closest = ( abs( Vx_in_cell - (double(lVx) * P->dVx) ) ) / P->dVx;
		for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
			Distance_Vx_array[iVx] = value_iVx[iVx] + (coeff_iVx[iVx] * Distance_Vx_closest);			
		}
// 		Distance_Vx_opposite = 1.0 - Distance_Vx_closest;
// 		half_shadow = int(float(size_shadow_Vx/2.0));
		
// 		for (int iVx=0; iVx<half_shadow; ++iVx){
// 			Distance_Vx_array[iVx] = Distance_Vx_closest + (iVx + (half_shadow-1));//abs(Distance_Vx_closest + (int(float(size_shadow_Vx/2.0)-1) - iVx));
// 		}
// 		for (int iVx=half_shadow; iVx<size_shadow_Vx; ++iVx){
// 			Distance_Vx_array[iVx] = Distance_Vx_opposite + (half_shadow-iVx);//abs(Distance_Vx_closest + (int(float(size_shadow_Vx/2.0)-1) - iVx));
// 		}
//   		Distance_Vx_array[0] = Distance_Vx_closest;
//   		Distance_Vx_array[1] = Distance_Vx_opposite;


// 		Weight_PP_matrix[0][0] = sqrt(2.0)-(sqrt(pow (Distance_X_array[0],2) + pow(Distance_Vx_array[0],2)));
// 		Weight_PP_matrix[0][1] = sqrt(2.0)-(sqrt(pow (Distance_X_array[0],2) + pow(Distance_Vx_array[1],2)));
// 		Weight_PP_matrix[1][0] = sqrt(2.0)-(sqrt(pow (Distance_X_array[1],2) + pow(Distance_Vx_array[0],2)));
// 		Weight_PP_matrix[1][1] = sqrt(2.0)-(sqrt(pow (Distance_X_array[1],2) + pow(Distance_Vx_array[1],2)));
		
		for (int iX=0; iX<size_shadow_x; ++iX){
			for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
				Weight_PP_matrix[iX][iVx] = largest_distance-(sqrt(pow (Distance_X_array[iX],2) + pow(Distance_Vx_array[iVx],2)));
			}		
		}

		Wtotal = 0.0;
		for (int iX=0; iX<size_shadow_x; ++iX){
			for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
				Wtotal += Weight_PP_matrix[iX][iVx];
			}
		}
	
		for (int iX=0; iX<size_shadow_x; ++iX){
			for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
				Weight_PP_matrix[iX][iVx] = Weight_PP_matrix[iX][iVx]/Wtotal;
			}
		}

		i_X_array[0] = lX;
		int rX = lX;
		for (int iX=1; iX<size_shadow_x; ++iX){
			i_X_array[iX] = rX;
			++rX;
		}
		
		i_Vx_array[0] = lVx;
		int rVx = P->VxUp[lVx];
		for (int iVx=1; iVx<size_shadow_Vx; ++iVx){
			i_Vx_array[iVx] = rVx;
			rVx = P->VxUp[rVx];
		}
		for (int iX=0; iX<size_shadow_x; ++iX){
			for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
			T_nD_1D (i_Fg_matrix[iX][iVx],i_X_array[iX],i_Vx_array[iVx]);
			}
		}
// 		T_nD_1D (i_Fg_matrix[0][0],i_X_array[0],i_Vx_array[0]);
// 		T_nD_1D (i_Fg_matrix[0][1],i_X_array[0],i_Vx_array[1]);
// 		T_nD_1D (i_Fg_matrix[1][0],i_X_array[1],i_Vx_array[0]);
// 		T_nD_1D (i_Fg_matrix[1][1],i_X_array[1],i_Vx_array[1]);
		//       printf("(%i,%i)-->%i  ==== ",lX,lVx,i_Fg_ll);

	
		for (int iX=0; iX<size_shadow_x; ++iX){
			for (int iVx=0; iVx<size_shadow_Vx; ++iVx){
				DF_g_1D->at(i_Fg_matrix[iX][iVx])		+= Weight_PP_matrix[iX][iVx] * (*particle_iter).Data_PP->at(P->fp);
				CNT        [i_Fg_matrix[iX][iVx]]		+= Weight_PP_matrix[iX][iVx];
				#ifdef print_EPP
					EPP_g_1D->at(i_Fg_matrix[iX][iVx])		+=      (*particle_iter).Data_PP->at(P->epp);
					CNT_EPP [i_Fg_matrix[iX][iVx]]			+= 1.0;
				#endif
				#ifdef print_MPP
					MPP_g_1D->at(i_Fg_matrix[iX][iVx])		+=      (*particle_iter).Data_PP->at(P->mpp);
					CNT_MPP [i_Fg_matrix[iX][iVx]]			+= 1.0;
				#endif
			}
		}
	}
	//------------------------------------------------------------------------------
	
	//---------------------------------------------------------------------------------
	for (int iX = 0; iX<=P->nX_cpu-1; iX++) {
		for (int iVx = 0; iVx<=P->nVx-1; iVx++) {
			T_nD_1D (i_Fg_ii,iX,iVx);			
			//----------------------------------------------------------------------
			if (CNT[i_Fg_ii]!=0) {
				DF_g_1D->at(i_Fg_ii) 		 = DF_g_1D->at(i_Fg_ii)/CNT[i_Fg_ii];
				} else {
				//cout<<"zero grid point of"<<P->Name<<": iX="<<iX<<"  iVx="<<iVx<<endl;
				DF_g_1D->at(i_Fg_ii)		 = 0.0;    
			}
			//----------------------------------------------------------------------
			
			//----------------------------------------------------------------------
			#ifdef print_EPP
				if (CNT_EPP[i_Fg_ii]!=0) {
					EPP_g_1D->at(i_Fg_ii) 	 = EPP_g_1D->at(i_Fg_ii)/CNT_EPP[i_Fg_ii];
				} else {
					EPP_g_1D->at(i_Fg_ii)	 = 0.0;
				}
			#endif
			#ifdef print_MPP
				if (CNT_MPP[i_Fg_ii]!=0) {
					MPP_g_1D->at(i_Fg_ii) 	 = MPP_g_1D->at(i_Fg_ii)/CNT_MPP[i_Fg_ii];
				} else {
					MPP_g_1D->at(i_Fg_ii)	 = 0.0; 
				}
			#endif
			//----------------------------------------------------------------------
		}
	}
	for (int iVx = 0; iVx<=P->nVx-1; ++iVx) {
		T_nD_1D (i_Fg_ii,P->nX_cpu,iVx);
		DF_g_SND[iVx] = DF_g_1D->at(i_Fg_ii);
		#ifdef print_EPP
			EPP_g_SND[iVx] = EPP_g_1D->at(i_Fg_ii);
		#endif
		#ifdef print_MPP
			MPP_g_SND[iVx] = MPP_g_1D->at(i_Fg_ii);
		#endif
	}
  }
  //***********************************************************************************************

  //***********************************************************************************************
  void DF_Boundary_correction_for_iX_0 () {
	fill(DF_g_RCV.begin(), DF_g_RCV.end(), 0.0);
	#ifdef print_EPP
		fill(EPP_g_RCV.begin(), EPP_g_RCV.end(), 0.0);
	#endif
	#ifdef print_MPP
		fill(MPP_g_RCV.begin(), MPP_g_RCV.end(), 0.0);
	#endif
    //1) to send data of ghost cell to right cpu ------------------------------
	#if flag_timing_of_run==2
	struct timeval timer1,timer2;
	gettimeofday(&timer1, 0);
	#endif
    
     MPI_Sendrecv(&DF_g_SND[0], P->nVx+1,  MPI_DOUBLE, Right_cpu[rank], 233,
 		 &DF_g_RCV[0],  P->nVx+1,  MPI_DOUBLE, Left_cpu[rank],  233,
  		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     
	#ifdef print_EPP
		MPI_Sendrecv(&EPP_g_SND[0], P->nVx+1,  MPI_DOUBLE, Right_cpu[rank], 234,
 		 &EPP_g_RCV[0],  P->nVx+1,  MPI_DOUBLE, Left_cpu[rank],  234,
  		 MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
	#endif
		
	#ifdef print_MPP
		MPI_Sendrecv(&MPP_g_SND[0], P->nVx+1,  MPI_DOUBLE, Right_cpu[rank], 235,
			&MPP_g_RCV[0],  P->nVx+1,  MPI_DOUBLE, Left_cpu[rank],  235,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	#endif
	#if flag_timing_of_run==2
	gettimeofday(&timer2, 0);
	if (rank==0 && P->Charge == -1.0){ P->file_Info << " !"<< fixed <<(timer2.tv_sec - timer1.tv_sec)+(timer2.tv_usec - timer1.tv_usec)/1e6<<"! ";}
	#endif
    
    for (int iVx = 0; iVx<=P->nVx-1; iVx++) {
	
        //---------------------------------------------------------------------
        if (DF_g_1D->at(iVx)!=0.0 && DF_g_RCV[iVx]!=0.0){
	  DF_g_1D->at(iVx) = (DF_g_1D->at(iVx) + DF_g_RCV[iVx])/2.0;
	} else {
	  DF_g_1D->at(iVx) =  DF_g_1D->at(iVx) + DF_g_RCV[iVx];
	}
        //---------------------------------------------------------------------
	
	
	//---------------------------------------------------------------------
	#ifdef print_EPP
		if ( EPP_g_1D->at(iVx)!=0 && EPP_g_RCV[iVx] !=0){             
			EPP_g_1D->at(iVx) = (EPP_g_1D->at(iVx) + EPP_g_RCV[iVx])/2.0;
		} else {             
			EPP_g_1D->at(iVx) =  EPP_g_1D->at(iVx) + EPP_g_RCV[iVx];
		}
		
	#endif
	
	#ifdef print_MPP
		if ( MPP_g_1D->at(iVx)!=0 && MPP_g_RCV[iVx] !=0){             
			MPP_g_1D->at(iVx) = (MPP_g_1D->at(iVx) + MPP_g_RCV[iVx])/2.0;
		} else {             
			MPP_g_1D->at(iVx) =  MPP_g_1D->at(iVx) + MPP_g_RCV[iVx];
		}
	#endif
	//---------------------------------------------------------------------

      }
    //-------------------------------------------------------------------------
  }
  //***********************************************************************************************
 
 
   //***********************************************************************************************
  void DF_Boundary_correction_for_iX_nX () {
	int i_Fg_ii;
	int iX_0 = 0;
	fill(DF_g_RCV.begin(), DF_g_RCV.end(), 0.0);
	for (int iVx = 0; iVx<=P->nVx-1; ++iVx) {
		T_nD_1D (i_Fg_ii,iX_0,iVx);
		DF_g_SND[iVx] = DF_g_1D->at(i_Fg_ii);
		#ifdef print_EPP
			EPP_g_SND[iVx] = EPP_g_1D->at(i_Fg_ii);
		#endif
		#ifdef print_MPP
			MPP_g_SND[iVx] = MPP_g_1D->at(i_Fg_ii);
		#endif		
	}
	#ifdef print_EPP
		fill(EPP_g_RCV.begin(), EPP_g_RCV.end(), 0.0);
	#endif
	#ifdef print_MPP
		fill(MPP_g_RCV.begin(), MPP_g_RCV.end(), 0.0);
	#endif
    //1) to send data of ghost cell to right cpu ------------------------------
	#if flag_timing_of_run==2
	struct timeval timer1,timer2;
	gettimeofday(&timer1, 0);
	#endif
    
     MPI_Sendrecv(&DF_g_SND[0], P->nVx+1,  MPI_DOUBLE, Left_cpu[rank], 233,
 		 &DF_g_RCV[0],  P->nVx+1,  MPI_DOUBLE, Right_cpu[rank],  233,
  		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     
	#ifdef print_EPP
		MPI_Sendrecv(&EPP_g_SND[0], P->nVx+1,  MPI_DOUBLE, Left_cpu[rank], 234,
 		 &EPP_g_RCV[0],  P->nVx+1,  MPI_DOUBLE, Right_cpu[rank],  234,
  		 MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
	#endif
		
	#ifdef print_MPP
		MPI_Sendrecv(&MPP_g_SND[0], P->nVx+1,  MPI_DOUBLE, Left_cpu[rank], 235,
			&MPP_g_RCV[0],  P->nVx+1,  MPI_DOUBLE, Right_cpu[rank],  235,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	#endif
	#if flag_timing_of_run==2
	gettimeofday(&timer2, 0);
	if (rank==0 && P->Charge == -1.0){ P->file_Info << " !"<< fixed <<(timer2.tv_sec - timer1.tv_sec)+(timer2.tv_usec - timer1.tv_usec)/1e6<<"! ";}
 	#endif
    
    for (int iVx = 0; iVx<=P->nVx-1; iVx++) {
	T_nD_1D (i_Fg_ii,P->nX_cpu,iVx);
	DF_g_1D->at(i_Fg_ii) 		= DF_g_RCV[iVx];
       
	//---------------------------------------------------------------------
	#ifdef print_EPP
		EPP_g_1D->at(i_Fg_ii) 	=   EPP_g_RCV[iVx];	
	#endif
	
	#ifdef print_MPP
		MPP_g_1D->at(iVx) 	=  MPP_g_RCV[iVx];
	#endif
	//---------------------------------------------------------------------

      }
    //-------------------------------------------------------------------------
  }
  //***********************************************************************************************
  
  //***********************************************************************************************
  void Num_Density_Calculator_zeor_order ()  {
      int i_Fg_ii;
      fill(Num_Den_g->begin(), Num_Den_g->end(), 0.0);
      for (int iX = 0; iX<=P->nX_cpu-1; iX++) {
	  for (int iVx = 0; iVx<=P->nVx-1; iVx++) {
	    T_nD_1D (i_Fg_ii,iX,iVx);

 		Num_Den_g->at(iX) = Num_Den_g->at(iX) + DF_g_1D->at(i_Fg_ii);  
	  }   
	}
      for (int iX = 0; iX<=P->nX_cpu-1; iX++) {
	Num_Den_g->at(iX) =  P->dVx * Num_Den_g->at(iX);
      }
  }
  //***********************************************************************************************
  
  //***********************************************************************************************
  void Num_Density_Calculator ()  {
      int i_Fg_ii;
      fill(Num_Den_g->begin(), Num_Den_g->end(), 0.0);
      for (int iX = 0; iX<=P->nX_cpu-1; iX++) {
	  for (int iVx = 0; iVx<=P->nVx-1; iVx++) {
	    T_nD_1D (i_Fg_ii,iX,iVx);
	    if (iVx%2 == 0.0) 
 		Num_Den_g->at(iX) = Num_Den_g->at(iX) + 2.0 * DF_g_1D->at(i_Fg_ii);  
	      else 
 		Num_Den_g->at(iX) = Num_Den_g->at(iX) + 4.0 * DF_g_1D->at(i_Fg_ii);   
	  }   
	}
      for (int iX = 0; iX<=P->nX_cpu-1; iX++) {
	Num_Den_g->at(iX) = (1.0/3.0) * P->dVx * Num_Den_g->at(iX);
      }
  }
  //***********************************************************************************************
 
  //***********************************************************************************************
  void Write_Dset_Num_Den_HDF5 (hid_t & Group_ID_S)  {
      stringstream DSet_name;
      DSet_name << "Num_Den_"<<P->Name;
    //-------------------------------------------------------------------------
   
    // 2) Create the dataspace for the dataset  -------------------------------
    P->Space_ID_Num_Den = H5Screate_simple(P->Rank_X, P->Dim_X , NULL); 
    P->Memory_ID_Num_Den = H5Screate_simple(P->Rank_X, P->Dim_X_Chunk, NULL);
    //-------------------------------------------------------------------------
   
    // 3) Create chunked dataset ----------------------------------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(P->PList_ID, P->Rank_X, P->Dim_X_Chunk);

    P->DSet_ID_Num_Den = H5Dcreate(Group_ID_S, DSet_name.str().c_str(), H5T_NATIVE_DOUBLE, P->Space_ID_Num_Den,
			  H5P_DEFAULT, P->PList_ID, H5P_DEFAULT);
    H5Pclose(P->PList_ID);
    //------------------------------------------------------------------------
    H5Sclose(P->Space_ID_Num_Den);
  
    // 4) Select hyperslab in the file----------------------------------------
    P->Space_ID_Num_Den = H5Dget_space(P->DSet_ID_Num_Den);
    P->status_hdf5 = H5Sselect_hyperslab(P->Space_ID_Num_Den, H5S_SELECT_SET, P->OffSet_X, P->Stride_X, P->Count_X, P->Block_X);
    //-----------------------------------------------------------------------

    // 5) Create property list for collective dataset write------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(P->PList_ID, H5FD_MPIO_COLLECTIVE);
    //-----------------------------------------------------------------------
    
    P->status_hdf5 = H5Dwrite(P->DSet_ID_Num_Den, H5T_NATIVE_DOUBLE, P->Memory_ID_Num_Den, P->Space_ID_Num_Den,
		      P->PList_ID, Num_Den_g->data());
   
    //H5Fflush(File_ID_S, H5F_SCOPE_GLOBAL ); 
    
    // Close/release resources    
    H5Dclose(P->DSet_ID_Num_Den);
    H5Pclose(P->PList_ID);
    
    H5Sclose(P->Space_ID_Num_Den);
    H5Sclose(P->Memory_ID_Num_Den);
  }
  //***********************************************************************************************
 
 
 //***********************************************************************************************
  void Write_Dset_Ene_kin_HDF5 (hid_t & Group_ID_S)  {
      stringstream DSet_name;
      DSet_name << "Ene_Kin_"<<P->Name;
    //-------------------------------------------------------------------------
   
    // 2) Create the dataspace for the dataset  -------------------------------
    P->Space_ID_Ene_kin = H5Screate_simple(P->Rank_X, P->Dim_X , NULL); 
    P->Memory_ID_Ene_kin = H5Screate_simple(P->Rank_X, P->Dim_X_Chunk, NULL);
    //-------------------------------------------------------------------------
   
    // 3) Create chunked dataset ----------------------------------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(P->PList_ID, P->Rank_X, P->Dim_X_Chunk);

    P->DSet_ID_Ene_kin = H5Dcreate(Group_ID_S, DSet_name.str().c_str(), H5T_NATIVE_DOUBLE, P->Space_ID_Ene_kin,
			  H5P_DEFAULT, P->PList_ID, H5P_DEFAULT);
    H5Pclose(P->PList_ID);
    //------------------------------------------------------------------------
    H5Sclose(P->Space_ID_Ene_kin);
  
    // 4) Select hyperslab in the file----------------------------------------
    P->Space_ID_Ene_kin = H5Dget_space(P->DSet_ID_Ene_kin);
    P->status_hdf5 = H5Sselect_hyperslab(P->Space_ID_Ene_kin, H5S_SELECT_SET, P->OffSet_X, P->Stride_X, P->Count_X, P->Block_X);
    //-----------------------------------------------------------------------

    // 5) Create property list for collective dataset write------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(P->PList_ID, H5FD_MPIO_COLLECTIVE);
    //-----------------------------------------------------------------------
    
    P->status_hdf5 = H5Dwrite(P->DSet_ID_Ene_kin, H5T_NATIVE_DOUBLE, P->Memory_ID_Ene_kin, P->Space_ID_Ene_kin,
		      P->PList_ID, Ene_kin_x_cpu->data());
   
    //H5Fflush(File_ID_S, H5F_SCOPE_GLOBAL ); 
    
    // Close/release resources    
    H5Dclose(P->DSet_ID_Ene_kin);
    H5Pclose(P->PList_ID);
    
    H5Sclose(P->Space_ID_Ene_kin);
    H5Sclose(P->Memory_ID_Ene_kin);
  }
  //***********************************************************************************************

 
 
	//***********************************************************************************************
	void Write_Dset_Contour_HDF5(hid_t & Group_ID_DF) {
		
		/*
		vector<double> * DF_g_1D_print;
		double Vx_min_print, Vx_max_print,Vx0;
		int nVx_print, i_counter,i_Fg_print;
		Vx_min_print = -2.0;
		Vx_max_print = 2.0;
		
		Vx_min_print = Vx_min_print*sqrt(1.0/P->norm_factor);
		Vx_max_print = Vx_max_print*sqrt(1.0/P->norm_factor);
		nVx_print = (Vx_max_print - Vx_min_print)/P->nVx;
		DF_g_1D_print = new vector<double> (P->nX*nVx_print);
		//==================================================================================
		i_counter = 0;
		for (int iX = 0; iX<=P->nX_cpu-1; iX++)  {
			for (int iVx = 0; iVx<=P->nVx-1; iVx++) {
				Vx0 = P->Vx_min + double(iVx) * P->dVx;
				if (Vx0>Vx_min_print && Vx0<Vx_max_print){
					T_nD_1D (i_Fg_print,iX,iVx);
					DF_g_1D_print->at(i_counter)=DF_g_1D->at(i_Fg_print);
					i_counter++;
				}
				
			}
		}
		*/
		//-------------------------------------------------------------------------
		//     #ifdef debug
		P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="********  Write Contour HDF5 ");
		//     #endif
		// 2) Create the dataspace for the dataset  -------------------------------
		P->Space_ID_DF = H5Screate_simple(P->Rank_X_Vx, P->Dim_X_Vx , NULL); 
		P->Memory_ID_DF = H5Screate_simple(P->Rank_X_Vx, P->Dim_X_Vx_Chunk, NULL);
		//-------------------------------------------------------------------------

		// 3) Create chunked dataset ----------------------------------------------
		P->PList_ID = H5Pcreate(H5P_DATASET_CREATE);
		H5Pset_chunk(P->PList_ID, P->Rank_X_Vx, P->Dim_X_Vx_Chunk);

		stringstream DSet_name;
		DSet_name << "DF_"<<P->Name; 

		P->DSet_ID_DF = H5Dcreate(Group_ID_DF, DSet_name.str().c_str(), H5T_NATIVE_DOUBLE, P->Space_ID_DF,
			  H5P_DEFAULT, P->PList_ID, H5P_DEFAULT);
		H5Pclose(P->PList_ID);
		//------------------------------------------------------------------------
		H5Sclose(P->Space_ID_DF);

		// 4) Select hyperslab in the file----------------------------------------
		P->Space_ID_DF = H5Dget_space(P->DSet_ID_DF);
		P->status_hdf5 = H5Sselect_hyperslab(P->Space_ID_DF, H5S_SELECT_SET, P->OffSet_X_Vx, P->Stride_X_Vx, P->Count_X_Vx, P->Block_X_Vx);
		//-----------------------------------------------------------------------

		// 5) Create property list for collective dataset write------------------
		P->PList_ID = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(P->PList_ID, H5FD_MPIO_COLLECTIVE);
		//-----------------------------------------------------------------------


		P->status_hdf5 = H5Dwrite(P->DSet_ID_DF, H5T_NATIVE_DOUBLE, P->Memory_ID_DF, P->Space_ID_DF,
		      P->PList_ID, DF_g_1D->data());

		//H5Fflush(File_ID_S, H5F_SCOPE_GLOBAL ); 

		// Close/release resources    
		H5Dclose(P->DSet_ID_DF);
		H5Pclose(P->PList_ID);

		H5Sclose(P->Space_ID_DF);
		H5Sclose(P->Memory_ID_DF);
	}
	//***********************************************************************************************

    //***********************************************************************************************
  void Write_Dset_EPP_HDF5(hid_t & Group_ID_DF) {
      //-------------------------------------------------------------------------
//     #ifdef debug
	P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="********  Write Contour HDF5 ");
//     #endif
    // 2) Create the dataspace for the dataset  -------------------------------
    P->Space_ID_DF = H5Screate_simple(P->Rank_X_Vx, P->Dim_X_Vx , NULL); 
    P->Memory_ID_DF = H5Screate_simple(P->Rank_X_Vx, P->Dim_X_Vx_Chunk, NULL);
    //-------------------------------------------------------------------------
   
    // 3) Create chunked dataset ----------------------------------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(P->PList_ID, P->Rank_X_Vx, P->Dim_X_Vx_Chunk);
    
    stringstream DSet_name;
    DSet_name << "EPP_"<<P->Name; 

    P->DSet_ID_DF = H5Dcreate(Group_ID_DF, DSet_name.str().c_str(), H5T_NATIVE_DOUBLE, P->Space_ID_DF,
			  H5P_DEFAULT, P->PList_ID, H5P_DEFAULT);
    H5Pclose(P->PList_ID);
    //------------------------------------------------------------------------
    H5Sclose(P->Space_ID_DF);
  
    // 4) Select hyperslab in the file----------------------------------------
    P->Space_ID_DF = H5Dget_space(P->DSet_ID_DF);
    P->status_hdf5 = H5Sselect_hyperslab(P->Space_ID_DF, H5S_SELECT_SET, P->OffSet_X_Vx, P->Stride_X_Vx, P->Count_X_Vx, P->Block_X_Vx);
    //-----------------------------------------------------------------------

    // 5) Create property list for collective dataset write------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(P->PList_ID, H5FD_MPIO_COLLECTIVE);
    //-----------------------------------------------------------------------

  
    
	#ifdef print_EPP
		P->status_hdf5 = H5Dwrite(P->DSet_ID_DF, H5T_NATIVE_DOUBLE, P->Memory_ID_DF, P->Space_ID_DF,
		      P->PList_ID, EPP_g_1D->data());
	#endif
   
    //H5Fflush(File_ID_S, H5F_SCOPE_GLOBAL ); 
    
    // Close/release resources    
    H5Dclose(P->DSet_ID_DF);
    H5Pclose(P->PList_ID);
    
    H5Sclose(P->Space_ID_DF);
    H5Sclose(P->Memory_ID_DF);
    }
  //***********************************************************************************************
  
     //***********************************************************************************************
  void Write_Dset_MPP_HDF5(hid_t & Group_ID_DF) {
      //-------------------------------------------------------------------------
//     #ifdef debug
	P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="********  Write Contour HDF5 ");
//     #endif
    // 2) Create the dataspace for the dataset  -------------------------------
    P->Space_ID_DF = H5Screate_simple(P->Rank_X_Vx, P->Dim_X_Vx , NULL); 
    P->Memory_ID_DF = H5Screate_simple(P->Rank_X_Vx, P->Dim_X_Vx_Chunk, NULL);
    //-------------------------------------------------------------------------
   
    // 3) Create chunked dataset ----------------------------------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(P->PList_ID, P->Rank_X_Vx, P->Dim_X_Vx_Chunk);
    
    stringstream DSet_name;
    DSet_name << "MPP_"<<P->Name; 

    P->DSet_ID_DF = H5Dcreate(Group_ID_DF, DSet_name.str().c_str(), H5T_NATIVE_DOUBLE, P->Space_ID_DF,
			  H5P_DEFAULT, P->PList_ID, H5P_DEFAULT);
    H5Pclose(P->PList_ID);
    //------------------------------------------------------------------------
    H5Sclose(P->Space_ID_DF);
  
    // 4) Select hyperslab in the file----------------------------------------
    P->Space_ID_DF = H5Dget_space(P->DSet_ID_DF);
    P->status_hdf5 = H5Sselect_hyperslab(P->Space_ID_DF, H5S_SELECT_SET, P->OffSet_X_Vx, P->Stride_X_Vx, P->Count_X_Vx, P->Block_X_Vx);
    //-----------------------------------------------------------------------

    // 5) Create property list for collective dataset write------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(P->PList_ID, H5FD_MPIO_COLLECTIVE);
    //-----------------------------------------------------------------------

  
    
	#ifdef print_MPP
		P->status_hdf5 = H5Dwrite(P->DSet_ID_DF, H5T_NATIVE_DOUBLE, P->Memory_ID_DF, P->Space_ID_DF,
		      P->PList_ID, MPP_g_1D->data());
	#endif
   
    //H5Fflush(File_ID_S, H5F_SCOPE_GLOBAL ); 
    
    // Close/release resources    
    H5Dclose(P->DSet_ID_DF);
    H5Pclose(P->PList_ID);
    
    H5Sclose(P->Space_ID_DF);
    H5Sclose(P->Memory_ID_DF);
    }
  //***********************************************************************************************

  
  //***********************************************************************************************
	void Moving_box_simulation (ParametersClass * P, const double& d_Time) {
		// moving frame
		P->X_max_Tot += P->velocity_box_simulation * d_Time;
		P->X_min_Tot += P->velocity_box_simulation * d_Time;
		P->X_min     += P->velocity_box_simulation * d_Time;
		P->X_max     += P->velocity_box_simulation * d_Time;
	}
 //***********************************************************************************************
 
 //***********************************************************************************************
  void RELOCATION_CPU(){
  
    //---------------------------------------------------------------
    double X_middle=  ((P->X_max + P->X_min)/2) + ((P->X_max_Tot - P->X_min_Tot)/2);
    
	//---------------------------------------------------------------
	/*int Number_of_data_in_PP, i_number_of_data_EPP, i_number_of_data_MPP;
	Number_of_data_in_PP = 3;
	i_number_of_data_MPP = 0;
	i_number_of_data_EPP = 0;
    
    
	#ifdef print_EPP
		Number_of_data_in_PP += 1;
		i_number_of_data_EPP = 3;
	#endif
	#ifdef print_MPP	
		Number_of_data_in_PP += 1;
		i_number_of_data_MPP = 3;
	#endif
	//this has to be after ifdef print_MPP, and yes EPP flag intefers with MPP, just to make sure if the two flags are on the numbering of data goes smoothly
	#ifdef print_EPP
		i_number_of_data_MPP = 4;
	#endif
        */
	//---------------------------------------------------------------
	
    if (X_middle>P->X_max_Tot) { //periodic boundary condition
      X_middle = P->X_min_Tot + (X_middle - P->X_max_Tot);
    }
    if (X_middle<P->X_min_Tot) { 
      X_middle=P->X_max_Tot + (X_middle-P->X_min_Tot);
    }
    
    bool X_middle_right=false, X_middle_left=false;
    if (X_middle > P->X_max) X_middle_right = true;
    if (X_middle < P->X_min)  X_middle_left = true;
    //---------------------------------------------------------------

    //     printf("start of relocation: rank %i, X_max %f,X_middle %f, X_middle_left %d, X_middle_right %d \n",
    //	   rank, P->X_max,X_middle,X_middle_left, X_middle_right);
    
    
    //---------------------------------------------------------------
    vector<double> * Data_snd_right, * Data_snd_left;
    vector<double> * Data_rcv_right, * Data_rcv_left;
    vector<double> Data_store_PP;
    int N_snd_right, N_rcv_left;
    int N_snd_left, N_rcv_right;
    
    int i, i_count;
    double X0, Vx0, Fp0, EPP0, MPP0;
    EPP0 = 0;
    MPP0 = 0;
    vector<int> Flag_Transfer_cpu_array(size);
    int Flag_Transfer_Tot, Flag_Transfer_cpu;
    
	#if flag_timing_of_run==2
	struct timeval timer1,timer2;
	#endif
    double Energy_kinetic;
    
    
    i_count = 0;
    Flag_Transfer_cpu = 1;
      
	do {
		#ifdef debug
		if (rank==0) {cout << "in relocation cpu =" << endl;}
		#endif
	
	/* 0) initializing the loop
		1) checking the phase points
		0) recognizing them
		1) adding their data to the buffer
		2) delete them from the storage
		2)sendrecv the data
		3)communication flag for cpus
		4) deleting the initial vectors*/
	
	//---------------------------------------------------------------
	Data_snd_right = new vector<double>;
	Data_snd_left  = new vector<double>;
	Data_store_PP.resize(P->number_of_data);
	
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("start of relocation iteration: cpu_flag %i, total flag %i, rank %i \n"
	//,Flag_Transfer_cpu, Flag_Transfer_Tot ,rank);

	
	i_count=i_count+1;
	
	N_snd_right = 0;
	N_rcv_left = 0;
	N_snd_left=0;
	N_rcv_right=0;     
	//---------------------------------------------------------------

	//---------------------------------------------------------------
	if (Flag_Transfer_cpu == 1) {       
		if ( X_middle_right == true) {
		for(Storage<PhasePointClass>::Iterator iter = PhasePoints->begin(); iter != PhasePoints->end(); ++iter) { 
		#ifdef debug
			if ( (*iter).Data_PP->at(P->x) >= P->X_max_Tot || (*iter).Data_PP->at(P->x) <= P->X_min_Tot  ) {
			cout << "out of boundary in X (relocation cpu) =" << (*iter).Data_PP->at(P->x) << "-- for ="<< P->Name <<" ---rank="<<rank<< endl;
		}
		#endif
		if ( (*iter).Data_PP->at(P->x) >= P->X_max || (*iter).Data_PP->at(P->x) < P->X_min  ) { //jenab <= is replaced by <, 2018-11-16
		if ((*iter).Data_PP->at(P->x) >= P->X_max &&  (*iter).Data_PP->at(P->x) <= X_middle ) {		
				for (int i_data=0; i_data<P->number_of_data; ++i_data){
					Data_snd_right->push_back((*iter).Data_PP->at(i_data));
					++N_snd_right;
				}
	// 			Data_snd_right->push_back((*iter).Data_PP->at(P->x));
	// 			++N_snd_right;
	// 			
	// 			Data_snd_right->push_back((*iter).Data_PP->at(P->vx));
	// 			++N_snd_right;
	// 			
	// 			Data_snd_right->push_back((*iter).Data_PP->at(P->fp));
	// 			++N_snd_right;
	// 			#ifdef print_EPP
	// 				Data_snd_right->push_back((*iter).EPP);
	// 				++N_snd_right;
	// 			#endif
	// 			
	// 			#ifdef print_MPP
	// 				Data_snd_right->push_back((*iter).MPP);
	// 				++N_snd_right;
	// 			#endif
			
	//  		PhasePointClass PhasePoint_delete = iter.pop();
			(*iter).delete_Data_PP();
			iter.pop();
		
		} else {
				for (int i_data=0; i_data<P->number_of_data; ++i_data){
					Data_snd_left->push_back((*iter).Data_PP->at(i_data));
					++N_snd_left;
				}
	// 			Data_snd_left->push_back((*iter).Data_PP->at(P->x));    
	// 			++N_snd_left;
	// 		
	// 			Data_snd_left->push_back((*iter).Data_PP->at(P->vx));	    
	// 			++N_snd_left;
	// 		
	// 			Data_snd_left->push_back((*iter).Data_PP->at(P->fp));
	// 			++N_snd_left;
	// 			#ifdef print_EPP
	// 				Data_snd_left->push_back((*iter).EPP);
	// 				++N_snd_left;
	// 			#endif
	// 			
	// 			#ifdef print_MPP
	// 				Data_snd_left->push_back((*iter).MPP);
	// 				++N_snd_left;
	// 			#endif
			
		
	//  		PhasePointClass PhasePoint_delete = iter.pop();
			(*iter).delete_Data_PP();
			iter.pop();
		}
		}
		}
		} 
		else if ( X_middle_left == true) {  
		for(Storage<PhasePointClass>::Iterator iter = PhasePoints->begin(); iter != PhasePoints->end(); ++iter) { 
		#ifdef debug
			if ( (*iter).Data_PP->at(P->x) >= P->X_max_Tot || (*iter).Data_PP->at(P->x) <= P->X_min_Tot  ) {
			cout << "out of boundary in X (relocation cpu) =" << (*iter).Data_PP->at(P->x) << "-- for ="<< P->Name <<" ---rank="<<rank<< endl;
		}
		#endif
		if ( (*iter).Data_PP->at(P->x) >= P->X_max || (*iter).Data_PP->at(P->x) < P->X_min  ) { //jenab <= is replaced by <, 2018-11-16
		if ((*iter).Data_PP->at(P->x) <= P->X_min &&  (*iter).Data_PP->at(P->x) >= X_middle ) { 
				for (int i_data=0; i_data<P->number_of_data; ++i_data){
					Data_snd_left->push_back((*iter).Data_PP->at(i_data));
					++N_snd_left;
				}
	// 			Data_snd_left->push_back((*iter).Data_PP->at(P->x));    
	// 			++N_snd_left;
	// 		
	// 			Data_snd_left->push_back((*iter).Data_PP->at(P->vx));	    
	// 			++N_snd_left;
	// 		
	// 			Data_snd_left->push_back((*iter).Data_PP->at(P->fp));
	// 			++N_snd_left;
	// 			#ifdef print_EPP
	// 			Data_snd_left->push_back((*iter).EPP);
	// 			++N_snd_left;
	// 			#endif
	// 			
	// 			#ifdef print_MPP
	// 			Data_snd_left->push_back((*iter).MPP);
	// 			++N_snd_left;
	// 			#endif
		
	//  		PhasePointClass PhasePoint_delete = iter.pop();
			(*iter).delete_Data_PP();
			iter.pop();
			
		} else {
				for (int i_data=0; i_data<P->number_of_data; ++i_data){
					Data_snd_right->push_back((*iter).Data_PP->at(i_data));
					++N_snd_right;
				}
		
	// 			Data_snd_right->push_back((*iter).Data_PP->at(P->x));
	// 			++N_snd_right;
	// 			
	// 			Data_snd_right->push_back((*iter).Data_PP->at(P->vx));
	// 			++N_snd_right;
	// 			
	// 			Data_snd_right->push_back((*iter).Data_PP->at(P->fp));
	// 			++N_snd_right;
	// 			
	// 			#ifdef print_EPP
	// 			Data_snd_right->push_back((*iter).EPP);
	// 			++N_snd_right;
	// 			#endif
	// 			
	// 			#ifdef print_MPP
	// 			Data_snd_right->push_back((*iter).MPP);
	// 			++N_snd_right;
	// 			#endif
			
	//  		PhasePointClass PhasePoint_delete = iter.pop();
			(*iter).delete_Data_PP();
			iter.pop();
		}
		}
		}
		}
	}
	//---------------------------------------------------------------

	//---------------------------------------------------------------
	// can be optmized by having even/odd send-and-recieve
	
	#if flag_timing_of_run==2
	gettimeofday(&timer1, 0);
	#endif
	
		MPI_Sendrecv(&N_snd_right, 1, MPI_INT, Right_cpu[rank], 144,
			&N_rcv_left,  1, MPI_INT, Left_cpu[rank],  144,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if (N_rcv_left!=0 ) {
		Data_rcv_left  = new vector<double>(N_rcv_left);
	}else {
		Data_rcv_left  = new vector<double>(1);;      
	}
	
	if (N_snd_right == 0) Data_snd_right->push_back(0);
	//if (N_rcv_left!=0) Data_rcv_left->assign(Data_rcv_left->size(),0);
	//if (N_rcv_left!=0) fill(Data_rcv_left->begin(), Data_rcv_left->end(), 0.0);
		MPI_Sendrecv(&Data_snd_right->at(0), N_snd_right, MPI_DOUBLE, Right_cpu[rank], 233,
			&Data_rcv_left->at(0),  N_rcv_left,  MPI_DOUBLE, Left_cpu[rank],  233,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	
		if (N_rcv_left!=0) {
			for(int i_num_pp=0;  i_num_pp<N_rcv_left-1; i_num_pp=i_num_pp+P->number_of_data) {
				for (int i_data=0; i_data<P->number_of_data; ++i_data){	
					Data_store_PP[i_data]= Data_rcv_left->at(i_num_pp +i_data);
				}
				PhasePointClass PhasePoint (Data_store_PP, P->number_of_data);
				PhasePoints->push_back(PhasePoint);
			}
		
	// 		//check if i should start from zero
	// 		for((i=0); (i<N_rcv_left-1); (i=i+Number_of_data_in_PP)) { 
	// 			if (rank==0 && i>N_rcv_left-10)cout <<"  i="<<i<<"  "<<N_rcv_left<<endl;
	// 			X0  = Data_rcv_left->at(i);
	// 			Vx0 = Data_rcv_left->at(i+1);
	// 			Fp0 = Data_rcv_left->at(i+2);
	// 			#ifdef print_EPP			
	// 			EPP0 = Data_rcv_left->at(i+i_number_of_data_EPP);
	// 			#endif
	// 			#ifdef print_MPP
	// 			MPP0 = Data_rcv_left->at(i+i_number_of_data_MPP);
	// 			#endif						
	// 			PhasePointClass PhasePoint (X0,Vx0,Fp0,EPP0,MPP0);
	// 			PhasePoints->push_back(PhasePoint);
	// 		}
		}
	//---------------------------------------------------------------
		MPI_Sendrecv(&N_snd_left,  1, MPI_INT, Left_cpu[rank],  377,
			&N_rcv_right, 1, MPI_INT, Right_cpu[rank], 377,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if (N_rcv_right!=0) {
		Data_rcv_right = new vector<double>(N_rcv_right);
	}else{
		Data_rcv_right = new vector<double>(1);
	}
	if (N_snd_left == 0) Data_snd_left->push_back(0);
		MPI_Sendrecv(&Data_snd_left->at(0),  N_snd_left,  MPI_DOUBLE, Left_cpu[rank],  610,
			&Data_rcv_right->at(0), N_rcv_right, MPI_DOUBLE, Right_cpu[rank], 610,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	#if flag_timing_of_run==2
	gettimeofday(&timer2, 0);
	if (rank==0 && P->Charge == -1.0){ P->file_Info << " !"<< fixed <<(timer2.tv_sec - timer1.tv_sec)+(timer2.tv_usec - timer1.tv_usec)/1e6<<"! ";}
	#endif
		if (N_rcv_right!=0) {
			for(int i_num_pp=0;  i_num_pp<N_rcv_right-1; i_num_pp=i_num_pp+P->number_of_data) { 
				for (int i_data=0; i_data<P->number_of_data; ++i_data){	
					Data_store_PP[i_data]= Data_rcv_right->at(i_num_pp +i_data);
				}
				PhasePointClass PhasePoint (Data_store_PP, P->number_of_data);
				PhasePoints->push_back(PhasePoint);
			}
	// 		for((i=0); (i<N_rcv_right-1); (i=i+Number_of_data_in_PP)) { //check if i should start from zero
	// 		X0  = Data_rcv_right->at(i);
	// 		Vx0 = Data_rcv_right->at(i+1);
	// 		Fp0 = Data_rcv_right->at(i+2);
	// 		#ifdef print_EPP
	// 			EPP0 = Data_rcv_right->at(i+i_number_of_data_EPP);
	// 		#endif
	// 		#ifdef print_MPP
	// 			MPP0 = Data_rcv_right->at(i+i_number_of_data_MPP);
	// 		#endif
	// 		
	// 		PhasePointClass PhasePoint (X0,Vx0,Fp0,EPP0,MPP0);
	// 		PhasePoints->push_back(PhasePoint);
	// 		}
		}
	//---------------------------------------------------------------
	/*printf("Time%f: ,  %i <---%i---  %i ---%i---> %i |||   %i ---%i--->  %i <---%i--- %i \n",
		P->R_Time, Left_cpu[rank],N_snd_left, rank,N_snd_right,Right_cpu[rank],
		Left_cpu[rank],N_rcv_left, rank,N_rcv_right,Right_cpu[rank]);
	vector<int> N_snd_right_array(size), N_snd_left_array(size), N_rcv_left_array(size), N_rcv_right_array(size) ;
	MPI_Gather(&N_snd_right, 1, MPI_INT,&N_snd_right_array.at(0), 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&N_snd_left, 1, MPI_INT,&N_snd_left_array.at(0), 1, MPI_INT, 0, MPI_COMM_WORLD);
		
	MPI_Gather(&N_rcv_right, 1, MPI_INT,&N_rcv_right_array.at(0), 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&N_rcv_left, 1, MPI_INT,&N_rcv_left_array.at(0), 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (rank==0) {
		int leaking=0;
		
		for (i=0; (i<=NODE); ++i){
		
		leaking = leaking
			+ N_snd_right_array.at(i)
			+ N_snd_left_array.at(i)
			- N_rcv_right_array.at(i)
			- N_rcv_left_array.at(i);
		
		}
		cout <<" leaking= "<<leaking<<endl;
	}*/
	
	
	//-------------------------------- set the flag for each cpu
	Flag_Transfer_cpu = 0;
	Flag_Transfer_Tot = 0;
	for(Storage<PhasePointClass>::Iterator iter = PhasePoints->begin(); iter != PhasePoints->end(); ++iter) { 
		if ( ((*iter).Data_PP->at(P->x) >= P->X_max) | ((*iter).Data_PP->at(P->x) <= P->X_min) ) {
		Flag_Transfer_cpu = 1;
		break; // this make sure by the first instnace of phase point out of the domain, the loop will be stopped.
		}
	}
	//--------------- gather, sum up, broadcast the flags -------------------------------------------------

	MPI_Gather(&Flag_Transfer_cpu, 1, MPI_INT,
			&Flag_Transfer_cpu_array.at(0), 1, MPI_INT,
			0, MPI_COMM_WORLD);
	if (rank==0) {
		Flag_Transfer_Tot= 0;
		for (i=0; (i<=NODE); ++i){
		Flag_Transfer_Tot += Flag_Transfer_cpu_array.at(i);
		//cout <<Flag_Transfer_Tot<<endl;
		}
	}
		MPI_Bcast( &Flag_Transfer_Tot, 1,  MPI_INT, 0, MPI_COMM_WORLD);
	//           P->Gather_Compare_Cast_INT(Flag_Transfer_cpu, Flag_Transfer_Tot,1,P->Debug_mes="--- (0): relocation cpu");
	
	//----------------------------------------------------------------------------------------------------

	delete Data_snd_right;
	delete Data_snd_left;
	delete Data_rcv_right;
	delete Data_rcv_left;
	} 
	while (Flag_Transfer_Tot >= 1);
	}  //  RELOCATION_CPU()
  //***********************************************************************************************  
  
 
  //***********************************************************************************************
  void Entropy_F2_Species(){
    Entropy_F2_cpu (S_F2_cpu);
    if (size>1){
        P->Sum_up_Double(S_F2_cpu,S_F2_Tot,1);
    }else{
        S_F2_Tot = S_F2_cpu;
    }    
  }    
  //***********************************************************************************************
  
  //***********************************************************************************************
  void Entropy_F2_cpu(double & S_F2_cpu) {
    vector<double> S_F2_temp(P->nX_cpu+1);
    fill(S_F2_temp.begin(), S_F2_temp.end(), 0.0);

    //Integration over v
    int i_Fg_ii,d_iVx,u_iVx,i_Fg_i0,i_Fg_id,i_Fg_iu;
    for(int iX=0; iX<=P->nX_cpu-1; iX++){
      for(int iVx=1; iVx<=P->nVx-2; iVx++){ //jenab : problem starting integration from 1  
	d_iVx = iVx-1;
	u_iVx = iVx+1;
	T_nD_1D (i_Fg_ii,iX,iVx);
	T_nD_1D (i_Fg_id,iX,d_iVx);
	T_nD_1D (i_Fg_iu,iX,u_iVx);
	
	S_F2_temp[iX] = S_F2_temp[iX] + (
		      (pow(DF_g_1D->at(i_Fg_id),2) + 4 * pow(DF_g_1D->at(i_Fg_ii),2) + pow(DF_g_1D->at(i_Fg_iu),2))
					/3)*P->dVx;
	  }
 	  int iVx0 = 0;
 	  T_nD_1D (i_Fg_i0,iX,iVx0);
	  S_F2_temp[iX] = S_F2_temp[iX] + (pow(DF_g_1D->at(i_Fg_i0),2)) * P->dVx;
	}

	//Integration over x
	S_F2_cpu = 0.0;
// 	double ghost_minus_2,ghost_minus_1,ghost_nX_cpu,ghost_nX_cpu_plus_1;
	P->Ghost_Grid_Points_Vector(S_F2_temp,P->nX_cpu,1);
	S_F2_cpu= S_F2_cpu + ( 
		      (P->minus_1+ 4 * S_F2_temp[0]+ S_F2_temp[1])
			  /3 )*P->dX;//CPU Boundary Correction
	
	S_F2_cpu= S_F2_cpu + ( 
		      (S_F2_temp[P->nX_cpu-1] + 4 *P->ghost_nX_cpu + P->nX_cpu_plus_1)
			  /3 )*P->dX;//CPU Boundary Correction
			  
			  
			  
	for(int iX=1; iX<=P->nX_cpu-2; iX++){  
	  S_F2_cpu= S_F2_cpu + ( 
			  (S_F2_temp[iX-1]+ 4 * S_F2_temp[iX]+ S_F2_temp[iX+1])
			  /3 )*P->dX;
	}
	 //-------------------------------------------------------------------
    
    }
  //***********************************************************************************************
      
  //***********************************************************************************************
  void Energy_kin_Species(){
    Energy_kin_cpu (Ene_kin_cpu);
    if (size>1){
         P->Sum_up_Double(Ene_kin_cpu,Ene_kin_Tot,1);
    }else{
        Ene_kin_Tot = Ene_kin_cpu;
    }   
  }    
  //***********************************************************************************************
  
  //***********************************************************************************************
  void Entropy_FlnF_Species(){      
    Entropy_FlnF_cpu (S_FlnF_cpu);
    if (size>1){
         P->Sum_up_Double(S_FlnF_cpu,S_FlnF_Tot,1); 
    }else{
        S_FlnF_Tot = S_FlnF_cpu;
    }     
  }    
  //***********************************************************************************************
  
  //***********************************************************************************************
  void Entropy_FlnF_cpu(double & S_FlnF_cpu) {
    vector<double> S_FlnF_temp(P->nX_cpu+1);
    fill(S_FlnF_temp.begin(), S_FlnF_temp.end(), 0);

    //Integration over v
    int i_Fg_ii,d_iVx,u_iVx,i_Fg_i0,i_Fg_id,i_Fg_iu;
    for(int iX=0; iX<=P->nX_cpu-1; iX++){
      for(int iVx=1; iVx<=P->nVx-2; iVx++){
	d_iVx = iVx-1;
	u_iVx = iVx+1;
	T_nD_1D (i_Fg_ii,iX,iVx);
	T_nD_1D (i_Fg_id,iX,d_iVx);
	T_nD_1D (i_Fg_iu,iX,u_iVx);
	
	if ( DF_g_1D->at(i_Fg_id)!=0 && DF_g_1D->at(i_Fg_ii)!=0 && DF_g_1D->at(i_Fg_iu)!=0){
	  S_FlnF_temp[iX] = S_FlnF_temp[iX] + ( (
		      (		DF_g_1D->at(i_Fg_id) * log (DF_g_1D->at(i_Fg_id)) ) +
		      ( 4.0 * 	DF_g_1D->at(i_Fg_ii) * log (DF_g_1D->at(i_Fg_ii)) ) +
		      (		DF_g_1D->at(i_Fg_iu) * log (DF_g_1D->at(i_Fg_iu)) ) 
					)/3.0)*P->dVx;
	 }
      }
      int iVx0 = 0;
      T_nD_1D (i_Fg_i0,iX,iVx0);
      if (DF_g_1D->at(i_Fg_i0)!=0) S_FlnF_temp[iX] = S_FlnF_temp[iX] + (DF_g_1D->at(i_Fg_i0)* log(DF_g_1D->at(i_Fg_i0))) *P->dVx;	
    }

	//Integration over x
	S_FlnF_cpu = 0.0;
	
	
// 	double ghost_minus_2,ghost_minus_1,ghost_nX_cpu,ghost_nX_cpu_plus_1;
	P->Ghost_Grid_Points_Vector(S_FlnF_temp,P->nX_cpu,1);
	S_FlnF_cpu= S_FlnF_cpu + ( 
		      (P->minus_1+ 4 * S_FlnF_temp[0]+ S_FlnF_temp[1])
			  /3 )*P->dX;//CPU Boundary Correction
	
	S_FlnF_cpu= S_FlnF_cpu + ( 
		      (S_FlnF_temp[P->nX_cpu-1] + 4 *P->ghost_nX_cpu + P->nX_cpu_plus_1)
			  /3 )*P->dX;//CPU Boundary Correction
			  
			  
			  
	for(int iX=1; iX<=P->nX_cpu-2; iX++){  
	  S_FlnF_cpu= S_FlnF_cpu + ( 
			  (S_FlnF_temp[iX-1]+ 4 * S_FlnF_temp[iX]+ S_FlnF_temp[iX+1])
			  /3 )*P->dX;
	}
	 //-------------------------------------------------------------------
    
    }
  //***********************************************************************************************
  
  //***********************************************************************************************
  void Energy_kin_cpu(double & Ene_kin_cpu) {

    
    vector<double> Ene_kin_temp(P->nX_cpu+1);
    fill(Ene_kin_temp.begin(), Ene_kin_temp.end(), 0);

    //Integration over v
    int i_Fg_ii,d_iVx,u_iVx,i_Fg_i0,i_Fg_id,i_Fg_iu;
    for(int iX=0; iX<=P->nX_cpu-1; iX++){
      for(int iVx=1; iVx<=P->nVx-2; iVx++){
	d_iVx = iVx-1;
	u_iVx = iVx+1;
	T_nD_1D (i_Fg_ii,iX,iVx);
	T_nD_1D (i_Fg_id,iX,d_iVx);
	T_nD_1D (i_Fg_iu,iX,u_iVx);
	
	Ene_kin_temp[iX] = Ene_kin_temp[iX] + ((                                    
                           4*( pow((iVx*P->dVx+P->Vx_min),2) * (DF_g_1D->at(i_Fg_ii)))            
                           + pow(((iVx-1)*P->dVx+P->Vx_min),2) * DF_g_1D->at(i_Fg_id)           
                           + pow(((iVx+1)*P->dVx+P->Vx_min),2) * DF_g_1D->at(i_Fg_iu)            
                                             )/3*P->dVx);	

	 
      }
      int iVx0 = 0;
      T_nD_1D (i_Fg_i0,iX,iVx0);
      Ene_kin_temp[iX] = Ene_kin_temp[iX] +  pow (P->Vx_min,2) *DF_g_1D->at(i_Fg_i0)  * P->dVx;   
    }
    for(int iX=0; iX<=P->nX_cpu-1; iX++){
        Ene_kin_x_cpu->at(iX) = 0.5*P->Mass * Ene_kin_temp[iX];
    }
    
	//Integration over x
	Ene_kin_cpu = 0.0;
	
	
// 	double ghost_minus_2,ghost_minus_1,ghost_nX_cpu,ghost_nX_cpu_plus_1;
	P->Ghost_Grid_Points_Vector(Ene_kin_temp,P->nX_cpu,1);
	Ene_kin_cpu= Ene_kin_cpu + ( 
		      (P->minus_1+ 4 * Ene_kin_temp[0]+ Ene_kin_temp[1])
			  /3 )*P->dX;//CPU Boundary Correction
	
	Ene_kin_cpu= Ene_kin_cpu + ( 
		      (Ene_kin_temp[P->nX_cpu-1] + 4 *P->ghost_nX_cpu + P->nX_cpu_plus_1)
			  /3 )*P->dX;//CPU Boundary Correction
			  
			  
			  
	for(int iX=1; iX<=P->nX_cpu-2; iX++){  
	  Ene_kin_cpu= Ene_kin_cpu + ( 
			  (Ene_kin_temp[iX-1]+ 4 * Ene_kin_temp[iX]+ Ene_kin_temp[iX+1])
			  /3 )*P->dX;
	
        
    }
    
	  Ene_kin_cpu= 0.5*P->Mass * Ene_kin_cpu;        
    
	 //-------------------------------------------------------------------
    
    }
  //***********************************************************************************************
  
  //***********************************************************************************************
  void Write_Entropy_Temporal_Techplot() {
      file_Time_S_F2   << P->R_Time <<"   "<< S_F2_Tot<< endl;
      file_Time_S_FlnF << P->R_Time <<"   "<< S_FlnF_Tot<< endl;
      file_Time_Ene_kin << P->R_Time <<"   "<< Ene_kin_Tot<< endl;
  } 
  //***********************************************************************************************

  //***********************************************************************************************
  void PhasePoints_Sorting_Storing(hid_t & Group_ID_DF){
    
 	hsize_t Dim_DFp[1],Dim_DFp_Chunk[1];
 	hsize_t Count_DFp[1], OffSet_DFp[1], Stride_DFp[1], Block_DFp[1];
    
    
	vector<int>  PP_marker(PhasePoints->getFill());
	vector<int>  Num_PP_X(P->nX_cpu);
	vector<int>  marker_integer(P->nX_cpu); 
	int i_counter, lX;
	stringstream DSet_name;
//     	Storage<PhasePointClass> * PhasePoints_temp;
    
	vector<double>  X_temp_array(PhasePoints->getFill());
	vector<double> Vx_temp_array(PhasePoints->getFill());
	vector<double> Fp_temp_array(PhasePoints->getFill());


	int  sorted_i;
// 	int N_PP_cpu = PhasePoints->getFill();
	
	vector<int> Num_PhasePoints_array(size);
	vector<int> Integer_PhasePoints_array(size);

	int Num_PhasePoints, temp_storage, total_number_phasepoints;
      
        
	fill(Num_PhasePoints_array.begin(), Num_PhasePoints_array.end(), 0);
	fill(Integer_PhasePoints_array.begin(), Integer_PhasePoints_array.end(), 0);
	
	fill(PP_marker.begin(), PP_marker.end(), 0);
	fill(Num_PP_X.begin(), Num_PP_X.end(), 0);
	fill(marker_integer.begin(), marker_integer.end(), 0);
    
	fill(X_temp_array.begin(), X_temp_array.end(), 0.0);
	fill(Vx_temp_array.begin(), Vx_temp_array.end(), 0.0);
	fill(Fp_temp_array.begin(), Fp_temp_array.end(), 0.0);
	i_counter = 0;
	
	//--------------------------------------------------------------------------------------------------------------------------------
	for(Storage<PhasePointClass>::Iterator particle_iter = PhasePoints->begin(); particle_iter != PhasePoints->end(); ++particle_iter) { //interpolation main loop
            lX  = int (( (*particle_iter).Data_PP->at(P->x)  - P->X_min  ) / P->dX  ); //lx can be used to find out of boundry error for cpus
            PP_marker [i_counter]  = lX ;
	    i_counter += 1;
	    Num_PP_X[lX] +=1 ;
	    
	}
	//--------------------------------------------------------------------------------------------------------------------------------
	
	i_counter = 0;
	
	for (int i = 0; i<=P->nX_cpu-1; i++) {
	    for (int j = 0; j<=i; j++) {
		marker_integer [i] += Num_PP_X[j];	
	    }
	marker_integer [i] += -1 ;// if you have for 3rd grid point = 5+6+3=14, but you mean 0, 1, 2, ..., 12, 13. 
//       cout <<" rank="<<rank<<"  nX="<<i<<"   number of pp="<< Num_PP_X[i]<<"  integer of pp="<<marker_integer [i]<<endl;
//       if (marker_integer [i]>PhasePoints->getFill()){
// 	  cout << marker_integer [i] <<" rank="<<rank<<"  nX="<< PP_marker [i_counter]<<"         ";
//       }
	}


    
	for(Storage<PhasePointClass>::Iterator particle_iter = PhasePoints->begin(); particle_iter != PhasePoints->end(); ++particle_iter) { //interpolation main loop        
	    sorted_i = marker_integer[PP_marker [i_counter]];
	    //if (sorted_i>PhasePoints->getFill()){
	    //	 cout << marker_integer [PP_marker [i_counter]]<<"=expcted, ... recieved"<<sorted_i <<" rank="<<rank<<"  nX="<< PP_marker [i_counter]<<" bigger than"<<PhasePoints->getFill()<<"        ";
	    //}
	
	    X_temp_array [sorted_i]  = (*particle_iter).Data_PP->at(P->x);
	    Vx_temp_array[sorted_i]  = (*particle_iter).Data_PP->at(P->vx);
	    Fp_temp_array[sorted_i]  = (*particle_iter).Data_PP->at(P->fp);
	    marker_integer[PP_marker [i_counter]] +=-1;
	    i_counter += 1; // should be after the marker_integer[PP_marker [i_counter]] +=1;
	}
	
	
      // to replace the phasepoints with their sorted one ---------------------------------------------------------------------------------------------	
//     PhasePoints_temp = new Storage<PhasePointClass>;
    
// 	for(Storage<PhasePointClass>::Iterator particle_iter = PhasePoints->begin(); particle_iter != PhasePoints->end(); ++particle_iter) { //interpolation main loop        
// 	  particle_iter.pop();
// 	}
//     cout <<PhasePoints->getFill()<<"  rank="<<rank<<endl;
//     for (int i_pp = 0; i_pp<N_PP_cpu; i_pp++) {
//       PhasePointClass PhasePoint (X_temp_array[i_pp],Vx_temp_array[i_pp],Fp_temp_array[i_pp]); 
//       PhasePoints->push_back(PhasePoint);
//       //PhasePoints_temp->push_back(PhasePoint);
//     }
//     if (PhasePoints->getFill() != PhasePoints_temp->getFill()) {cout<<"errrrrror in rank"<<rank<<"   "<<PhasePoints->getFill()<<"="<< PhasePoints_temp->getFill();}
//     PhasePoints = PhasePoints_temp;
//     cout <<"finished the job at rank"<<rank<<endl;
//     delete PhasePoints_temp;
      // -------------------------------------------- ---------------------------------------------------------------------------------------------	

    
    //===========================================================================================      
    Num_PhasePoints = PhasePoints->getFill();
    
    //-----------------------------------------------------------------------------------------------
    MPI_Gather(&Num_PhasePoints, 1, MPI_INT,&Num_PhasePoints_array[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank==0) {
	    total_number_phasepoints=0;
	    for (int i=0; (i<=NODE); ++i){
		total_number_phasepoints +=  Num_PhasePoints_array[i];
	    } 
	    temp_storage = 0;
	    Integer_PhasePoints_array [0] = 0;
	    for (int i=1; (i<=NODE); ++i){
		temp_storage +=  Num_PhasePoints_array[i-1];
		Integer_PhasePoints_array [i] = temp_storage;
	    } 
     }

	MPI_Bcast( &Integer_PhasePoints_array[0], size,  MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast( &total_number_phasepoints, 1,  MPI_INT, 0, MPI_COMM_WORLD);

/*
	vector<double> * test;
	
	
 	if (rank == 0){test  = new vector<double>(10);}else if (rank == 1){test  = new vector<double>(10);}else if (rank == 2){test  = new vector<double>(10);}else{test  = new vector<double>(11);}
 	fill(test->begin(), test->end(), (rank+1)*111.0);
   	Dim_DFp [0]=41; if (rank == 0){Dim_DFp_Chunk [0] = 10;}else if (rank == 1){Dim_DFp_Chunk [0] = 10;}else if (rank == 2){Dim_DFp_Chunk [0] = 10;}else{ Dim_DFp_Chunk [0] = 11;}
   	if (rank == 0){OffSet_DFp [0] = 0;}else if (rank == 1){OffSet_DFp [0] = 10;}else if (rank == 2){OffSet_DFp [0] = 20;}else{OffSet_DFp [0] = 30;}
  	
  	
    
 	if (rank == 0){test  = new vector<double>(8);}else if (rank == 1){test  = new vector<double>(7);}else if (rank == 2){test  = new vector<double>(13);}else{test  = new vector<double>(4);}
	fill(test->begin(), test->end(), (rank+1)*111.0);
  	Dim_DFp [0]=32; if (rank == 0){Dim_DFp_Chunk [0] = 8;}else if (rank == 1){Dim_DFp_Chunk [0] = 7;}else if (rank == 2){Dim_DFp_Chunk [0] = 13;}else{ Dim_DFp_Chunk [0] = 4;}
  	if (rank == 0){OffSet_DFp [0] = 0;}else if (rank == 1){OffSet_DFp [0] = 8;}else if (rank == 2){OffSet_DFp [0] = 15;}else{OffSet_DFp [0] = 28;}

  */  	
	Dim_DFp [0]=total_number_phasepoints;
  	Dim_DFp_Chunk [0] = Num_PhasePoints;
	OffSet_DFp [0] = Integer_PhasePoints_array[rank];
	
	
	Block_DFp [0] = Dim_DFp_Chunk [0];
	Count_DFp[0]= 1;
	Stride_DFp[0]=1;
// 	cout << fixed <<" rank = "<<rank<<"   from "<<OffSet_DFp [0]<<"   end at "<<OffSet_DFp [0]+Dim_DFp_Chunk [0]-1<<endl;
  
	
      //============================================================================================
    // 2) Create the dataspace for the dataset  -------------------------------
	P->Space_ID_DFp_Num = H5Screate_simple(P->Rank_X, P->Dim_X , NULL); 
	P->Memory_ID_DFp_Num = H5Screate_simple(P->Rank_X, P->Dim_X_Chunk, NULL);
      //-------------------------------------------------------------------------
   
      // 3) Create chunked dataset ----------------------------------------------
	P->PList_ID = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(P->PList_ID, P->Rank_X, P->Dim_X_Chunk);
    
	
	DSet_name << "Num_"<<P->Name; 

	P->DSet_ID_DFp_Num = H5Dcreate(Group_ID_DF, DSet_name.str().c_str(), H5T_NATIVE_DOUBLE, P->Space_ID_DFp_Num,
			  H5P_DEFAULT, P->PList_ID, H5P_DEFAULT);
	H5Pclose(P->PList_ID);
    //------------------------------------------------------------------------
    H5Sclose(P->Space_ID_DFp_Num);
  
    // 4) Select hyperslab in the file----------------------------------------
    P->Space_ID_DFp_Num = H5Dget_space(P->DSet_ID_DFp_Num);
    P->status_hdf5 = H5Sselect_hyperslab(P->Space_ID_DFp_Num, H5S_SELECT_SET, P->OffSet_X, P->Stride_X, P->Count_X, P->Block_X);
    //-----------------------------------------------------------------------

    // 5) Create property list for collective dataset write------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(P->PList_ID, H5FD_MPIO_COLLECTIVE);
    //-----------------------------------------------------------------------

  
    
    
    P->status_hdf5 = H5Dwrite(P->DSet_ID_DFp_Num, H5T_NATIVE_INT, P->Memory_ID_DFp_Num, P->Space_ID_DFp_Num,
		      P->PList_ID, Num_PP_X.data());
    
    

    
    
   
    //H5Fflush(File_ID_S, H5F_SCOPE_GLOBAL ); 
   
    // Close/release resources    
    H5Dclose(P->DSet_ID_DFp_Num);
    
    H5Pclose(P->PList_ID);
    
    H5Sclose(P->Space_ID_DFp_Num);
    H5Sclose(P->Memory_ID_DFp_Num);
  
   
   
    
    //============================================================================================
      //-------------------------------------------------------------------------
      //#ifdef debug
	P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="********  Write DFp HDF5 ");
      //#endif
	
    // 2) Create the dataspace for the dataset  -------------------------------
	P->Space_ID_DFp_x = H5Screate_simple(P->Rank_DFp, Dim_DFp , NULL); 
	P->Memory_ID_DFp_x = H5Screate_simple(P->Rank_DFp, Dim_DFp_Chunk, NULL);
      //-------------------------------------------------------------------------
   
      // 3) Create chunked dataset ----------------------------------------------
	P->PList_ID = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(P->PList_ID, P->Rank_DFp, Dim_DFp_Chunk);
	DSet_name.str("");
	DSet_name << "X_"<<P->Name; 

	P->DSet_ID_DFp_x = H5Dcreate(Group_ID_DF, DSet_name.str().c_str(), H5T_NATIVE_DOUBLE, P->Space_ID_DFp_x,
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); //one before the last one should be H5P_DEFAULT, so irregular chunk are written correctly, 19 Sep 2016, 3 days
	H5Pclose(P->PList_ID);
    //------------------------------------------------------------------------
    H5Sclose(P->Space_ID_DFp_x);
  
    // 4) Select hyperslab in the file----------------------------------------
    P->Space_ID_DFp_x = H5Dget_space(P->DSet_ID_DFp_x);
    P->status_hdf5 = H5Sselect_hyperslab(P->Space_ID_DFp_x, H5S_SELECT_SET, OffSet_DFp, Stride_DFp, Count_DFp, Block_DFp);
    //-----------------------------------------------------------------------

    // 5) Create property list for collective dataset write------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(P->PList_ID, H5FD_MPIO_COLLECTIVE);
    //-----------------------------------------------------------------------

  
    
    P->status_hdf5 = H5Dwrite(P->DSet_ID_DFp_x, H5T_NATIVE_DOUBLE, P->Memory_ID_DFp_x, P->Space_ID_DFp_x,
		      P->PList_ID, X_temp_array.data());

   
    //H5Fflush(File_ID_S, H5F_SCOPE_GLOBAL ); 
    // Close/release resources    
    H5Dclose(P->DSet_ID_DFp_x);
 
    H5Pclose(P->PList_ID);
    
    H5Sclose(P->Space_ID_DFp_x);
    H5Sclose(P->Memory_ID_DFp_x);
    
    
    
    //============================================================================================
      //-------------------------------------------------------------------------
      //#ifdef debug
	P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="********  Write DFp HDF5 ");
      //#endif
	
    // 2) Create the dataspace for the dataset  -------------------------------
	P->Space_ID_DFp_x = H5Screate_simple(P->Rank_DFp, Dim_DFp , NULL); 
	P->Memory_ID_DFp_x = H5Screate_simple(P->Rank_DFp, Dim_DFp_Chunk, NULL);
      //-------------------------------------------------------------------------
   
      // 3) Create chunked dataset ----------------------------------------------
	P->PList_ID = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(P->PList_ID, P->Rank_DFp, Dim_DFp_Chunk);
	DSet_name.str("");
	DSet_name << "Vx_"<<P->Name; 

	P->DSet_ID_DFp_x = H5Dcreate(Group_ID_DF, DSet_name.str().c_str(), H5T_NATIVE_DOUBLE, P->Space_ID_DFp_x,
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); //one before the last one should be H5P_DEFAULT, so irregular chunk are written correctly, 19 Sep 2016, 3 days
	H5Pclose(P->PList_ID);
    //------------------------------------------------------------------------
    H5Sclose(P->Space_ID_DFp_x);
  
    // 4) Select hyperslab in the file----------------------------------------
    P->Space_ID_DFp_x = H5Dget_space(P->DSet_ID_DFp_x);
    P->status_hdf5 = H5Sselect_hyperslab(P->Space_ID_DFp_x, H5S_SELECT_SET, OffSet_DFp, Stride_DFp, Count_DFp, Block_DFp);
    //-----------------------------------------------------------------------

    // 5) Create property list for collective dataset write------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(P->PList_ID, H5FD_MPIO_COLLECTIVE);
    //-----------------------------------------------------------------------

  
    
    P->status_hdf5 = H5Dwrite(P->DSet_ID_DFp_x, H5T_NATIVE_DOUBLE, P->Memory_ID_DFp_x, P->Space_ID_DFp_x,
		      P->PList_ID, Vx_temp_array.data());

   
    //H5Fflush(File_ID_S, H5F_SCOPE_GLOBAL ); 
    // Close/release resources    
    H5Dclose(P->DSet_ID_DFp_x);
 
    H5Pclose(P->PList_ID);
    
    H5Sclose(P->Space_ID_DFp_x);
    H5Sclose(P->Memory_ID_DFp_x);
    
    //============================================================================================
      //-------------------------------------------------------------------------
      //#ifdef debug
	P->Print_Debug(P->line_NO=__LINE__,P->Debug_mes="********  Write DFp HDF5 ");
      //#endif
	
    // 2) Create the dataspace for the dataset  -------------------------------
	P->Space_ID_DFp_x = H5Screate_simple(P->Rank_DFp, Dim_DFp , NULL); 
	P->Memory_ID_DFp_x = H5Screate_simple(P->Rank_DFp, Dim_DFp_Chunk, NULL);
      //-------------------------------------------------------------------------
   
      // 3) Create chunked dataset ----------------------------------------------
	P->PList_ID = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(P->PList_ID, P->Rank_DFp, Dim_DFp_Chunk);
	DSet_name.str("");
	DSet_name << "Fp_"<<P->Name; 

	P->DSet_ID_DFp_x = H5Dcreate(Group_ID_DF, DSet_name.str().c_str(), H5T_NATIVE_DOUBLE, P->Space_ID_DFp_x,
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); //one before the last one should be H5P_DEFAULT, so irregular chunk are written correctly, 19 Sep 2016, 3 days
	H5Pclose(P->PList_ID);
    //------------------------------------------------------------------------
    H5Sclose(P->Space_ID_DFp_x);
  
    // 4) Select hyperslab in the file----------------------------------------
    P->Space_ID_DFp_x = H5Dget_space(P->DSet_ID_DFp_x);
    P->status_hdf5 = H5Sselect_hyperslab(P->Space_ID_DFp_x, H5S_SELECT_SET, OffSet_DFp, Stride_DFp, Count_DFp, Block_DFp);
    //-----------------------------------------------------------------------

    // 5) Create property list for collective dataset write------------------
    P->PList_ID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(P->PList_ID, H5FD_MPIO_COLLECTIVE);
    //-----------------------------------------------------------------------

  
    
    P->status_hdf5 = H5Dwrite(P->DSet_ID_DFp_x, H5T_NATIVE_DOUBLE, P->Memory_ID_DFp_x, P->Space_ID_DFp_x,
		      P->PList_ID, Fp_temp_array.data());

   
    //H5Fflush(File_ID_S, H5F_SCOPE_GLOBAL ); 
    // Close/release resources    
    H5Dclose(P->DSet_ID_DFp_x);
 
    H5Pclose(P->PList_ID);
    
    H5Sclose(P->Space_ID_DFp_x);
    H5Sclose(P->Memory_ID_DFp_x);
  }
  //***********************************************************************************************
   
 
  };
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#endif //SpeciesClass_Included
