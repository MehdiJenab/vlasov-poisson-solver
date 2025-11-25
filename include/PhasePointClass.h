#ifndef PhasePointClass_Included
#define PhasePointClass_Included

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
class  PhasePointClass{  

public:
	
	vector<double> * Data_PP;
	
	// Data_PP = [Fp, X, X0, Vx, Vx0, EPP, MPP]
	int fp = 0, x = 1, x0 = 2, vx=3, vx0=4;
	#ifdef print_EPP
	int epp=5; 
	#endif
	#ifdef print_MPP
	int mpp=6; 
	#endif


	//=====================================================================
	PhasePointClass(vector<double> & data_in, int& number_of_data  ){
		Data_PP =  new vector<double>;
		Data_PP->resize(number_of_data);
		for (int i=0; i<number_of_data; ++i){
			Data_PP->at(i) = data_in[i];
		}
	}
	

// 	~PhasePointClass(){ // don't use deconstructor since the Data_pp are stored in the storage class as PhasePoints
// 		delete Data_PP;
// 	}

	
	void delete_Data_PP(){
		delete Data_PP;
	}
	
	//*********************************************************************
	void Xpush_new (int& P_x_out, int& P_x_in, int& P_v_in, const double& d_Time, ParametersClass * P) {
		Data_PP->at(P_x_out) = Data_PP->at(P_x_in) + Data_PP->at(P_v_in) * (d_Time); //	X = X + Vx * (d_Time); //Xpush		
		Boundary_Condition_function_new(P_x_out,P);
	}
	//*********************************************************************
	//*********************************************************************
	void Xpush (ParametersClass * P, const double& d_Time) {
		Data_PP->at(x) += Data_PP->at(vx) * (d_Time); //	X += Vx * (d_Time); //Xpush		
		Boundary_Condition_function(P);
	}
	//*********************************************************************
  
	//*********************************************************************
	void Xpush_trapezoidal (ParametersClass * P, const double& d_Time) {
		double X1 = Data_PP->at(x);
		Data_PP->at(x) += Data_PP->at(vx) * (d_Time);	// X += Vx * (d_Time); //Xpush
		Data_PP->at(x) = (X1+Data_PP->at(x))/2.0; 	// X = (X1+X)/2.0;
		Boundary_Condition_function(P);
	}
	//*********************************************************************
	
	
	//*********************************************************************
	void Xpush_over_X0_by_Vxhalf (ParametersClass * P, const double& d_Time) {
		Data_PP->at(x) = Data_PP->at(x0) + Data_PP->at(vx0) * (d_Time);
		Boundary_Condition_function(P);
	}
	//*********************************************************************
	
	//*********************************************************************
	void average_V_V0 (int& P_vx_out, int& P_v_in_1, int& P_v_in_2,  ParametersClass * P) {
		Data_PP->at(P_vx_out) = ( Data_PP->at(P_v_in_1) + Data_PP->at(P_v_in_2) )/2.0;
	}
	//*********************************************************************
	
	
	
	//*********************************************************************
	void Boundary_Condition_function (ParametersClass * P) {
		#ifdef flag_boundary_condition_open
		boundary_condition_open(P);
		#endif
		#ifdef flag_boundary_condition_periodic
		boundary_condition_periodic(P);
		#endif
		
	}
	
	//*********************************************************************
	void Boundary_Condition_function_new (int & P_x_out, ParametersClass * P) {
		#ifdef flag_boundary_condition_open
		boundary_condition_open(P);
		#endif
		#ifdef flag_boundary_condition_periodic
		boundary_condition_periodic_new(P_x_out,P);
		#endif
		
	}
	
	//*********************************************************************
	void Calculate_Vxhalf_Store_in_Vx0(){
		Data_PP->at(vx0) = (Data_PP->at(vx0) + Data_PP->at(vx))/2.0;
	}
	
	
	//*********************************************************************
	
	//*********************************************************************
	void Boundary_Condition_Vx_new (int& P_vx_out, ParametersClass * P, Storage<PhasePointClass>::Iterator & particle_iter) {
		if (Data_PP->at(P_vx_out)<P->Vx_min || Data_PP->at(P_vx_out)>P->Vx_max) {
			//PhasePointClass PhasePoint_delete = particle_iter.pop();
			delete Data_PP;
			particle_iter.pop(); // CAUSION I am not sure if this is synced, probably it is better to have this out
						// of here in the future Phase_Point_Function class
		}
		
	}
	
	//*********************************************************************
	void Boundary_Condition_Vx (ParametersClass * P, Storage<PhasePointClass>::Iterator & particle_iter) {
		if (Data_PP->at(vx)<P->Vx_min || Data_PP->at(vx)>P->Vx_max) {
			//PhasePointClass PhasePoint_delete = particle_iter.pop();
			delete Data_PP;
			particle_iter.pop(); // CAUSION I am not sure if this is synced, probably it is better to have this out
						// of here in the future Phase_Point_Function class
		}
		
	}
	//*********************************************************************
	void Vpush_trapezoidal (ParametersClass * P, Storage<PhasePointClass>::Iterator & particle_iter,const double& d_Time, const double& Ex ) {
		double Vx1 = Data_PP->at(vx);
		Data_PP->at(vx) += (P->Charge/P->Mass)*(Ex*d_Time);
		Data_PP->at(vx) = (Data_PP->at(vx)+Vx1)/2.0;
// 		Vx = (Vx+Vx1)/2.0;
		Boundary_Condition_Vx(P,particle_iter);

	}
	//*********************************************************************

	
	//*********************************************************************
	void Vpush_over_Vx0 (ParametersClass * P, Storage<PhasePointClass>::Iterator & particle_iter,const double& d_Time, const double& Ex ) {
		Data_PP->at(vx) = Data_PP->at(vx0) + (P->Charge/P->Mass)*(Ex*d_Time); 
		Boundary_Condition_Vx(P,particle_iter);
	}
	//*********************************************************************
	
	
	//*********************************************************************
	void Vpush (ParametersClass * P, Storage<PhasePointClass>::Iterator & particle_iter,const double& d_Time, const double& Ex ) {
		Data_PP->at(vx) += (P->Charge/P->Mass)*(Ex*d_Time); 
		Boundary_Condition_Vx(P,particle_iter);

	}
	//*********************************************************************
  
	//*********************************************************************
	void Vpush_new (int& P_vx_out, int P_vx_in, ParametersClass * P, Storage<PhasePointClass>::Iterator & particle_iter,const double& d_Time, const double& Ex ) {
		Data_PP->at(P_vx_out) = Data_PP->at(P_vx_in) + (P->Charge/P->Mass)*(Ex*d_Time); 
		Boundary_Condition_Vx_new(P_vx_out,P,particle_iter);

	}
	//*********************************************************************
	
	
	//*********************************************************************
	void Save_X0_Vx0 () {
		Data_PP->at(x0) = Data_PP->at(x);
		Data_PP->at(vx0) = Data_PP->at(vx);
	}
	//*********************************************************************
	void Maxwellian_distribution_function(ParametersClass * P){
// 		double Vx_on_grid;
// 		Vx_on_grid = ((float(int(Data_PP->at(vx) - P->Vx_min)/ P->dVx))*P->dVx) + P->Vx_min; // this is very expensive
		Data_PP->at(fp) = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) * 
			exp( - 0.5 * pow(Data_PP->at(vx),2) * P->norm_factor); //energy_kinetic/Temperature
	}
	
	//*********************************************************************
	void boundary_condition_open(ParametersClass * P){		
		if (Data_PP->at(x) >= P->X_max_Tot) {
			while (Data_PP->at(x)>=P->X_max_Tot) { 
				Data_PP->at(x)=P->X_min_Tot + (Data_PP->at(x)-P->X_max_Tot);
			}
			
			Maxwellian_distribution_function(P);
			#ifdef print_EPP
			Data_PP->at(epp)=((P->Mass/2.0) * pow(Data_PP->at(vx),2)) / P->Temp;//4000.0;//+1000.0*P->R_Time;
			#endif
			
			#ifdef print_MPP
			Data_PP->at(mpp)=0.0;//4000.0;
			#endif
		}
		if (Data_PP->at(x) < P->X_min_Tot) {
			while (Data_PP->at(x)<P->X_min_Tot) { 
				Data_PP->at(x)=P->X_max_Tot + (Data_PP->at(x)-P->X_min_Tot);
			}
			Maxwellian_distribution_function(P);
			
			#ifdef print_EPP
			Data_PP->at(epp)=((P->Mass/2.0) * pow(Data_PP->at(vx),2)) / P->Temp;//5000.0+1000.0*P->R_Time;
			#endif
			
			#ifdef print_MPP
			Data_PP->at(mpp)=0.0;//5000.0;
			#endif
		}
	}
	//*********************************************************************

	//*********************************************************************
	void boundary_condition_periodic_new(int& P_x_out, ParametersClass * P){
		while (Data_PP->at(P_x_out)>=P->X_max_Tot) { //periodic boundary condition
				Data_PP->at(P_x_out)=P->X_min_Tot + (Data_PP->at(P_x_out)-P->X_max_Tot);
		}
		while (Data_PP->at(P_x_out)<P->X_min_Tot) { 
				Data_PP->at(P_x_out)=P->X_max_Tot + (Data_PP->at(P_x_out)-P->X_min_Tot);
		}
	}
	//*********************************************************************	
	
	//*********************************************************************
	void boundary_condition_periodic(ParametersClass * P){
		while (Data_PP->at(x)>=P->X_max_Tot) { //periodic boundary condition
				Data_PP->at(x)=P->X_min_Tot + (Data_PP->at(x)-P->X_max_Tot);
		}
		while (Data_PP->at(x)<P->X_min_Tot) { 
				Data_PP->at(x)=P->X_max_Tot + (Data_PP->at(x)-P->X_min_Tot);
		}
// 		while (X>=P->X_max_Tot) { X=P->X_min_Tot + (X-P->X_max_Tot);}//periodic boundary condition
// 		while (X<P->X_min_Tot) { X=P->X_max_Tot + (X-P->X_min_Tot);}
	}
	//*********************************************************************	
	//*********************************************************************
// 	void boundary_condition_reflective(){ // this is not developed ... 
	
			// reflective boundary for phase points
		// 	if (X>P->X_max_Tot) {
		// 		X=P->X_max_Tot - (X- (int(X/P->X_lenght)*P->X_lenght));
		// 		Vx =  (-(Vx + 12.0))-12 ;
		// 		Fp = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) * 
		// 		exp( - ((P->Mass/2.0) * pow(Vx+12,2)) / P->Temp); //energy_kinetic/Temperature
		// 		#ifdef print_EPP
		// 			EPP=4000.0;
		// 		#endif
		// 		#ifdef print_MPP
		// 			MPP=4000.0;
		// 		#endif
		// 	}
		// 	if (X<P->X_min_Tot) {
		// 		X=P->X_min_Tot + ((int(X/P->X_lenght)*P->X_lenght) - X);
		//  		Vx = (-(Vx + 12.0))-12 ;
		// 		Fp = (P->Den_Ratio_b) * sqrt(1.0/(2.0*M_PI)) * sqrt(P->norm_factor) * 
		// 		exp(- ( (P->Mass/2.0) * pow(Vx+12,2) ) /P->Temp); //energy_kinetic/Temperature
		// 		#ifdef print_EPP
		// 			EPP=5000.0;
		// 		#endif
		// 		#ifdef print_MPP
		// 			MPP=5000.0;
		// 		#endif
		// 	}
// 	}
		
  };
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
#endif //PhasePointClass_Included
