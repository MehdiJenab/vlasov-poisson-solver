
// Saturday 4 August 2018, 2018-08-04, 12:42
// Umeå University
// Umeå, Sweden


//  #define debug_initial
//  #define debug
//  #define debug_Poisson


#define long_time_run


// to produce the reflected population distribution function 
// = 0 simple reflection
// = 1 fully mixed version
// #define flag_reflected_distribution_full_mixture_mehtod
#define flag_reflected_distribution_simple_cut

// one the following needs to be active when inserting a hole of trapped population inside the simulation box. 
#define single_beta_for_trapped_DF
// #define multiple_beta_for_trapped_DF // this needs a file called Jenab_Distribution_Function.json

#define flag_Poisson_Solver 0
//-------------------------------------

// #define flag_boundary_condition_open
#define flag_boundary_condition_periodic
// boundary_condition = 0 --> Periodic
// boundary_condition = 1 --> Open

// #define Co_moving_frame

//-------------------------------------
#define print_output
#ifdef print_output
    #define HDF5_Usage
    #ifdef HDF5_Usage
	#define print_Spatial_HDF5
	#define print_Contour_HDF5
	#ifdef print_Contour_HDF5 
			#define print_FPP
			#define print_EPP
// 			#define print_MPP
	#endif
	#define print_energy_kinetic_in_X_dimesion
	//#define print_Temporal_HDF5 //not working yet (8-Feb-2016)
    #endif
    #define Techplot_Usage
    #ifdef Techplot_Usage
	#define print_Temporal_Techplot
    #endif
#endif



// Timing flags =0 total time, number of phase points, =1 add detail of timing of each function
#define flag_timing_of_run 0



  

      
     

