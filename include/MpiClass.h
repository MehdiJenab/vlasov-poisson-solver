#ifndef MpiClass_Included
#define MpiClass_Included
#include "naming.h"

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
class MpiClass {
		//4) parallel parameters
public:
	int size, rank,NODE;
	MPI_Comm  COMM_CART;
	int * Left_cpu;
	int * Right_cpu;
	

	// --------------------------------------------------------------------------	
	void Initialize_MPI(){	
		// Initialize the MPI environment
		MPI_Init(NULL, NULL);

		// Get the number of processes
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		// Get the rank of the process
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		// Get the name of the processor
		char processor_name[MPI_MAX_PROCESSOR_NAME];
		int name_len;
		MPI_Get_processor_name(processor_name, &name_len);
		NODE=size-1;
		
		//---- Cartesian Coordiantes ----------------------------------------------------------
		Left_cpu = new int [size];
		Right_cpu = new int [size];
		int coords;
		
		int Cart_Dim [1]= {size};
		int PERIODS[1]= {1}; //1=true
		
		
		MPI_Cart_create(MPI_COMM_WORLD,1,Cart_Dim,PERIODS,false,&COMM_CART);
		MPI_Cart_coords(COMM_CART, rank, 1,&coords);
		#ifdef debug
		cout <<"my rank is="<<rank<<"  and my coord is="<<coords<<endl;
		#endif
		MPI_Cart_shift(COMM_CART,0,1, &Left_cpu[rank], &Right_cpu[rank] );

		
		#ifdef debug
		printf("left = %i,  rank: %i, righ = %i, \n", Left_cpu[rank],rank, Right_cpu[rank]);
		#endif
	
	}
	//------------------------------------------------------------------------------------ 
	
	//------------------------------------------------------------------------------------ 
	void Finalize_MPI(){
		 MPI_Finalize();
	} 
	//------------------------------------------------------------------------------------ 
};
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#endif // MpiClass_Included
