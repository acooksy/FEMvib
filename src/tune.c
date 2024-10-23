	#include <iostream>
	#include <fstream>
	#include "nr3.h"
	#include "ludcmp.h"
	#include "krig.h"
	#include "pointSelector2.h"
	#include "parameter.h"
	#include <cstdlib>
	
using namespace std;

int main(int argc, char** argv)
{

Int dim = atoi(argv[1]);
Doub rparam = atof(argv[2]);
Int numP = atoi(argv[3]);

	// call the Kriging function:
	 krigcall(dim,rparam,numP); 

  return 0;

}

