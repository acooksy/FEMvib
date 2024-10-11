//
//  main.cpp
//  GMAT
//
//  Created by PETER ZAJAC on 10/4/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "nr3.h"
#include "ludcmp.h"
#include "qrdcmp.h"
#include "krig.h"
#include "pointSelector2.h"
#include "geomswap.h"
#include "svd.h"
#include "eigen_sym.h"
#include "testRot.h"
#include "mins.h"
#include "mins_ndim.h"
#include "roots_multidim.h"
#include "quasinewton.h"
#include "interp_rbf.h"
#include "helpf.h"
#include "checkinversion.h"
#include "newFunc.h"
#include "eckart.h"
#include "imls2.h"
#include "jakobAssembly.h"
#include "gmatUtility.h"


using namespace std;

int main (int argc, const char * argv[])
{
    
    // This part will not contain a code since it needs to be automated it will only assemble the G-Matrix
    // by calling the function
    gmatF();
    
    
    
    return 0;
}

