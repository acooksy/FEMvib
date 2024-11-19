//
//  gmatUtility.h
//  GMAT
//
//  Created by PETER ZAJAC on 10/1/2011.
//  Copyright 2011 SDSU. All rights reserved.
//

void gmatF()
{
    int i,j,k,l,dim,natom,numStructures;
    int g=1;
    ifstream indata;
   
    
    // form of input: 
    // 1.line: dimension of internal coordinates, number of atoms, number of structures
    // 1. Geometry is the reference with the masses of the atoms at the beginning
    // then it is just the internal coordinates followed by the cartesians
    
    // Examples in 2D internals "CO2" : 
    /*
    2 3 36
    16. 0.0000000000 0.0000000000 -1.1603899217
    12. 0.0000000000 0.0000000000 0.0000000445
    16. 0.0000000000 0.0000000000 1.1603898883
    1.074 0.501
    0.0000000000 0.0000000000 -0.8656905735
    0.0000000000 0.0000000000 0.2083094265
    0.0000000000 0.0000000000 0.7093094265
    1.084 0.684
    0.0000000000 0.0000000000 -0.9385832974
    0.0000000000 0.0000000000 0.1454167026
    0.0000000000 0.0000000000 0.8294167026
    */
    
    
    
    // read in the reference and all the internal coordinates along with 
     indata.open("/home/acooksy/femvib/results/coordinates");
    if (!indata) 
            {
                std::cerr << "Error: file could not be opened\n"; exit(1);
            }
    
        indata >> dim >> natom >> numStructures;
        MatDoub xreference(1,3*natom); // this will hold the reference geometry
        VecDoub masses(natom); // vector of masses of individual atoms
        MatDoub xgeom(numStructures,natom*3); // matrix of individual geometries
        MatDoub internalcoords(numStructures,dim); // matrix of internal coordinates
    
        for (i=0;i<natom;i++)
            {
                l=3*i;
                cout.precision(10);
                cout.setf(ios::fixed,ios::floatfield);
                indata >> masses[i] >> xreference[0][0+l] >> xreference[0][1+l] >> xreference[0][2+l];
            } // this loop read in the reference geometries and the masses
        l=0;
        if (dim == 3) // check if dimensionality of internal coordinates == 3
            {
                for (j=0; j<numStructures; j++) // this loop cycles through each Structure
                {
                    indata >> internalcoords[j][0] >> internalcoords[j][1] >> internalcoords[j][2];  
                    for (k=0; k<natom; k++) // within the Structure this loop cycles through individual coorinates
                        {
                            l=3*k;
                            indata >> xgeom[j][0 + l] >> xgeom[j][1 + l] >> xgeom[j][2 + l];
                        }
                }
            }
    
        if (dim == 2) // or if dimensionality of internal coordinates == 2 
                {
                    for (j=0; j<numStructures; j++) // this loop cycles through each geometry
                    {
                        indata >> internalcoords[j][0] >> internalcoords[j][1];  
                        for (k=0; k<natom; k++) // within the geometry this loop cycles through individual coorinates
                        {
                            l=3*k;
                            indata >> xgeom[j][0 + l] >> xgeom[j][1 + l] >> xgeom[j][2 + l];
                        }
                    }
                }
    
    indata.close(); // close the file   
    
    // pre-build Jacobian matrix:
    MatDoub jacob(numStructures,dim*3*natom);
    MatDoub Gmat(numStructures,dim*dim);


       if (natom == 3 && dim == 3) 
    {
        gmatAnalytical(Gmat, internalcoords, masses, numStructures);
        printGmat(Gmat, internalcoords, dim, numStructures);
        exit(0);
    } 
    
    Int st=1;
    Swapper(xreference, natom, st); 
    Swapper(xgeom, natom, numStructures); 
    
    prerot(xgeom, xreference, natom, numStructures);
    
    
    // At this stage of the program we have read in all the necessary data to perform the gmatrix transformation
    // Step no 1 is to get rid of translational and rotational motion by impossing Eckart Conditions:
    
    // correct for translation:
    translation(g,natom,xreference,masses);
    translation(numStructures,natom,xgeom,masses);

    
    //correct for rotation there are two options still in testing phase:
    rotation(xreference, numStructures, natom, xgeom, masses);
    
    // checking if the rotational condition is satisfied (commented for now):
    // rottest(xgeom, xreference, masses, numStructures, natom);
    
    
    // ==== By observing the data I noticed that some of the structures have an exchanged sign ====== //
    // ==== this makes the gradients at many points extremly large and inaccurate =================== //
    // ===== the routine below takes care of the signs and plays around with the data =============== //
     playdata(xgeom,natom,numStructures,xreference);

    // Assemble the Jakobian Matrix
    jakobAssembly(jacob, internalcoords, xgeom, numStructures, natom, dim);

    
    // Calculate G-Matrix
    calcGmat(jacob, masses, dim, numStructures, natom, Gmat);
 
    
    // print the Gmatrix:
      printGmat(Gmat, internalcoords, dim, numStructures);
    

    // Program finished
}
