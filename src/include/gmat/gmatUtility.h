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
    Int molecule_type=0,coordinates_type=0;
    long int zeros=0;
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
    
//   indata.open("/home/acooksy/FEMvib/results/coordinates");  
   indata.open("./results/coordinates");  
    
    // read in the reference and all the internal coordinates along with 
    // indata.open("coordinates");
    // indata.open("/Users/peter/Projects/GMAT/GMAT/filein.txt");
    
    //  indata.open("/Users/peter/Projects/GMAT/GMAT/BLBA"); // artificial H2O, two bond lengths - excellent results
   //   indata.open("/Users/peter/Projects/GMAT/GMAT/BLBA2"); // artificial H2O, BL BA these resuts seem also very good.
     // indata.open("/Users/peter/Projects/GMAT/GMAT/HNOarti"); // WORKS
     // indata.open("/Users/peter/Projects/GMAT/GMAT/rotating_rod"); // hno looks good too
   //   indata.open("/Users/peter/Projects/GMAT/GMAT/coordinates_sin"); // hno looks good too
    //  indata.open("/Users/peter/Projects/GMAT/GMAT/hno"); // hno looks good too
  //    indata.open("/Users/peter/Projects/GMAT/GMAT/2D_introt_model_constG.xyz");
  //   indata.open("/Users/peter/Projects/GMAT/GMAT/2D_introt_model_varyG.xyz"); 
    //    indata.open("/Users/peter/Projects/GMAT/GMAT/coords"); // hno looks good too
    //    indata.open("/Users/peter/Projects/GMAT/GMAT/coords_amrit"); // ======= principal by itself works
    //  indata.open("/Users/peter/Projects/GMAT/GMAT/coordinates_feb_2013"); // hno looks good too
    //  indata.open("/Users/peter/Projects/GMAT/GMAT/coordinates_april2013");
    // indata.open("/Users/peter/Projects/GMAT/GMAT/coordinate_methyl_formate");
   //  indata.open("/Users/peter/Projects/GMAT/GMAT/BEFarti"); // new code works
   //  indata.open("/Users/peter/Projects/GMAT/GMAT/h2o_coordinates_cooksy"); // new code works
   //  indata.open("/Users/peter/Projects/GMAT/GMAT/test_1D"); // new code works
   
   // indata.open("/Users/peter/Projects/GMAT/GMAT/water_3D_rmin_rplus_small"); // new code works
    
    //  indata.open("/Users/peter/Projects/GMAT/GMAT/BEFartiN"); // new code works
    //  indata.open("/Users/peter/Projects/GMAT/GMAT/BEFartiM"); // bef2 
    //  indata.open("/Users/peter/Projects/GMAT/GMAT/vo"); // this is water with angle and bond
    // indata.open("/Users/peter/Projects/GMAT/GMAT/pps.txt");
    //  indata.open("/Users/peter/Projects/GMAT/GMAT/coordinatesCO2"); // new code works
   //  indata.open("/Users/peter/Projects/GMAT/GMAT/structuredco2");
   //  indata.open("/Users/peter/Projects/GMAT/GMAT/coordinatesH2O"); // aritificial -- cannot test
   //  indata.open("/Users/peter/Projects/GMAT/GMAT/h2o2d"); // new code works
   //  indata.open("/Users/peter/Projects/GMAT/GMAT/coordinatesFHF"); // artificial -- cannot test
    // indata.open("/Users/peter/Projects/GMAT/GMAT/coordintesH2Oartificial");
    // indata.open("/Users/peter/Projects/GMAT/GMAT/waterGaussian"); // new code works
    //  indata.open("/Users/peter/Projects/GMAT/GMAT/vodauholRAD"); // testing case
    // indata.open("/Users/peter/Projects/GMAT/GMAT/bef2"); // this is with Degrees - that does not work
    // indata.open("/Users/peter/Projects/GMAT/GMAT/bef2rad");
    // indata.open("/Users/peter/Projects/GMAT/GMAT/bef2rad2");
    // indata.open("/Users/peter/Projects/GMAT/GMAT/bef2extended.txt"); // -- this seem like a bad geometry too
   //  indata.open("/Users/peter/Projects/GMAT/GMAT/coordArt2D");
    // indata.open("/Users/peter/Projects/GMAT/GMAT/largeMass");
    // indata.open("/Users/peter/Projects/GMAT/GMAT/largemass2D200");
    // indata.open("/Users/peter/Projects/GMAT/GMAT/fhfsys"); //  -- bad geometry structure
    // indata.open("/Users/peter/Projects/GMAT/GMAT/vodauholDEG");
    // indata.open("/Users/peter/Projects/GMAT/GMAT/narrowH2O"); // new code works
    // indata.open("/Users/peter/Projects/GMAT/GMAT/water2Ddenser"); // new code works
    //  indata.open("/Users/peter/Projects/GMAT/GMAT/hno"); // new code works
    // indata.open("/Users/peter/Projects/GMAT/GMAT/mascheck");
    // indata.open("/Users/peter/Projects/GMAT/GMAT/artificial.txt"); // new code works
    // indata.open("/Users/peter/Projects/GMAT/GMAT/coordRandom");
    
    //indata.open("/Users/peter/Projects/GMAT/GMAT/test_3D_NOV2013");
    //indata.open("/Users/peter/Projects/GMAT/GMAT/test_3D_NOV2013_RPhiR");
    //indata.open("/Users/peter/Projects/GMAT/GMAT/test_3D_NOV2013_PhiRR");
    
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
    
        if (dim == 1) // or if dimensionality of internal coordinates == 2
        {
            for (j=0; j<numStructures; j++) // this loop cycles through each geometry
            {
                indata >> internalcoords[j][0];
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
    
    // these variables are defined solely for printing purposes and testing:
       Int cols=3*natom, rows=numStructures;
    
    
    
    
    check_molecule(xgeom, natom, numStructures, molecule_type, coordinates_type, zeros);
    // if the system has less then 4 atoms and dim = 3 there is no need fr numerical estimation:
    if (natom == 3 && dim == 3)
    {
        if (coordinates_type)
        {
            gmatAnalytical(Gmat, internalcoords, xgeom, masses, numStructures, natom);
           // cout << "hello";
            printGmat(Gmat, internalcoords, dim, numStructures);
            exit(0);
        }
        
        else
        {
        translation(numStructures,natom,xgeom,masses);
        Swapper(xgeom, natom, numStructures);
        prerot(xgeom, xreference, natom, numStructures);
        rotation(xreference, numStructures, natom, xgeom, masses);
        playdata(xgeom,natom,numStructures,xreference);
       // printStuff(xgeom, rows, cols);
        jakobAssembly(jacob, internalcoords, xgeom, numStructures, natom, dim);
        checkJac(jacob, dim, numStructures, natom);
        calcGmat(jacob, masses, dim, numStructures, natom, Gmat);
        printGmat(Gmat, internalcoords, dim, numStructures);
        
        
        exit(0);
        }
    }
    
    
    

   if (molecule_type)
    {
        
    
    //  translation(st,natom,xreference,masses);
        translation(numStructures,natom,xgeom,masses);
    //  Swapper(xreference, natom, st);
    Swapper(xgeom, natom, numStructures);
    // cout << "\n swapper: \n";
    // printStuff(xgeom, rows, cols);
    
    prerot(xgeom, xreference, natom, numStructures);
    //  cout << " translation above\n";
    // printStuff(xgeom, rows, cols);

    
    //correct for rotation there are two options still in testing phase:
    rotation(xreference, numStructures, natom, xgeom, masses);
    
    
   //  printStuff(xgeom, rows, cols);
   //  cout << "\n rotation above: \n";
    

 
    
    // if printing is needed uncomment the two lines below:
    // printStuff(xgeom, rows, cols);
    
    // checking if the rotational condition is satisfied (commented for now):
    // rottest(xgeom, xreference, masses, numStructures, natom);
    
    // print the reference structure:
    
    // ==== By observing the data I noticed that some of the structures have an exchanged sign ====== //
    // ==== this makes the gradients at many points extremly large and inaccurate =================== //
    // ===== the routine below takes care of the signs and plays around with the data =============== //
     playdata(xgeom,natom,numStructures,xreference);
     // cout << " nakonci zrotovane: \n";
     // printStuff(xgeom, rows, cols);
    }

    
    

   else
    {
        //   translation(st,natom,xreference,masses);
        translation(numStructures,natom,xgeom,masses);
        rot_princ_axis(numStructures, natom, xgeom, masses, zeros);
     
        // printStuff(xgeom, rows, cols);
        //  cout << "\n";
    }



    // atom_exchange(xgeom, natom, numStructures, masses);
    //    printStuff(xgeom, rows, cols);    
    
    /* //  ===== TESTING THE GEOMETRY ===========================
    // for testing purposes I am printing intermediate results
    ofstream myfile;
    char filename[100]="//Users//peter//Projects//GMAT//GMAT//out.txt";
    myfile.open (filename);
    for (i=0;i<numStructures;i++)
    {
        for (j=0; j<natom; j++)
        {
            l=3*j;
            myfile << xgeom[i][0+l] << " " << xgeom[i][1+l] << " " << xgeom[i][2+l] << "\n";
        }
        myfile << "\n";
       
    }
    myfile.close();
    // ============ TESTING FINISHED - PROGRAM WORKS =================================
     */ 

	//    new routine without using the reference:
//    trans_rot(xgeom, internalcoords, masses, numStructures, dim, natom);
    
//    printStuff(xgeom, rows, cols);
   
    
//    rottest(xgeom, xreference, masses, numStructures, natom);
    
    
    // Assemble the Jakobian Matrix
    jakobAssembly(jacob, internalcoords, xgeom, numStructures, natom, dim);

    // checkJac(jacob, dim, numStructures, natom);
    
    // If testing is needed print out the Jakobian Matrix, uncomment the two lines below:
    // Int rowss = numStructures, colss = natom*3*dim;
    // printStuff(jacob, rowss, colss);
    
    
    
    // Calculate G-Matrix
    calcGmat(jacob, masses, dim, numStructures, natom, Gmat);
 
   
    // Check ratio of masses
   //  checkmasses(Gmat,dim,numStructures,natom,masses);			
 
    // print the Gmatrix:
     printGmat(Gmat, internalcoords, dim, numStructures);
    

    // Program finished
}
