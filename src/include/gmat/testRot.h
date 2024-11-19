//
//  testRot.h
//  GMAT
//
//  Created by PETER ZAJAC on 11/9/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

// here I will test rotation:

void rottest(MatDoub &xgeoms, MatDoub &ref, VecDoub &mas, Int &str, Int &at)
{
    int i,j,k,l;
    MatDoub temp(at,3);
    MatDoub temRef(at,3);
    VecDoub res(3);
    double totalmass;
    
    for (j=0; j<at; j++) 
    {
        l=3*j;
        temRef[j][0] = ref[0][0+l];
        temRef[j][1] = ref[0][1+l];
        temRef[j][2] = ref[0][2+l];
        
    }
    
    
    for (i=0;i<str;i++)
    {
        totalmass = 0.0;
        for (j=0; j<at; j++) // create a temporary geometry storage for the given structure
        {
            l=3*j;    
            temp[j][0] = xgeoms[i][0+l];
            temp[j][1] = xgeoms[i][1+l];
            temp[j][2] = xgeoms[i][2+l];
        }
        
        // initialize res to 0
        for (k=0; k<3; k++) 
        {
            res[k]=0;
        }

        // compute res as a vector product weighted by the masses
        for (k=0; k<at; k++) 
        {
            res[0] = res[0] + mas[k]*( temRef[k][1]*temp[k][2] - temRef[k][2]*temp[k][1]);
            res[1] = res[1] + mas[k]*( temRef[k][2]*temp[k][0] - temRef[k][0]*temp[k][2]);
            res[2] = res[2] + mas[k]*( temRef[k][0]*temp[k][1] - temRef[k][1]*temp[k][0]);
        }
        
        totalmass = res[0] + res[1] + res[2];
        cout << res[0] << " " << res[1] << " " << res[2] << "     " << totalmass << "\n";
    } // end the loop for Structures
}

