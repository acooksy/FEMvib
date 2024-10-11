//
//  checkinversion.h
//  GMAT
//
//  Created by PETER ZAJAC on 10/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//


// MatDoub &refx, Int &numberofStructures, Int &numberofAtoms, MatDoub &xgeoms, VecDoub &mass
void checkinv(MatDoub &ref, MatDoub &geom, Int &natom)
{
    // ref is numatoms X 3
    // geom is numatoms X 3
    VecInt imax(3);
    VecDoub xmax(3);    
    Doub controla;
    int i,j,k;
    
    for (i=0; i<3; i++)
    {
        imax[i] = 0;
        xmax[i] = abs(ref[0][i]);
        for (j=1; j<natom; j++)
        {
            if (abs(ref[j][i]) > xmax[i]) 
            {
                imax[i] = j;
                xmax[i] = abs(ref[j][i]);
            }
        }
        
    }
    
    for (i=0;i<3; i++) 
    {
        k=imax[i];
        controla = geom[k][i]*ref[k][i];
        if (controla < 0.0) 
        {
            for (j=0; j<natom; j++) 
            {
                geom[j][i] = -geom[j][i];
            }
        }
                                        
        
    }
    
    
}

void normalize(VecDoub &n)
{
    
    double p;
    double pi=3.1415926535897932384626433;
    int i,j;
    
   // cout << n[0] << " " << n[1] << " " << n[2] << "\n";
    p=0.0;
    for (i=0;i<3;i++)
    {
        p = p + n[i]*n[i];
    }
    p = sqrt(p);
    
    for (j=0; j<3; j++) 
    {
        n[j] = n[j]/p;
    }
   // cout << p << n[0] << " " << n[1] << " " << n[2] << "\n";
    
  //  cout << n[3] << "\n";
    while (abs(n[3]) > pi) 
    {
        if (n[3] < 0.0) 
        {
            n[3] = n[3] + pi;
        }
        if (n[3] > pi)
        {
            n[3] = n[3] - pi;
        }
    }
    
    if (n[3] < 0.0)
    {
        n[3] = n[3] + pi;
    }
   // cout << n[3] << "\n\n";
    
    
}
