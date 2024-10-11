//
//  jakobAssembly.h
//  GMAT
//
//  Created by PETER ZAJAC on 10/7/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

void jakobAssembly(MatDoub &jakob, MatDoub &internals, MatDoub &geometries, Int &numstructures, Int &numatoms, Int &dim)
{
   /*
    jakob matrix is (numstructs , dim*3*natom)
    internals (numstructs , dim)
    geometries (numstructs , 3*natom )
    */
    
    int i,j,k;
    Doub cutoff=1.0,diff=0.0;
    
    // MatDoub xtemp(numstructures,dim); // this matrix will hold the temporary grid
    VecDoub ytemp(numstructures); // this will hold the functional values for temp grid 
    VecDoub pointTemp(dim);
    VecDoub pointTempXplus(dim);
    VecDoub pointTempX2plus(dim);
    VecDoub pointTempXminus(dim);
    VecDoub pointTempX2minus(dim);
    VecDoub pointTempYplus(dim);
    VecDoub pointTempY2plus(dim);
    VecDoub pointTempYminus(dim);
    VecDoub pointTempY2minus(dim);
    VecDoub pointTempZplus(dim);
    VecDoub pointTempZ2plus(dim);
    VecDoub pointTempZminus(dim);
    VecDoub pointTempZ2minus(dim);

    
    VecDoub deriAlongCoord(dim);
   // Doub delta=0.001;
    VecDoub deltas(dim);
  
    if (dim == 1)
    {
        //findCut(internals, numstructures, dim, cutoff, deltas);
        
        for (i=0; i<3*numatoms; i++) // main loop over all coordinates
        {
            for (j=0; j<numstructures; j++)
            {
                ytemp[j] = geometries[j][i];
                
            }
            
            for (j=0; j<numstructures; j++)
            {
                
                pointTemp[0] = internals[j][0];
                
                imls2(internals, ytemp, pointTemp, numstructures, dim, deriAlongCoord);
                
                jakob[j][i] = deriAlongCoord[0];
                
            }
        }
    }
    
    
    if (dim == 2)
    {
        //findCut(internals, numstructures, dim, cutoff, deltas);
   
        for (i=0; i<3*numatoms; i++) // main loop over all coordinates
            {
                for (j=0; j<numstructures; j++) 
                    {
                        ytemp[j] = geometries[j][i];
          
                    }   
        
                for (j=0; j<numstructures; j++) 
                {
            
                    for (k=0; k<dim; k++) 
                        {
                            pointTemp[k] = internals[j][k];
                        }
                    
                    imls2(internals, ytemp, pointTemp, numstructures, dim, deriAlongCoord);

            /*  // tu zakomentovat      
                    for (k=0; k<dim; k++) 
                    {
                        if (abs(deriAlongCoord[k]) < 0.5) {deriAlongCoord[k] = 0.0;}
                    }
              // tu zakomentovat
            */
                    for (k=0; k<dim; k++) 
                        {
                            jakob[j][dim*i + k] = deriAlongCoord[k];
                        }
        
                }
            }
    }

    
    if (dim == 3) 
    {
       // findCut(internals, numstructures, dim, cutoff, deltas);
        
        for (i=0; i<3*numatoms; i++) // main loop over all coordinates
        {
            for (j=0; j<numstructures; j++) 
            {
                ytemp[j] = geometries[j][i];
            }   
         /*
            // check whether the difference is significant:
            diff=0.0;
            for (j=1; j<numstructures; j++)
            {
                diff += abs(ytemp[j-1] - ytemp[j]);
            }
            
            diff=diff/(numstructures-1);
            
            if (diff < 0.03) // coordinate is not important
            {
                for (j=0; j<numstructures; j++)
                {
                    ytemp[j] = 0.0;
                }
            }
          */  
            for (j=0; j<numstructures; j++) 
            {

                
                for (k=0; k<dim; k++) 
                {
                    pointTemp[k] = internals[j][k];
                }
                
                imls2(internals, ytemp, pointTemp, numstructures, dim, deriAlongCoord); 
               
              
           /*     if (numatoms == 3 && dim == 3 && j>0)
                {
                    for (k=0; k<dim; k++) 
                    {
                        if (abs(deriAlongCoord[k]) > 1.0) {deriAlongCoord[k] = jakob[j-1][dim*i+k];}
                    }
                
                }
             */ 
 
                
                for (k=0; k<dim; k++) 
                {
                    jakob[j][dim*i + k] = deriAlongCoord[k];
                }
                
                
            }
            
        }
       // cout << "atom done \n";
    }
    
    // end of routine
}


void checkJac(MatDoub &jakobian, Int &dim, Int &numstructs, Int &natom)
{
    
    int i,j;
    double a,b;
    
    
    for (i=0;i<3*dim*natom;i++)
    {
        for (j=0; j<numstructs; j++)
        {
            a += jakobian[j][i];
        }
        
        a = a / numstructs;
        b = 15*a;
        
        for (j=0; j<numstructs; j++)
        {
            if (abs(jakobian[j][i]) > b)
            {
                jakobian[j][i] = 7*a;
            }
        }
        
        
    }
    
    
}



void calcGmat(MatDoub &jakobian, VecDoub &masses, Int &dim, Int &numstructs, Int &natom, MatDoub &Gmat)

{
    // Gmat is size structures * (dim^2)
    Int p,r,s,i,j,k,l;
    Doub x,y,z,det;
    MatDoub ginv(dim,dim);
    MatDoub gmatTemp(dim,dim); 
    MatDoub rhs(dim,dim);
    l= dim * 3;
    
    
    for (i=0; i<dim; i++)
    {
        for (j=0; j<dim; j++) 
        {
            rhs[i][j] = 0.0;
        }
    }
    
    for (i=0; i<dim; i++) {
        rhs[i][i] = 1.0;
    }
    
    for (p=0; p<numstructs; p++)  // loop through all the structures
    {
        // == SET GINV TO ZERO ======
        for (r=0; r<dim; r++) 
        {
            for (s=0; s<dim; s++) 
            {
                ginv[r][s] = 0.0;
            }
        }
        
        

        // === FINISH SETTING GINV TO ZERO
        
        
        // ASSEMBLE GINV using the JACOBIAN MATRIX AND THE MASSES
        for (r=0; r<dim; r++) 
        {
            for (s=0; s<dim; s++) 
            {
                for (i=0; i<natom; i++) 
                {
                    for (k=0; k<3; k++) 
                    {
                        ginv[r][s] = ginv[r][s] + masses[i]*jakobian[p][i*l + r + k*dim]*jakobian[p][i*l + s + k*dim];
       //                 cout << masses[i] << " times " << jakobian[p][i*l + r + k*dim] << " times " << jakobian[p][i*l + s + k*dim] << " plus ";
                    }   
                   
                }
           //      cout << " equals:  " << ginv[r][s] << " " << r << " " << s << "\n";
                
            }
          //  cout << "\n";
            
        } // ASSEBMLY FINISHED
     
        /*
     
        x = ginv[0][0]*(ginv[1][1]*ginv[2][2] - ginv[1][2]*ginv[2][1]);
        y = ginv[0][1]*(ginv[1][0]*ginv[2][2] - ginv[1][2]*ginv[2][0]);
        z = ginv[0][2]*(ginv[1][0]*ginv[2][1] - ginv[1][1]*ginv[2][0]);
        det = x - y + z;
        
        
        cout << det << "\n";
        
        if (det < 0.1)
        {
            Doub max=ginv[0][0];
            for (r=0; r<dim; r++)
            {
                for (s=0; s<dim; s++)
                {
                    if (ginv[r][s] > max)
                    {
                        max= ginv[r][s];
                    }
                }
            }
            
            for (r=0; r<dim; r++)
            {
                for (s=0; s<dim; s++)
                {
                    if (ginv[r][s] == max)
                    {
                        ginv[r][s] = ginv[r][s] - 0.2*ginv[r][s];
                    }
                }
            }
        }
        */
	 
        if (dim > 1)
        {
           SVD decompose(ginv);
            decompose.solve(rhs, gmatTemp);
            
            // write the temporary gmatrix into the main big Gmatrix structure
            for (r=0; r<dim; r++)
            {
                for (s=0; s<dim; s++)
                {
                    Gmat[p][r*dim + s] = gmatTemp[r][s];
                }
            }
        }
        else
        {
            
            Gmat[p][0] = 1/ginv[0][0];
            
        }
        // inversion finished
        
        
        

        
    } // structure loop finished
    
}





