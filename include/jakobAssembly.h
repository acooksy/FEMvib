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
    Doub cutoff;
    
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
    
    Doub xp, xm ,yp, ym, zp, zm;
    Doub x2p, x2m ,y2p, y2m, z2p, z2m;
    
    VecDoub deriAlongCoord(dim);
   // Doub delta=0.001;
    VecDoub deltas(dim);
    
    if (dim == 2) 
    {
        findCut(internals, numstructures, dim, cutoff, deltas);
   
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
                    
                    imls2(internals, ytemp, pointTemp, numstructures, dim, cutoff, deriAlongCoord);            
                    
             /*       for (k=0; k<dim; k++) 
                    {
                        if (abs(deriAlongCoord[k]) < 0.0001) {deriAlongCoord[k] = 0.0;}
                    }
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
        findCut(internals, numstructures, dim, cutoff, deltas);
        
        for (i=0; i<3*numatoms; i++) // main loop over all coordinates
        {
            for (j=0; j<numstructures; j++) 
            {
                ytemp[j] = geometries[j][i];
            }   
            
        //    Shep_interp shepard(internals,ytemp);
        //    Powvargram vgram(internals,ytemp);
        //    Krig<Powvargram> kr(internals,ytemp,vgram);
            
            for (j=0; j<numstructures; j++) 
            {
     /*           for (k=0; k<dim; k++) // if using imls this line starts the comment
                {
                    if (k == 0)
                    {
                        pointTempXplus[k] = internals[j][k] + deltas[0];
                        pointTempX2plus[k] = internals[j][k] + 2*deltas[0];
                        pointTempXminus[k] = internals[j][k] - deltas[0];
                        pointTempXminus[k] = internals[j][k] - 2*deltas[0];
                    }
                    if (k == 1)
                    {
                        pointTempYplus[k] = internals[j][k] + deltas[1];
                        pointTempY2plus[k] = internals[j][k] + 2*deltas[1];
                        pointTempYminus[k] = internals[j][k] - deltas[1];
                        pointTempY2minus[k] = internals[j][k] - 2*deltas[1];
                    }
                    if (k == 2)
                    {
                        pointTempZplus[k] = internals[j][k] + deltas[2];
                        pointTempZ2plus[k] = internals[j][k] + 2*deltas[2];
                        pointTempZminus[k] = internals[j][k] - deltas[2];
                        pointTempZ2minus[k] = internals[j][k] - 2*deltas[2];
                    }
                    if (k !=2 )
                    {
                    pointTempZplus[k] = internals[j][k];
                    pointTempZminus[k] = internals[j][k];
                    pointTempZ2plus[k] = internals[j][k];
                    pointTempZ2minus[k] = internals[j][k];
                    }
                    if (k != 1) 
                    {
                    pointTempYplus[k] = internals[j][k];
                    pointTempYminus[k] = internals[j][k];                    
                    pointTempY2plus[k] = internals[j][k];
                    pointTempY2minus[k] = internals[j][k];                    
                        
                    }
                    if (k != 0)
                    {
                    pointTempXplus[k] = internals[j][k];
                    pointTempXminus[k] = internals[j][k];
                    pointTempX2plus[k] = internals[j][k];
                    pointTempX2minus[k] = internals[j][k];    
                    }
                }
                // imls2(internals, ytemp, pointTemp, numstructures, dim, cutoff, deriAlongCoord);
                xp = kr.interp(pointTempXplus);
                xm = kr.interp(pointTempXminus);
                yp = kr.interp(pointTempYplus);
                ym = kr.interp(pointTempYminus);
                zp = kr.interp(pointTempZplus);
                zm = kr.interp(pointTempZminus);
                x2p = kr.interp(pointTempX2plus);
                x2m = kr.interp(pointTempX2minus);
                y2p = kr.interp(pointTempY2plus);
                y2m = kr.interp(pointTempY2minus);
                z2p = kr.interp(pointTempZ2plus);
                z2m = kr.interp(pointTempZ2minus);
                deriAlongCoord[0] = (-x2p + 8*xp - 8*xm + x2m)/(12*deltas[0]);
                deriAlongCoord[1] = (-y2p + 8*yp - 8*ym + y2m)/(12*deltas[1]);
                deriAlongCoord[2] = (-z2p + 8*zp - 8*zm + z2p)/(12*deltas[2]);  // comment ends here
          */
           /*     deriAlongCoord[0] = (xp - internals[j][0])/0.1;
                deriAlongCoord[1] = (yp - internals[j][1])/0.1;
                deriAlongCoord[2] = (zp - internals[j][2])/0.1;
            */   
                
                for (k=0; k<dim; k++) 
                {
                    pointTemp[k] = internals[j][k];
                }
                
                imls2(internals, ytemp, pointTemp, numstructures, dim, cutoff, deriAlongCoord); 
               
             /* 
                if (numatoms == 3 && dim == 3)
                {
                    for (k=0; k<dim; k++) 
                    {
                        if (abs(deriAlongCoord[k]) < 0.01) {deriAlongCoord[k] = 0.0;}
                    }
                
                }
              */
              //  cout << xp << " " << xm << " " << yp << " " << ym << " " << zp << " " << zm << "\n"; 
                
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


void calcGmat(MatDoub &jakobian, VecDoub &masses, Int &dim, Int &numstructs, Int &natom, MatDoub &Gmat)

{
    // Gmat is size structures * (dim^2)
    Int p,r,s,i,j,k,l;
    
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
                     //   cout << masses[i] << " times " << jakobian[p][i*l + r + k*dim] << " times " << jakobian[p][i*l + s + k*dim] << " plus \n";
                    }   
                   
                }
                // cout << "\n equals:  " << ginv[r][s] << "\n\n";
                
            }
            
        } // ASSEBMLY FINISHED
     
        
    // invert ginv to obtain Gmat
      /*  LUdcmp alu(ginv);
        alu.inverse(gmatTemp);
       */ 
           SVD decompose(ginv);
            decompose.solve(rhs, gmatTemp);
        
      /*  
        for (int h=0;h<dim;h++)
        {
            for (int u=0; u<dim; u++) 
            {
                cout << ginv[h][u] << " ";
            }
            cout << " \n";
        }
        cout << "\n";
       */ 
         
        /* cout << ginv[0][0] << " " << ginv[0][1] << " " <<; 
            cout << ginv[1][0] << " " << ginv[1][1] << "\n";
            cout << "\n";
       */
        // inversion finished 
        
        
        
        // write the temporary gmatrix into the main big Gmatrix structure
        for (r=0; r<dim; r++)
        {
            for (s=0; s<dim; s++) 
            {
                Gmat[p][r*dim + s] = gmatTemp[r][s];
            }
        }
        
    }
    
}





