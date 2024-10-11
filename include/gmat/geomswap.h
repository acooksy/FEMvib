//
//  geomswap.h
//  GMAT
//
//  Created by PETER ZAJAC on 11/2/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

// this function will swap all geometries if needed - not all the QCH programs are smart


void Swapper(MatDoub &xgeoms, Int &at, Int &str)
{
    // geoms( str , at*3 ) 
    // str - number of structures
    // at number of atoms
    
    MatDoub temp(at,3);
    int i,j,l,p;
    Doub sumX, sumY, sumZ, t;
    
    
    for (i=0;i<str;i++)
    {
        // set the sums along each coordinate equal to zero
        sumX = 0.0;
        sumY = 0.0;
        sumZ = 0.0;
        
        for (j=0; j<at; j++) // create a temporary geometry storage for the given structure
        {
            l=3*j;    
            temp[j][0] = xgeoms[i][0+l]; // x coordinates
            temp[j][1] = xgeoms[i][1+l]; // y coordinates
            temp[j][2] = xgeoms[i][2+l]; // z coordinates
        }

        // determine the sums for each coordinate:
        for (p=0; p<at; p++) 
        {
            sumX = sumX + abs(temp[p][0]);
            sumY = sumY + abs(temp[p][1]);
            sumZ = sumZ + abs(temp[p][2]);
        }
   
       // if statements for vib_rot_spectro
//	cout << i + 1 << " ";
//	if (sumZ == 0 && sumY == 0 && sumX == 0){cout << "problem\n"; exit(1);}
//	if (sumZ == 0 && sumY == 0 && sumX != 0){cout << "1\n";} // xz is 1
//	if (sumZ == 0 && sumY != 0 && sumX == 0){cout << "2\n";} // yz is 2
//	if (sumZ == 0 && sumY != 0 && sumX != 0){cout << "1\n";} // xz is 1
//	if (sumZ != 0 && sumY == 0 && sumX == 0){cout << "0\n";} // nothing is 0
//	if (sumZ != 0 && sumY != 0 && sumX == 0){cout << "0\n";}
//	if (sumZ != 0 && sumY != 0 && sumX != 0){cout << "0\n";}
//	if (sumZ != 0 && sumY == 0 && sumX != 0){cout << "3\n";} // xy is 3




 
        // THE IF STATEMENTS:
        if (sumZ == 0) 
        {
            if (sumY == 0) // z and y are 0
            {
                // swap x and z
                for (p=0; p<at; p++) 
                {
                    t = temp[p][0];
                    temp[p][0] = 0;
                    temp[p][2] = t;
                }
            }
            if (sumX == 0) // z and x are 0
            {
                // swap y and z
                for (p=0; p<at; p++) 
                {
                    t = temp[p][1];
                    temp[p][1] = 0;
                    temp[p][2] = t;
                }
            }
            else // only z is zero 
            {
                // swap x and z, keep y
                for (p=0; p<at; p++) 
                {
                    t = temp[p][0];
                    temp[p][0] = 0;
                    temp[p][2] = t;
                } 
            }
            
         } // end of if z == 0 statement
        
        else // beginning of z != 0 statement
        {
            if (sumY == 0 && sumX != 0)
            {
                // swap x and y
                for (p=0; p<at; p++) 
                {
                    t = temp[p][0];
                    temp[p][0] = 0;
                    temp[p][1] = t;
                } 
            }
        }

    
    
        for (j=0; j<at; j++)
        {
            l=3*j;    
            xgeoms[i][0+l] = temp[j][0];
            xgeoms[i][1+l] = temp[j][1];
            xgeoms[i][2+l] = temp[j][2];
        }    
        
        
    } // end of structure loop
    
}


void prerot(MatDoub &xg, MatDoub &ref, Int &nat, Int &nstr)
{
    
    // xg size is (nstr , 3*nat)
    int i,j,k,l,ind;
    MatDoub m1(nat,3),m2(nat,3),m3(nat,3); // these will hold the rotated structures
    MatDoub temp(nat,3),comp(nat,3); // comp will store the reference structur; temp individual structures
    Doub min;
    VecDoub sum(8);
  /*
    for (j=0; j<nat; j++) // create a temporary geometry storage for the reference
    {
        l=3*j;    
        comp[j][0] = ref[0][0+l];
        comp[j][1] = ref[0][1+l];
        comp[j][2] = ref[0][2+l];
    }
    */
    // cycle through all the structures and check the signs
    for (i=1;i<nstr; i++)
    {
        for (j=0; j<nat; j++) // create a temporary geometry storage for the reference
        {
            l=3*j;
            comp[j][0] = xg[i-1][0+l];
            comp[j][1] = xg[i-1][1+l];
            comp[j][2] = xg[i-1][2+l];
        }
        
        for (j=0; j<nat; j++) // create a temporary geometry storage for the previous structure
        {
            l=3*j;
            temp[j][0] = xg[i][0+l];
            temp[j][1] = xg[i][1+l];
            temp[j][2] = xg[i][2+l];
        }
        
        for (k=0;k<8;k++) {sum[k]=0.0; } // initialize the sums
        
        // calculate sum corresponding to the current structure:
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[0] = sum[0] + abs(temp[a][c] - comp[a][c]);
            }
        }
        
        // calculate sum corresponding to the geometry multiplied by "-1"
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[1] = sum[1] + abs(-temp[a][c] - comp[a][c]);
            }
        }
        
        // === Z AXIS =========
        // rotate the geometry of i-th structure along z-axis 180 deg
        
        for (k=0; k<nat; k++ ) 
        {
            m1[k][0] = -temp[k][0]; 
            m1[k][1] = -temp[k][1];
            m1[k][2] =  temp[k][2];
        } // rotation along z-axis of i-th structure finished
        
        // calculate sum corresponding to the rotated geometry 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[2] = sum[2] + abs(m1[a][c] - comp[a][c]);
            }
        }    
        
        // calculate sum corresponding to the rotated geometry multiplied by "-1" 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[3] = sum[3] + abs(-m1[a][c] - comp[a][c]);
            }
        }        
        
        // ===== X AXIS ======
        // rotate the geometry of i-th structure along x-axis 180 deg
        
        for (k=0; k<nat; k++ ) 
        {
            m2[k][0] =  temp[k][0]; 
            m2[k][1] = -temp[k][1];
            m2[k][2] = -temp[k][2];// + 2*temp[k][2];
        } // rotation along z-axis of i-th structure finished
        
        // calculate sum corresponding to the rotated geometry 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[4] = sum[4] + abs(m2[a][c] - comp[a][c]);
            }
        }    
        
        // calculate sum corresponding to the rotated geometry multiplied by "-1" 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[5] = sum[5] + abs(-m2[a][c] - comp[a][c]);
            }
        }  
        
        // ===== Y AXIS ======
        // rotate the geometry of i-th structure along y-axis 180 deg
        
        for (k=0; k<nat; k++ ) 
        {
            m3[k][0] = -temp[k][0]; // + 2*temp[k][0]; 
            m3[k][1] =  temp[k][1];
            m3[k][2] = -temp[k][2];// + 2*temp[k][2];
        } // rotation along z-axis of i-th structure finished
        
        // calculate sum corresponding to the rotated geometry 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[6] = sum[6] + abs(m3[a][c] - comp[a][c]);
            }
        }    
        
        // calculate sum corresponding to the rotated geometry multiplied by "-1" 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[7] = sum[7] + abs(-m3[a][c] - comp[a][c]);
            }
        } 
        
        // find out which sum is the smallest
        ind=0;
        min = sum[0];
        for (k=1;k<8;k++)
        {
            if (sum[k] < min)
            {
                ind = k;
                min = sum[k];
            }
        }
        
   //     cout << min << "\n";
        
        ind = ind + 1;
       
 //	cout << i+1 << " "; 
        switch (ind) 
        {
            case 1:
	//	cout << "1\n";
                break;
                
            case 2:
	
	//	cout << "2\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = -temp[a][c];     
                    }
                }
                break;
                
            case 3:
	//	cout << "3\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = m1[a][c];     
                    }
                }
                break;
                
            case 4:
	//	cout << "4\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = -m1[a][c];     
                    }
                }
                break;
                
            case 5:
	//	cout << "5\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = m2[a][c];     
                    }
                }
                break;    
                
            case 6:
	//	cout << "6\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = -m2[a][c];     
                    }
                }
                break;
                
            case 7:
	//	cout << "7\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = m3[a][c];     
                    }
                }
                break;
                
            case 8:
	//	cout << "8\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = -m3[a][c];     
                    }
                }
                break;
                
            default: cout << "index has not been determined - error in structures\n";
                break;
        }
        
        
        // return the geometry in the new form to the main geometry matrix
        for (k=0; k<nat; k++) 
        {
            l=k*3;
            xg[i][0+l] = temp[k][0];
            xg[i][1+l] = temp[k][1];
            xg[i][2+l] = temp[k][2];
            
        }
        
        
        
        
    } // and of structure loop
    
    
    
    
    
}


