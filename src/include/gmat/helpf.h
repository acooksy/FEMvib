//
//  helpf.h
//  ILMS
//
//  Created by PETER ZAJAC on 9/28/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

void check_molecule(MatDoub &xq, Int &nat, Int &nstr, Int &decision, Int &coordinates, long int &zs)
{
        int i,j,k,l,index=1,h,m,n,o;
        VecInt coefs_for_atoms(nat);
        MatDoub temp(nat,3);
        MatDoub comp(nat,3);
    zs=0;
    m=0;
    for (i=0; i<nstr; i++)
    {
    
        for (j=0; j<nat; j++) // create a temporary geometry storage for the given structure
        {
            l=3*j;
            temp[j][0] = xq[i][0+l];
            temp[j][1] = xq[i][1+l];
            temp[j][2] = xq[i][2+l];
        }
        
        // Determine if the structure has zeros along the coordinates for each atom:
        k=0;
        for (n=0; n<nat; n++)
        {
            h=0;
            for (j=0; j<3; j++)
            {
                if (abs(temp[n][j]) < 0.0001)
                {
                    h++;
                }
            }
        
            if (h>0)
            {
                k++;
                zs++;
            }
        }
        
        if (k >= nat)
        {
            m++;
        }
    }
  //  cout << m;
    if (m==nstr)
    {
        decision = 1;
    }
    
    
    
    
    for (i=0; i<3*nat; i++)
    {
        for (j=0; j<nstr; j++)
        {
            if (xq[j][i] < 0)
            {
                coordinates = 1;
            }
        }
    }
    
}

void exchange(MatDoub &xg, Int &nat, MatDoub &ref, Int print) //, MatDoub &Internals)
{
    // xg size is (nstr , 3*nat)
    int i,j,k,l,ind;
    MatDoub m1(nat,3),m2(nat,3),m3(nat,3); // these will hold the rotated structures
    MatDoub m90x(nat,3),m90y(nat,3), m90z(nat,3);
    MatDoub temp(nat,3),comp(nat,3);
    Doub min,min90;
    VecDoub sum(6);
    
    for (j=0; j<nat; j++) // create a temporary geometry storage for the reference
    {
        l=3*j;
        comp[j][0] = ref[0][0+l];
        comp[j][1] = ref[0][1+l];
        comp[j][2] = ref[0][2+l];
    }
    
    
    for (j=0; j<nat; j++) // create a temporary geometry storage for the previous structure
    {
        l=3*j;
        temp[j][0] = xg[0][0+l];
        temp[j][1] = xg[0][1+l];
        temp[j][2] = xg[0][2+l];
    }
    
    for (k=0;k<6;k++) {sum[k]=0.0; } // initialize the sums
    
    if (print)
    {
        for (j=0; j<nat; j++) // create a temporary geometry storage for the previous structure
        {
            cout << temp[j][0] << " " << temp[j][1] <<  " " << temp[j][2] << " ";
        }
        cout << "\n";
        
    }
    
    
    
    
    // calculate sum corresponding to the current structure:
    for (int a=0; a<nat; a++) // loop through atoms
    {
        sum[0] = sum[0] + abs(temp[a][0] - comp[a][0]);
    }
    
    for (int a=0; a<nat; a++) // loop through atoms
    {
        sum[1] = sum[1] + abs(temp[a][0] + comp[a][0]);
    }
    
    if (sum[1] < sum[0] )
    {
        for (int a=0;a<nat;a++)
        {
            temp[a][0] = -temp[a][0];
        }
    }
    
    
    for (int a=0; a<nat; a++) // loop through atoms
    {
        sum[2] = sum[2] + abs(temp[a][1] - comp[a][1]);
    }
    
    for (int a=0; a<nat; a++) // loop through atoms
    {
        sum[3] = sum[3] + abs(temp[a][1] + comp[a][1]);
    }
    
    if (sum[3] < sum[2] )
    {
        for (int a=0;a<nat;a++)
        {
            temp[a][1] = -temp[a][1];
        }
    }
    
    for (int a=0; a<nat; a++) // loop through atoms
    {
        sum[4] = sum[4] + abs(temp[a][2] - comp[a][2]);
    }
    
    for (int a=0; a<nat; a++) // loop through atoms
    {
        sum[5] = sum[5] + abs(temp[a][2] + comp[a][2]);
    }
    
    if (sum[5] < sum[4] )
    {
        for (int a=0;a<nat;a++)
        {
            temp[a][2] = -temp[a][2];
        }
    }
    
    if (print)
    {
        for (int h=0;h<6;h++)
        {
            cout << sum[h] << " ";
        }
        cout << "\n";
        for (j=0; j<nat; j++) // create a temporary geometry storage for the previous structure
        {
            cout << comp[j][0] << " " << comp[j][1] <<  " " << comp[j][2] << " ";
        }
        cout << "\n";
        for (j=0; j<nat; j++) // create a temporary geometry storage for the previous structure
        {
            cout << temp[j][0] << " " << temp[j][1] <<  " " << temp[j][2] << " ";
        }
        cout << "\n";
    }
    
    
    // return the geometry in the new form to the main geometry matrix
    for (k=0; k<nat; k++)
    {
        l=k*3;
        xg[0][0+l] = temp[k][0];
        xg[0][1+l] = temp[k][1];
        xg[0][2+l] = temp[k][2];
        
    }
    
    // end of function
}



void check_alignment(MatDoub &pred, Int &nat, MatDoub &ref)
{
    int i,j,k,l,index=1,h;
    VecInt coefs_for_atoms(nat);
    MatDoub temp(nat,3);
    MatDoub comp(nat,3);
    double aa,bb,cc;
    double min;
    int ind;
    
    
    for (j=0; j<nat; j++) // create a temporary geometry storage for the given structure
    {
        l=3*j;
        temp[j][0] = pred[0][0+l];
        temp[j][1] = pred[0][1+l];
        temp[j][2] = pred[0][2+l];
    }
    
    for (j=0; j<nat; j++) // create a temporary geometry storage for the previous structure
    {
        l=3*j;
        comp[j][0] = ref[0][0+l];
        comp[j][1] = ref[0][1+l];
        comp[j][2] = ref[0][2+l];
    }
    
    // Determine if the structure has zeros along the coordinates for each atom:
    k=0;
    for (i=0; i<nat; i++)
    {
        h=0;
        for (j=0; j<3; j++)
        {
            if (temp[i][j] == 0.0)
            {
                h++;
            }
        }
        if (h>0)
        {
            k++;
        }
        
    }
    
    // MatDoub test_xyz(nat,3);
    MatDoub test_xzy(nat,3);
    MatDoub test_yxz(nat,3);
    MatDoub test_yzx(nat,3);
    MatDoub test_zxy(nat,3);
    MatDoub test_zyx(nat,3);
    VecDoub sum(6);
    
    for (i=0; i<6; i++)
    {
        sum[i] = 0.0;
    }
    
    // play around with the coordinates to find the best alignment.
    if (k >= nat)
    {
        // test all 6 posibilities of interchanging the coordinates
        
        // sum[0] holds the unchanged coordinates XYZ
        for (i=0;i<nat;i++)
        {
            sum[0] = sum[0] + abs(abs(comp[i][0])-abs(temp[i][0]));
            sum[0] = sum[0] + abs(abs(comp[i][1])-abs(temp[i][1]));
            sum[0] = sum[0] + abs(abs(comp[i][2])-abs(temp[i][2]));
            
        }
        
        // sum[1] holds the XZY
        for (i=0;i<nat;i++)
        {
            sum[1] = sum[1] + abs(abs(comp[i][0])-abs(temp[i][0]));
            sum[1] = sum[1] + abs(abs(comp[i][1])-abs(temp[i][2]));
            sum[1] = sum[1] + abs(abs(comp[i][2])-abs(temp[i][1]));
        }
        
        // sum[2] holds the YXZ
        for (i=0;i<nat;i++)
        {
            sum[2] = sum[2] + abs(abs(comp[i][0])-abs(temp[i][1]));
            sum[2] = sum[2] + abs(abs(comp[i][1])-abs(temp[i][0]));
            sum[2] = sum[2] + abs(abs(comp[i][2])-abs(temp[i][2]));
        }
        
        // sum[3] holds the YZX
        for (i=0;i<nat;i++)
        {
            sum[3] = sum[3] + abs(abs(comp[i][0])-abs(temp[i][1]));
            sum[3] = sum[3] + abs(abs(comp[i][1])-abs(temp[i][2]));
            sum[3] = sum[3] + abs(abs(comp[i][2])-abs(temp[i][0]));
        }
        
        // sum[4] holds the ZXY
        for (i=0;i<nat;i++)
        {
            sum[4] = sum[4] + abs(abs(comp[i][0])-abs(temp[i][2]));
            sum[4] = sum[4] + abs(abs(comp[i][1])-abs(temp[i][0]));
            sum[4] = sum[4] + abs(abs(comp[i][2])-abs(temp[i][1]));
        }
        
        // sum[5] holds the ZYX
        for (i=0;i<nat;i++)
        {
            sum[5] = sum[5] + abs(abs(comp[i][0])-abs(temp[i][2]));
            sum[5] = sum[5] + abs(abs(comp[i][1])-abs(temp[i][1]));
            sum[5] = sum[5] + abs(abs(comp[i][2])-abs(temp[i][0]));
        }
        
        
        // evalaute the sums and make a decision:
        // find out which sum is the smallest
        //    for (i=0;i<6;i++)
        //    {cout << sum[i] << " ";}
        //    cout << "\n";
        
        ind=0;
        min = sum[0];
        for (k=1;k<6;k++)
        {
            if (sum[k] < min)
            {
                ind = k;
                min = sum[k];
            }
        }
        
        ind = ind + 1;
        //  cout << ind << "\n";
        switch (ind)
        {
            case 1: // XYZ
                for (int a=0; a<nat; a++)
                {
                    temp[a][0] = temp[a][0];
                    temp[a][1] = temp[a][1];
                    temp[a][2] = temp[a][2];
                }
                
                break;
                
            case 2: // XZY
                for (int a=0; a<nat; a++)
                {
                    temp[a][0] = temp[a][0];
                    aa = temp[a][1];
                    temp[a][1] = temp[a][2];
                    temp[a][2] = aa;
                }
                break;
                
            case 3: // YXZ
                for (int a=0; a<nat; a++)
                {
                    temp[a][2] = temp[a][2];
                    aa = temp[a][0];
                    temp[a][0] = temp[a][1];
                    temp[a][1] = aa;
                }
                break;
                
            case 4: // YZX
                for (int a=0; a<nat; a++)
                {
                    aa = temp[a][0];
                    bb = temp[a][1];
                    cc = temp[a][2];
                    temp[a][0] = bb;
                    temp[a][1] = cc;
                    temp[a][2] = aa;
                }
                break;
                
            case 5: // ZXY
                for (int a=0; a<nat; a++)
                {
                    aa = temp[a][0];
                    bb = temp[a][1];
                    cc = temp[a][2];
                    temp[a][0] = cc;
                    temp[a][1] = aa;
                    temp[a][2] = bb;
                }
                break;
                
            case 6: // ZYX
                for (int a=0; a<nat; a++)
                {
                    aa = temp[a][0];
                    bb = temp[a][1];
                    cc = temp[a][2];
                    temp[a][0] = cc;
                    temp[a][1] = bb;
                    temp[a][2] = aa;
                }
                break;
                
                
                
            default: cout << "index has not been determined - error in structures\n";
                break;
        }
        
        
        // write the data back into the array:
        
        // return the geometry in the new form to the main geometry matrix
        for (k=0; k<nat; k++)
        {
            l=k*3;
            pred[0][0+l] = temp[k][0];
            pred[0][1+l] = temp[k][1];
            pred[0][2+l] = temp[k][2];
            
        }
        
        
        
        // end of function
        
        
        
        
    }
    
    
    
    
    
    // function ends here
}



// findCut finds the cut off parameter for a particular grid
void findCut(MatDoub &xmat, Int &numpoints, Int &ndim, Doub &cut, VecDoub &delta)
{
	Int i;
	Doub maxX, maxY, maxZ;
	Doub minX, minY, minZ;
	Doub xcord,ycord,zcord;
    Doub grid,grid2D;
    
    grid = numpoints/50;


    if (ndim == 3)
    {
        maxX=xmat[0][0];
        minX=xmat[0][0];
        maxY=xmat[0][1];
        minY=xmat[0][1];
        maxZ=xmat[0][2];
        minZ=xmat[0][2];
        
        for (i=1; i<numpoints; i++) 
        {
            //check x
            if (xmat[i][0] > maxX ){maxX = xmat[i][0];}
            if (xmat[i][0] < minX ){minX = xmat[i][0];}
            
            //check y
            if (xmat[i][1] > maxY ){maxY = xmat[i][1];}
            if (xmat[i][1] < minY ){minY = xmat[i][1];}			
            
            //check z
            if (xmat[i][2] > maxZ ){maxZ = xmat[i][2];}
            if (xmat[i][2] < minZ ){minZ = xmat[i][2];}			
            
            // for loop end
        }   
         
        xcord = (maxX - minX)/grid;
        ycord = (maxY - minY)/grid;
        zcord = (maxZ - minZ)/grid;
        
        delta[0] = (maxX - minX)/50;        
        delta[1] = (maxY - minY)/50;
        delta[2] = (maxZ - minZ)/50;
        
        cut = sqrt(xcord * xcord + ycord * ycord + zcord * zcord);
        
    }
    
    if (ndim == 2)
    {
        maxX=xmat[0][0];
        minX=xmat[0][0];
        maxY=xmat[0][1];
        minY=xmat[0][1];
        
        for (i=1; i<numpoints; i++) 
        {
            //check x
            if (xmat[i][0] > maxX ){maxX = xmat[i][0];}
            if (xmat[i][0] < minX ){minX = xmat[i][0];}
            
            //check y
            if (xmat[i][1] > maxY ){maxY = xmat[i][1];}
            if (xmat[i][1] < minY ){minY = xmat[i][1];}			
            
            // for loop end
        }
        
        grid2D = 1.0;
        xcord = (maxX - minX)/grid2D;
        ycord = (maxY - minY)/grid2D;
        
        delta[0] = (maxX - minX)/50;        
        delta[1] = (maxY - minY)/50;
        
        cut = sqrt(xcord* xcord + ycord * ycord);
        // cout << "grid is " << grid << "\n";
    }    
    
    if (ndim == 1)
    {
        maxX=xmat[0][0];
        minX=xmat[0][0];
        
        for (i=1; i<numpoints; i++) 
        {
            //check x
            if (xmat[i][0] > maxX ){maxX = xmat[i][0];}
            if (xmat[i][0] < minX ){minX = xmat[i][0];}
            
            // for loop end
        }
        
       cut = (maxX - minX)/3;
        
    }  
	// function end
}

void transpose(MatDoub &inmat, MatDoub &outmat, Int &rows, Int &cols)
{
    int i,j;
    for(i=0;i<rows;i++)
    {
        for(j=0;j<cols;j++)
        {
            outmat[j][i]=inmat[i][j];
        }
    }
    
}

void matmult(MatDoub &a, MatDoub &b, MatDoub &c, Int &arow, Int &acol, Int &brow, Int &bcol)
{
    int i,j,k;
    Doub sum;
    // cmatrix will be arow x bcol
    
    for (i=0; i<arow; i++) 
    {
        for (j=0; j<bcol; j++) 
        {   
            sum = 0;
            for (k=0; k<acol; k++) 
            {
                sum += a[i][k]*b[k][j];
            }
            c[i][j] = sum;
        }
    }
    
    
    
}



void printStuff(MatDoub &matTpPrint, Int &rows, Int &columns)
{
    int i,j;
    
    
    for (i=0; i<rows; i++) 
    {
        for (j=0; j<columns; j++) 
        {
            cout << setw(15) << matTpPrint[i][j] << " ";
        }
        cout << "\n";
    } 
    
}

void printGmat(MatDoub &gmatrix, MatDoub &internals, Int &dim, Int &numstructs)
{
    int i,j,k;
    
    
    // print the 1D matrix
    if (dim == 1)
    {
        for (i=0; i<numstructs; i++)
        {
            cout << internals[i][0] << " ";
            cout << setw(15) << gmatrix[i][0];
            cout << "\n";
            
        }
    }
    
    
    // print the 2D matrix
    if (dim == 2)
    {
        for (i=0; i<numstructs; i++) 
        {
            cout << internals[i][0] << " " << internals[i][1] << " ";
            for (j=0; j<dim*dim; j++) 
            {
                cout << setw(15) << gmatrix[i][j] << " ";
            }
            cout << "\n";
            
        }
    }
    
    // print the   Gmatrix in 3D
    if (dim == 3)
    {
        for (i=0; i<numstructs; i++) 
        {
            cout << internals[i][0] << " " << internals[i][1] << " " << internals[i][2] << " ";
            for (j=0; j<dim*dim; j++) 
            {
                cout << setw(25) << gmatrix[i][j];
             //   cout  << " " << j+1 << " ";
            }
            cout << "\n";
            
        }
    }
}

// analytical gmatrix representation:
void gmatAnalytical(MatDoub &gmat, MatDoub &internals, MatDoub &xg, VecDoub &mas, Int &structures, Int nat)
{
    int i,j,k,l,h1=0,h2=0,h3=0;
    double angle_1, angle_2, angle_3;
    MatDoub temp(nat,3);
    MatDoub comp(nat,3);
    double aa,bb,cc;
    double int_1,int_2,int_3;
    double min;
    int ind;
    double temp_m_1, temp_m_2, temp_m_3;
    
    for (j=0; j<nat; j++) // create a temporary geometry storage for the given structure
    {
        l=3*j;
        temp[j][0] = xg[0][0+l];
        temp[j][1] = xg[0][1+l];
        temp[j][2] = xg[0][2+l];
    }
    
    int_1 = internals[0][0];
    int_2 = internals[0][1];
    int_3 = internals[0][2];
   
    // if angle_1 -- central carbon is mass[1] and no modification is needed
    angle_1=acos(
    (
    (temp[1][0]-temp[0][0])*(temp[1][0]-temp[2][0]) +
    (temp[1][1]-temp[0][1])*(temp[1][1]-temp[2][1]) +
    (temp[1][2]-temp[0][2])*(temp[1][2]-temp[2][2])
    ) /
    (
     sqrt(pow((temp[1][0]-temp[0][0]),2)+pow((temp[1][1]-temp[0][1]),2) + pow((temp[1][2]-temp[0][2]),2))*
     sqrt(pow((temp[2][0]-temp[1][0]),2)+pow((temp[2][1]-temp[1][1]),2) + pow((temp[2][2]-temp[1][2]),2))
     )
    );
    
    // if angle_2 -- central carbon is mass[2];
    angle_2=acos(
                 (
                  (temp[2][0]-temp[1][0])*(temp[2][0]-temp[0][0]) +
                  (temp[2][1]-temp[1][1])*(temp[2][1]-temp[0][1]) +
                  (temp[2][2]-temp[1][2])*(temp[2][2]-temp[0][2])
                  ) /
                 (
                  sqrt(pow(temp[1][0]-temp[2][0],2)+pow(temp[1][1]-temp[2][1],2) + pow(temp[1][2]-temp[2][2],2))*
                  sqrt(pow(temp[2][0]-temp[0][0],2)+pow(temp[2][1]-temp[0][1],2) + pow(temp[2][2]-temp[0][2],2))
                  )
                 );
    
    // if angle_3 -- central carbon is mass[0];
    angle_3=acos(
                 (
                  (temp[0][0]-temp[1][0])*(temp[0][0]-temp[2][0]) +
                  (temp[0][1]-temp[1][1])*(temp[0][1]-temp[2][1]) +
                  (temp[0][2]-temp[1][2])*(temp[0][2]-temp[2][2])
                  ) /
                 (
                  sqrt(pow(temp[1][0]-temp[0][0],2)+pow(temp[1][1]-temp[0][1],2) + pow(temp[1][2]-temp[0][2],2))*
                  sqrt(pow(temp[2][0]-temp[0][0],2)+pow(temp[2][1]-temp[0][1],2) + pow(temp[2][2]-temp[0][2],2))
                  )
                 );
    
   // cout << angle_1 << " " << angle_2 << " " << angle_3 << "\n";
    
    // Adjusting masses:
    
    if (abs(angle_2 - int_1) < 0.005 || abs(angle_2 - int_2) < 0.005 || abs(angle_2 - int_3) < 0.005 )
    {
        // mass[2] is central
        temp_m_1=mas[0];
        temp_m_2=mas[1];
        temp_m_3=mas[2];
        
        mas[1] = temp_m_3;
        mas[2] = temp_m_2;
    }
    
    if (abs(angle_3 - int_1) < 0.005 || abs(angle_3 - int_2) < 0.005 || abs(angle_3 - int_3) < 0.005 )
    {
        // mass[0] is central
        temp_m_1=mas[0];
        temp_m_2=mas[1];
        temp_m_3=mas[2];
        
        mas[0] = temp_m_2;
        mas[1] = temp_m_1;
    }
    
    
    // r12
    aa = sqrt( pow((temp[0][0] - temp[1][0]),2) + pow((temp[0][1] - temp[1][1]),2) + pow((temp[0][2] - temp[1][2]),2) );
    // r13
    bb = sqrt( pow((temp[0][0] - temp[2][0]),2) + pow((temp[0][1] - temp[2][1]),2) + pow((temp[0][2] - temp[2][2]),2) );
    // r32
    cc = sqrt( pow((temp[2][0] - temp[1][0]),2) + pow((temp[2][1] - temp[1][1]),2) + pow((temp[2][2] - temp[1][2]),2) );
    

    
    if (abs(int_1 - aa) < 0.001 || abs(int_1 - bb) < 0.001 || abs(int_1 - cc) < 0.001)
    {
        h1=1;
    }
    if (abs(int_2 - aa) < 0.001 || abs(int_2 - bb) < 0.001 || abs(int_2 - cc) < 0.001)
    {
        h2=1;
    }
        
    
    
    
    if (h1 == 1 && h2 == 1)
    {
        // angle is the thrird coordinate
        for (i=0; i<structures; i++)
        {
            gmat[i][0] = 1/mas[0] + 1/mas[1]; // 1,1
            gmat[i][1] = 1/mas[1]*cos(internals[i][2]); // 1,2
            gmat[i][2] = -1/mas[1]*1/internals[i][1]*sin(internals[i][2]); // 1,3
            gmat[i][3] = gmat[i][1]; // 2,1
            gmat[i][4] = 1/mas[1] + 1/mas[2]; // 2,2
            gmat[i][5] = -1/mas[1]*1/internals[i][0]*sin(internals[i][2]); // 2,3
            gmat[i][6] = gmat[i][2]; // 3,1
            gmat[i][7] = gmat[i][5]; // 3,2
            gmat[i][8] = 1/mas[0]*1/(internals[i][0]*internals[i][0]) +
                        1/mas[2]*1/(internals[i][1]*internals[i][1]) +
                    1/mas[1]*(1/(internals[i][0]*internals[i][0]) + 
                              1/(internals[i][1]*internals[i][1]) +
                              -2*cos(internals[i][2])/(internals[i][0]*internals[i][1]));
        
        }
    }
    
    else if (h1 == 1 && h2 ==0)
    {
        // angle is the second coordinate
        for (i=0; i<structures; i++)
        {
            gmat[i][0] = 1/mas[0] + 1/mas[1]; // 1,1
            gmat[i][1] =  -1/mas[1]*1/internals[i][2]*sin(internals[i][1]); // 1,2
            gmat[i][2] = 1/mas[1]*cos(internals[i][1]); // 1,3
   
            gmat[i][3] = gmat[i][1]; // 2,1
            gmat[i][4] = 1/mas[0]*1/(internals[i][0]*internals[i][0]) +
            1/mas[2]*1/(internals[i][2]*internals[i][2]) +
            1/mas[1]*(1/(internals[i][0]*internals[i][0]) +
                      1/(internals[i][2]*internals[i][2]) +
                      -2*cos(internals[i][1])/(internals[i][0]*internals[i][2])); // 2,2
            gmat[i][5] = -1/mas[1]*1/internals[i][0]*sin(internals[i][1]); // 2,3
            
            gmat[i][6] = gmat[i][2]; // 3,1
            gmat[i][7] = gmat[i][5]; // 3,2
            gmat[i][8] = 1/mas[1] + 1/mas[2]; // 3,3
        }
        
    }
    else if (h1 == 0 && h2 ==1)
    {
        // angle is the first coordinate
        for (i=0; i<structures; i++)
        {
            gmat[i][0] = 1/mas[0]*1/(internals[i][1]*internals[i][1]) +
            1/mas[2]*1/(internals[i][2]*internals[i][2]) +
            1/mas[1]*(1/(internals[i][1]*internals[i][1]) +
                      1/(internals[i][2]*internals[i][2]) +
                      -2*cos(internals[i][0])/(internals[i][1]*internals[i][2])); // 1,1
            gmat[i][1] = -1/mas[1]*1/internals[i][2]*sin(internals[i][0]); // 1,2
            gmat[i][2] = -1/mas[1]*1/internals[i][1]*sin(internals[i][0]); // 1,3
           
            
            gmat[i][3] = gmat[i][1]; // 2,1
            gmat[i][4] = 1/mas[0] + 1/mas[1]; // 2,2
            gmat[i][5] = 1/mas[1]*cos(internals[i][0]); // 2,3
            
            
            gmat[i][6] = gmat[i][2]; // 3,1
            gmat[i][7] = gmat[i][5]; // 3,2
            gmat[i][8] = 1/mas[1] + 1/mas[2]; // 3,3
            
        }
        
    }
    else
    {
        // default if none of the above work:
        for (i=0; i<structures; i++)
        {
            gmat[i][0] = 1/mas[0] + 1/mas[1]; // 1,1
            gmat[i][1] = 1/mas[1]*cos(internals[i][2]); // 1,2
            gmat[i][2] = -1/mas[1]*1/internals[i][1]*sin(internals[i][2]); // 1,3
            gmat[i][3] = gmat[i][1]; // 2,1
            gmat[i][4] = 1/mas[1] + 1/mas[2]; // 2,2
            gmat[i][5] = -1/mas[1]*1/internals[i][0]*sin(internals[i][2]); // 2,3
            gmat[i][6] = gmat[i][2]; // 3,1
            gmat[i][7] = gmat[i][5]; // 3,2
            gmat[i][8] = 1/mas[0]*1/(internals[i][0]*internals[i][0]) +
            1/mas[2]*1/(internals[i][1]*internals[i][1]) +
            1/mas[1]*(1/(internals[i][0]*internals[i][0]) +
                      1/(internals[i][1]*internals[i][1]) +
                      -2*cos(internals[i][2])/(internals[i][0]*internals[i][1]));
            
        }
        
        
    }
    
    
    // function end
}

void playdata(MatDoub &xg, Int &nat, Int &nstr, MatDoub &ref) //, MatDoub &Internals)
{
    // xg size is (nstr , 3*nat)
    int i,j,k,l,ind;
    MatDoub m1(nat,3),m2(nat,3),m3(nat,3); // these will hold the rotated structures
    MatDoub m90x(nat,3),m90y(nat,3), m90z(nat,3);
    MatDoub temp(nat,3),comp(nat,3);
    Doub min,min90;
    VecDoub sum(8);
    
    MatDoub predStruct(1,nat*3); // used as temp storage to send stuff to playdata
    MatDoub poStruct(1,nat*3); // used as temp storage to send stuff to playdata

    for (j=0; j<nat; j++) // create a temporary geometry storage for the reference
    {
        l=3*j;    
        comp[j][0] = ref[0][0+l];
        comp[j][1] = ref[0][1+l];
        comp[j][2] = ref[0][2+l];
    }
 
    // cycle through all the structures and check the signs
    for (i=1;i<nstr; i++)
    {
        
     /*   for (j=0; j<nat; j++) // create a temporary geometry storage for the reference
        {
            l=3*j;
            comp[j][0] = xg[i-1][0+l];
            comp[j][1] = xg[i-1][1+l];
            comp[j][2] = xg[i-1][2+l];
        }
      */
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
        /*
        for (int o=0; o<nat; o++)
        {
            for (int oo=0;oo<3;oo++)
            {
                cout << comp[o][oo] << " ";
            }
        }
        cout << "\n";
        */
        /*
        for (int o=0; o<nat; o++)
        {
            for (int oo=0;oo<3;oo++)
            {
            cout << temp[o][oo] << " ";
            }
        }
        cout << "\n";
         */
        // calculate sum corresponding to the geometry multiplied by "-1"
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[1] = sum[1] + abs(-temp[a][c] - comp[a][c]);
            }
        }
        /*
        for (int o=0; o<nat; o++)
        {
            for (int oo=0;oo<3;oo++)
            {
                cout << -temp[o][oo] << " ";
            }
        }
        cout << "\n";
        */
        // === Z AXIS =========
        // rotate the geometry of i-th structure along z-axis 180 deg
        
        for (k=0; k<nat; k++ ) 
        {
            m1[k][0] = -temp[k][0]; 
            m1[k][1] = -temp[k][1];
            m1[k][2] =  temp[k][2];
        } // rotation along z-axis of i-th structure finished
       /*
        for (int o=0; o<nat; o++)
        {
            for (int oo=0;oo<3;oo++)
            {
                cout << m1[o][oo] << " ";
            }
        }
        cout << "\n";
        */
        // calculate sum corresponding to the rotated geometry 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[2] = sum[2] + abs(m1[a][c] - comp[a][c]);
            }
        }    
        /*
        for (int o=0; o<nat; o++)
        {
            for (int oo=0;oo<3;oo++)
            {
                cout << -m1[o][oo] << " ";
            }
        }
        cout << "\n";
         */
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
        /*
        for (int o=0; o<nat; o++)
        {
            for (int oo=0;oo<3;oo++)
            {
                cout << m2[o][oo] << " ";
            }
        }
        cout << "\n";
        */
        // calculate sum corresponding to the rotated geometry multiplied by "-1" 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[5] = sum[5] + abs(-m2[a][c] - comp[a][c]);
            }
        }
        /*
        for (int o=0; o<nat; o++)
        {
            for (int oo=0;oo<3;oo++)
            {
                cout << -m2[o][oo] << " ";
            }
        }
        cout << "\n";
        */
        // ===== Y AXIS ======
        // rotate the geometry of i-th structure along y-axis 180 deg
        
        for (k=0; k<nat; k++ ) 
        {
            m3[k][0] = -temp[k][0]; // + 2*temp[k][0]; 
            m3[k][1] =  temp[k][1];
            m3[k][2] = -temp[k][2];// + 2*temp[k][2];
        } // rotation along z-axis of i-th structure finished
        /*
        for (int o=0; o<nat; o++)
        {
            for (int oo=0;oo<3;oo++)
            {
                cout << m3[o][oo] << " ";
            }
        }
        cout << "\n";
        */
         // calculate sum corresponding to the rotated geometry
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[6] = sum[6] + abs(m3[a][c] - comp[a][c]);
            }
        }
        /*
        for (int o=0; o<nat; o++)
        {
            for (int oo=0;oo<3;oo++)
            {
                cout << -m3[o][oo] << " ";
            }
        }
        cout << "\n";
        */
        // calculate sum corresponding to the rotated geometry multiplied by "-1" 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[7] = sum[7] + abs(-m3[a][c] - comp[a][c]);
            }
        } 
        
        
        
        /*
         for (i=0; i<natom; i++ ) 
         {
         rotstruct[i][0] =   -xgeom[i][1];
         rotstruct[i][1] =    xgeom[i][0];
         rotstruct[i][2] =    xgeom[i][2];
         }
        */
        // === Z AXIS =========
        // rotate the geometry of i-th structure along z-axis 180 deg
     /*   
        for (k=0; k<nat; k++ ) 
        {
            m90z[k][0] = -temp[k][1]; 
            m90z[k][1] =  temp[k][0];
            m90z[k][2] =  temp[k][2];
        } // rotation along z-axis of i-th structure finished
        
        // calculate sum corresponding to the rotated geometry 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[8] = sum[8] + abs(m90z[a][c] - comp[a][c]);
            }
        }    
        
        // calculate sum corresponding to the rotated geometry multiplied by "-1" 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[9] = sum[9] + abs(-m90z[a][c] - comp[a][c]);
            }
        }        
        
        // ===== X AXIS ======
        // rotate the geometry of i-th structure along x-axis 180 deg
        
        for (k=0; k<nat; k++ ) 
        {
            m90x[k][0] =  temp[k][0]; 
            m90x[k][1] = -temp[k][2];
            m90x[k][2] =  temp[k][1];// + 2*temp[k][2];
        } // rotation along z-axis of i-th structure finished
        
        // calculate sum corresponding to the rotated geometry 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[10] = sum[10] + abs(m90x[a][c] - comp[a][c]);
            }
        }    
        
        // calculate sum corresponding to the rotated geometry multiplied by "-1" 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[11] = sum[11] + abs(-m90x[a][c] - comp[a][c]);
            }
        }  
        
        // ===== Y AXIS ======
        // rotate the geometry of i-th structure along y-axis 180 deg
        
        for (k=0; k<nat; k++ ) 
        {
            m3[k][0] =  temp[k][2]; // + 2*temp[k][0]; 
            m3[k][1] =  temp[k][1];
            m3[k][2] = -temp[k][0];// + 2*temp[k][2];
        } // rotation along z-axis of i-th structure finished
        
        // calculate sum corresponding to the rotated geometry 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[12] = sum[12] + abs(m90y[a][c] - comp[a][c]);
            }
        }    
        
        // calculate sum corresponding to the rotated geometry multiplied by "-1" 
        for (int a=0; a<nat; a++) // loop through atoms
        {
            for (int c=0; c<3; c++) // loop through coordinates
            {
                sum[13] = sum[13] + abs(-m90y[a][c] - comp[a][c]);
            }
        } 
    */    
        // ===== 90 degree finished
        
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
     /*
        for (int h=0;h<8;h++)
        {
            cout << sum[h] << " ";
        }
       */
    //    cout << "\n" << min << " " << ind+1 << "\n";
        ind = ind + 1;
    //   cout << i+1 << " " << Internals[i][0] << " " << Internals[i][1] << " " << Internals[i][2] << " "; 
        switch (ind) 
        {
            case 1:
	//	cout << "1\n";
                for (int a=0; a<nat; a++)
                {
                    for (int c=0;c<3;c++)
                    {
                            temp[a][c] = -temp[a][c];
                      //  temp[a][c] = temp[a][c];
                    }
                }

                break;
                
            case 2:
	//	cout << "2\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                    //    temp[a][c] = -temp[a][c];
                        temp[a][c] = temp[a][c];
                    }
                }
                break;
                
            case 3:
	//	cout << "3\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                     //   temp[a][c] = m1[a][c];
                        temp[a][c] = -m1[a][c];
                    }
                }
                break;

            case 4:
	//	cout << "4\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                      //  temp[a][c] = -m1[a][c];
                       temp[a][c] = m1[a][c];
                    }
                }
                break;
                
            case 5:
	//	cout << "5\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                    //    temp[a][c] = m2[a][c];
                        temp[a][c] = -m2[a][c];
                    }
                }
                break;    
                
            case 6:
	//	cout << "6\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                     //  temp[a][c] = -m2[a][c];
                         temp[a][c] = m2[a][c];
                    }
                }
                break;

            case 7:
	//	cout << "7\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                     //   temp[a][c] = m3[a][c];
                        temp[a][c] = -m3[a][c];
                    }
                }
                break;
            
            case 8:
	//	cout << "8\n";
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                       // temp[a][c] = -m3[a][c];
                        temp[a][c] = m3[a][c];
                    }
                }
                break;

         /*   case 9:
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = m90z[a][c];     
                    }
                }
                break;
                
            case 10:
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = -m90z[a][c];     
                    }
                }
                break;
                
            case 11:
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = m90x[a][c];     
                    }
                }
                break;
                
            case 12:
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = -m90x[a][c];     
                    }
                }
                break;
                
            case 13:
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = m90y[a][c];     
                    }
                }
                break;
                
            case 14:
                for (int a=0; a<nat; a++) 
                {
                    for (int c=0;c<3;c++)
                    {
                        temp[a][c] = -m90y[a][c];     
                    }
                }
                break;
                
            */    
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
       // cout << "\n";
        
        // ======= ADDED ======= Oct 2 2013 ==========
        Int p = 0;
        for (k=0; k<3*nat; k++)
        {
            predStruct[0][k] = xg[i-1][k];
            poStruct[0][k] = xg[i][k];
        }
        //  exchange(poStruct, numberofAtoms, predStruct, p);
        //  playdata(poStruct, numberofAtoms, Play_data_const, predStruct);
        check_alignment(poStruct, nat, predStruct);
        exchange(poStruct, nat, predStruct, p);
        // write the structures back:
        for (k=0; k<3*nat; k++)
        {
            xg[i][k] = poStruct[0][k];
        }
        
        
        
        // ==================== Oct 2 2013 ============
        
        

    } // and of structure loop
    
    // end of function
}



