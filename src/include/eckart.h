//
//  eckart.h
//  GMAT
//
//  Created by PETER ZAJAC on 10/4/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

void translation(int &numberofStructures, int &numberofAtoms, MatDoub &xgeoms, VecDoub &mass)
{
    int i,j,k,l,m,p;
    MatDoub temp(numberofAtoms,3);
    VecDoub origin(3); // center of mass
    double totalmass;

    
    for (i=0;i<numberofStructures;i++)
    {
        totalmass = 0.0;
        for (j=0; j<numberofAtoms; j++) // create a temporary geometry storage for the given structure
        {
            l=3*j;    
            temp[j][0] = xgeoms[i][0+l];
            temp[j][1] = xgeoms[i][1+l];
            temp[j][2] = xgeoms[i][2+l];
        }
        
        // Calculate the center of Mass:
        for (k=0; k<numberofAtoms; k++) 
        {
            for (m=0; m<3; m++) 
            {
                origin[m] = origin[m] + mass[k]*temp[k][m];
            }
            totalmass = totalmass + mass[k];
        }
        
        for (p=0; p<3; p++) 
        {
            origin[p] = origin[p]/totalmass;
        }
        
        // reverse the process and write the new geometry into the original array
        for (j=0; j<numberofAtoms; j++)
        {
            l=3*j;    
            xgeoms[i][0+l] = temp[j][0] - origin[0];
            xgeoms[i][1+l] = temp[j][1] - origin[1];
            xgeoms[i][2+l] = temp[j][2] - origin[2];
        }
        
        
        
        
    } // end the loop for Structures
}

void rotation(MatDoub &refx, Int &numberofStructures, Int &numberofAtoms, MatDoub &xgeoms, VecDoub &mass)
{
    int i,j,l,k;
    
    MatDoub temp(numberofAtoms,3);
    MatDoub tempxref(numberofAtoms,3);
    VecDoub origin(3); // center of mass
    Doub tol=1e-9, cs,sn,sp,cc;
    Doub ValueAtMin, total=0.0, total2=0.0 ,minimum,delta=0.2;
    Doub ValueAtMin2, ValueAtMin3;
    VecDoub point(4);
    VecDoub pmin(4);
    VecDoub vp(3);
    VecDoub pmin2(4);
    
    pmin[0] = 1.0;
    pmin[1] = 0.0;
    pmin[2] = 0.0;
    pmin[3] = 0.1;
    
    for (j=0; j<numberofAtoms; j++) // create a temporary geometry storage for the reference structure
    {
        l=3*j;    
        tempxref[j][0] = refx[0][0+l];
        tempxref[j][1] = refx[0][1+l];
        tempxref[j][2] = refx[0][2+l];
    }  
    
    // know loop through all the structures and correct them for rotation
    for (i=0;i<numberofStructures;i++)
    {
        
        for (j=0; j<numberofAtoms; j++) // create a temporary geometry storage for the given structure
        {
            l=3*j;    
            temp[j][0] = xgeoms[i][0+l];
            temp[j][1] = xgeoms[i][1+l];
            temp[j][2] = xgeoms[i][2+l];
        }

        // here I need a function to do the rotation - minimization
        
        functionRot ffn(tempxref,temp,numberofAtoms,mass);
        
        minimum = 1.0;
        Int iter = 100;
        
        // cycle through the first direction x
        for (int yy=0; yy<4; yy++) 
        {
            pmin2[0] = 0.0;
            pmin2[1] = 0.0;
            pmin2[2] = 1.0;
            pmin2[3] = 0.1 + 0.1*yy;
            dfpmin(pmin2, tol, iter, ValueAtMin, ffn);
            if (ValueAtMin < minimum ) 
            {
                minimum = ValueAtMin;
                for (int hh=0; hh<4; hh++) 
                {
                    point[hh] = pmin2[hh];
                }
            }
        
        } // cycling finished
       // cout << minimum << "  ";
        
        // cycle through the second diretion
        for (int yy=0; yy<5; yy++) 
        {
            pmin2[0] = 0.0;
            pmin2[1] = 1.0;
            pmin2[2] = 0.0;
            pmin2[3] = 0.1 + 0.1*yy;
            dfpmin(pmin2, tol, iter, ValueAtMin2, ffn);
            
            if (ValueAtMin2 < minimum ) 
            {
                minimum = ValueAtMin2;
                for (int hh=0; hh<4; hh++) 
                {
                    point[hh] = pmin2[hh];
                }
            }
            
        }  // finished cycling through the second direction 
         // cout << minimum << "  ";   
         
        for (int yy=0; yy<5; yy++) 
        {
        
            pmin2[0] = 1.0;
            pmin2[1] = 0.0;
            pmin2[2] = 0.0;
            pmin2[3] = 0.1 + 0.1*yy;
            dfpmin(pmin2, tol, iter, ValueAtMin3, ffn);
            if (ValueAtMin3 < minimum ) 
            {
                minimum = ValueAtMin3;
                for (int hh=0; hh<4; hh++) 
                {
                    point[hh] = pmin2[hh];
                }
            }
        
        }
      
       
       // cout << minimum << "      " << point[0] << " " << point[1] << " " << point[2] << " " << point[3] << " \n" ;
            
        
         normalize(point);    
       // cout << minimum << "  ";
       // cout << point[0] << " " << point[1] << " " << point[2] << " " << point[3] << "   \n "; // << ffn(point) << "\n";
       // cout << "\n";
          
        
        
        // now rotate the structure according to the point
        cs = cos(point[3]);
        sn = sin(point[3]);
        for (k=0; k<numberofAtoms; k++ ) 
        {
            sp = 0.0;
            for (j=0; j<3; j++) 
            {
                sp = sp + temp[k][j]*point[j];
            }
            cc = (1.0-cs)*sp;
            vp[0] = point[1]*temp[k][2] - point[2]*temp[k][1];
            vp[1] = point[2]*temp[k][0] - point[0]*temp[k][2];
            vp[2] = point[0]*temp[k][1] - point[1]*temp[k][0];
            temp[k][0] = cs*temp[k][0] + cc*point[0] + sn*vp[0];
            temp[k][1] = cs*temp[k][1] + cc*point[1] + sn*vp[1];
            temp[k][2] = cs*temp[k][2] + cc*point[2] + sn*vp[2];
        } // rotation of i-th structure finished
        
        // check inversion
      //  checkinv(tempxref, temp, numberofAtoms);
        
        
        for (k=0; k<numberofAtoms; k++) 
        {
            l=k*3;
            xgeoms[i][0+l] = temp[k][0];
            xgeoms[i][1+l] = temp[k][1];
            xgeoms[i][2+l] = temp[k][2];

        }
        
    } // cycle through the structures consinues take the i+1 structure
  
   // cout << " two totals: " << total2 << " \n "; // << total2 << "\n"; 
}

void rot_princ_axis(Int &numberofStructures, Int &numberofAtoms, MatDoub &xgeoms, VecDoub &mass)
{
    int i,j,k,l;
    MatDoub temp(numberofAtoms,3);
    MatDoub A(3,3);  // this matrix is used as MI tensor
    MatDoub eig_vec(3,3); // this will hold the eigenvectors
    Doub mxx,myy,mzz,mxy,mxz,myz;
    Doub xx,yy,zz;
    
    for (i=0;i<numberofStructures;i++)
    {
        for (j=0; j<numberofAtoms; j++) // create a temporary geometry storage for the given structure
        {
            l=3*j;    
            temp[j][0] = xgeoms[i][0+l];
            temp[j][1] = xgeoms[i][1+l];
            temp[j][2] = xgeoms[i][2+l];
        }
        
        // set the moments to zero:
        mxx=0.0; mxy=0.0; myy=0.0;
        mxz=0.0; myz=0.0; mzz=0.0;

        // Calculate the Moment of Inertia tensor:
        for (k=0; k<numberofAtoms; k++) 
        {
            mxx = mxx + mass[k]*(temp[k][1]*temp[k][1] + temp[k][2]*temp[k][2]);
            myy = myy + mass[k]*(temp[k][0]*temp[k][0] + temp[k][2]*temp[k][2]);
            mzz = mzz + mass[k]*(temp[k][0]*temp[k][0] + temp[k][1]*temp[k][1]);
            mxy = mxy - mass[k]*temp[k][0]*temp[k][1];
            mxz = mxz - mass[k]*temp[k][0]*temp[k][2];
            myz = myz - mass[k]*temp[k][1]*temp[k][2];
        }
        
        // set up the moment of inertia tensor:
        A[0][0] = mxx; A[1][1] = myy; A[2][2] = mzz;
        A[0][1] = mxy; A[1][0] = mxy; A[0][2] = mxz;
        A[2][0] = mxz; A[1][2] = myz; A[2][1] = myz;
        
        Jacobi jac(A); // constructor for the jacobi transformation
        eig_vec  = jac.v; // assign the eigenvectors to eig_vec
         
        mxy = 0.0; mxz = 0.0; myz = 0.0;
        for (k=0;k<numberofAtoms; k++)
        {
            xx = temp[k][0]*eig_vec[0][0] + temp[k][1]*eig_vec[1][0] + temp[k][2]*eig_vec[2][0];
            yy = temp[k][0]*eig_vec[0][1] + temp[k][1]*eig_vec[1][1] + temp[k][2]*eig_vec[2][1];
            zz = temp[k][0]*eig_vec[0][2] + temp[k][1]*eig_vec[1][2] + temp[k][2]*eig_vec[2][2];
            mxy = mxy - mass[k]*xx*yy;
            mxz = mxz - mass[k]*xx*zz;
            myz = myz - mass[k]*yy*zz;
            temp[k][0] = xx;
            temp[k][1] = yy;
            temp[k][2] = zz;
            
        }
      
        // this is to double check if it works fine
        /*
        if ( (abs(mxy)>1e-4) || (abs(mxz)>1e-4) || (abs(myz)>1e-4))
        {
            cout << "Structure: " << i << " principal axis problem\n";
        }
        */
        // reverse the process and write the new geometry into the original array
        for (j=0; j<numberofAtoms; j++)
        {
            l=3*j;    
            xgeoms[i][0+l] = temp[j][0];
            xgeoms[i][1+l] = temp[j][1];
            xgeoms[i][2+l] = temp[j][2];
        }
        
    } // end the loop for Structures
    
    
    
}

