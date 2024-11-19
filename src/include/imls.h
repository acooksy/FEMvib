//
//  imls.h
//  ILMS
//
//  Created by PETER ZAJAC on 9/28/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//


double imls(MatDoub &xmat, VecDoub &yvec, VecDoub &intPoint, Int &npt, Int &ndim, Doub &cut)
{
    Doub result = 0; // here the final result is stored
    Doub distance = 0; // the distance between the interpolated point and the point on the grid
    Doub sfunc = 0, dc; // the scalling function for weight function
    Doub eps = 1.e-12; // small parameter for weighting function
    int i,j,k,p;
    
    // 1D case scenario => our xmat is a vector
    if (ndim == 1)
    {
	Int npts=10;
        Int rowb=6; 
	MatDoub pts(npts,ndim);
	VecDoub yvecnew(npts);
	howmany(xmat,yvec, intPoint, ndim, npt, npts, pts,yvecnew);
        
	MatDoub B(npts,rowb);
        // assemble the B matrix representing the function:
        //  a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5
        for (i=0; i<npts; i++) 
        {
            B[i][0] = 1.0;
            B[i][1] = pts[i][0];
            B[i][2] = pts[i][0]*pts[i][0];
            B[i][3] = pts[i][0]*pts[i][0]*pts[i][0];
            B[i][4] = pts[i][0]*pts[i][0]*pts[i][0]*pts[i][0];
            B[i][5] = pts[i][0]*pts[i][0]*pts[i][0]*pts[i][0]*pts[i][0];
            
        } // Matrix B assembled
        
        
        MatDoub W(npts,npts);
        // assemble the weight matrix
        for (i=0; i<npts; i++) 
        {
            distance = abs(intPoint[0] - pts[i][0]);
            dc = distance;
                sfunc = (1 - dc*dc*dc*dc*dc*dc);
            W[i][i] = sfunc * exp (-(distance*distance) ) / ( pow(distance,8) + eps);
            
        }
        // assembly of W matrix finished

        // solve the following system: B' * W * B * a' = B' * W * f'
            // B' * W first: B' is 1 x npt and W is npt npt
            MatDoub Btrans(rowb,npts);
            MatDoub A1(rowb,npts);
            MatDoub A(rowb,rowb);
        
            transpose(B,Btrans,npts,rowb);
            matmult(Btrans,W,A1,rowb,npts,npts,npts);
            matmult(A1,B,A,rowb,npts,npts,rowb);
        // so the left hand side is stored in A, now we need rhs:
        
            VecDoub y(rowb);
            MatDoub y1(rowb,npts);
        
            matmult(Btrans,W,y1,rowb,npts,npts,npts);
            // vector matrix multiplication
        for (k=0; k<rowb; k++) 
        {
            y[k] = 0;
            for (i=0; i<npts; i++) 
            {
                y[k] += y1[k][i] * yvecnew[i];
            }
        }    
    
        // know we need to solve the nonlinear equation A*a=y, where
        // small a represents the coefficients A is the LHS matrix and y is RHS
            SVD decom(A);
            // initialize the vector in which solution is stored:
            VecDoub a(rowb);
            // solve the system:
            decom.solve(y,a);
        
        // finally evaluate the function B with its coeficients
        result = 0;
        for (i=0; i<rowb; i++)
        {
            result+=a[i]*pow(intPoint[0],i);
        }

    }
    // =========    1D case is done ====================
    
    
    // 2D case scenario => our xmat has 2 columns
    if (ndim == 2) 
    {
        
            MatDoub B(npt,8);
            VecDoub weights(npt);
            Int rowb = 8;
            // assemble the B matrix representing the function:
            //  a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y
            for (i=0; i<npt; i++) 
            {
                B[i][0] = 1.0;
                B[i][1] = xmat[i][0];
                B[i][2] = xmat[i][1];
                B[i][3] = xmat[i][0]*xmat[i][0];
                B[i][4] = xmat[i][1]*xmat[i][1];
                B[i][5] = xmat[i][0]*xmat[i][1];
                B[i][6] = xmat[i][0]*xmat[i][0]*xmat[i][0];
                B[i][7] = xmat[i][1]*xmat[i][1]*xmat[i][1];
            } // Matrix B assembled
            
            
            // assemble the weight matrix
            for (i=0; i<npt; i++) 
            {
                distance = sqrt ( 
                                 (intPoint[0] - xmat[i][0])*(intPoint[0] - xmat[i][0]) + 
                                 (intPoint[1] - xmat[i][1])*(intPoint[1] - xmat[i][1])
                                );
                dc = distance; // dc is just a dummy variable to simplify equations
                if (dc > 1) 
                {
                    sfunc = 0;
                }
                else
                {
                    sfunc = (1 - dc*dc*dc*dc*dc*dc);
                }
                weights[i] = sfunc * exp (-(distance*distance) ) / ( pow(distance,8) + eps);
                
            }
            
            MatDoub W(npt,npt);
            for (i=0; i<npt; i++) 
            {
                W[i][i] = weights[i];
            }
            // assembly of W matrix finished
            
            // solve the following system: B' * W * B * a' = B' * W * f'
            // B' * W first: B' is 1 x npt and W is npt npt
            MatDoub Btrans(rowb,npt);
            MatDoub A1(rowb,npt);
            MatDoub A(rowb,rowb);
            
            transpose(B,Btrans,npt,rowb);
            matmult(Btrans,W,A1,rowb,npt,npt,npt);
            matmult(A1,B,A,rowb,npt,npt,rowb);
            // so the left hand side is stored in A, now we need rhs:
            
            VecDoub y(rowb);
            MatDoub y1(rowb,npt);
            
            matmult(Btrans,W,y1,rowb,npt,npt,npt);
            // vector matrix multiplication
            for (k=0; k<rowb; k++) 
            {
                y[k] = 0;
                for (i=0; i<npt; i++) 
                {
                    y[k] += y1[k][i] * yvec[i];
                }
            }    
            
            // know we need to solve the nonlinear equation A*a=y, where
            // small a represents the coefficients A is the LHS matrix and y is RHS
            SVD decom(A);
            // initialize the vector in which solution is stored:
            VecDoub a(rowb);
            // solve the system:
            decom.solve(y,a);
            
            // finally evaluate the function B with its coeficients
            result = 0;
            result = a[0] + a[1]*intPoint[0] + a[2]*intPoint[1] + a[3]*intPoint[0]*intPoint[0] +
                     a[4]*intPoint[1]*intPoint[1] + a[5]*intPoint[0]*intPoint[1] +
                     a[6]*intPoint[0]*intPoint[0]*intPoint[0] + a[7]*intPoint[1]*intPoint[1]*intPoint[1];
            
    
    }
    // =========    2D case is done ====================
    
    
     
    // 3D case scenario => our xmat has 3 columns 
    if (ndim == 3) 
    {
	Int npts=30;
        Int rowb=10; 
	MatDoub pts(npts,ndim);
	VecDoub yvecnew(npts);
	howmany(xmat,yvec, intPoint, ndim, npt, npts, pts,yvecnew);

        MatDoub B(npts,rowb);
        // assemble the B matrix representing the function:
	for (i=0; i<npts; i++) 
        {
            B[i][0] = 1.0; // a0
            B[i][1] = pts[i][0]; // a1*x
            B[i][2] = pts[i][1]; // a2*y
            B[i][3] = pts[i][2]; // a3*z
            B[i][4] = pts[i][0]*pts[i][0]; // a4*x^2
            B[i][5] = pts[i][1]*pts[i][1]; // a5*y^2
            B[i][6] = pts[i][2]*pts[i][2]; // a6*z^2
            B[i][7] = pts[i][0]*pts[i][1]; // a7*x*y
            B[i][8] = pts[i][0]*pts[i][2]; // a8*x*z
            B[i][9] = pts[i][1]*pts[i][2]; // a9*y*z

        } // Matrix B assembled
        
        // assemble the weight matrix
        MatDoub W(npts,npts);
        for (p=0; p<npts; p++) 
        {
            distance = sqrt ( 
                             (intPoint[0] - pts[p][0])*(intPoint[0] - pts[p][0]) + 
                             (intPoint[1] - pts[p][1])*(intPoint[1] - pts[p][1]) +
                             (intPoint[2] - pts[p][2])*(intPoint[2] - pts[p][2])
                             );

            dc = distance; // dc is just a dummy variable to simplify equations
            sfunc = (1 - dc*dc*dc*dc*dc*dc);
            W[p][p] = sfunc * exp (-(distance*distance) ) / ( pow(distance,8) + eps);
        }
        // assembly of W matrix finished
        
        // solve the following system: B' * W * B * a' = B' * W * f'
        // B' * W first: B' is 1 x npt and W is npt npt
        MatDoub Btrans(rowb,npts);
        MatDoub A1(rowb,npts);
        MatDoub A(rowb,rowb);
        
        transpose(B,Btrans,npts,rowb);
        matmult(Btrans,W,A1,rowb,npts,npts,npts);
        matmult(A1,B,A,rowb,npts,npts,rowb);
        // so the left hand side is stored in A, now we need rhs:
        
        VecDoub y(rowb);
        MatDoub y1(rowb,npts);
        
        matmult(Btrans,W,y1,rowb,npts,npts,npts);
        for (k=0; k<rowb; k++) 
        {
            y[k] = 0;
            for (i=0; i<npts; i++) 
            {
                y[k] += y1[k][i] * yvecnew[i];
            }
        }    
        
        // know we need to solve the nonlinear equation A*a=y, where
        // small a represents the coefficients A is the LHS matrix and y is RHS
        SVD decom(A);
        // initialize the vector in which solution is stored:
        VecDoub a(rowb);
        // solve the system:
        decom.solve(y,a);
        
        // finally evaluate the function B with its coeficients
        result = a[0] + a[1]*intPoint[0] + a[2]*intPoint[1] + a[3]*intPoint[2] +
                 a[4]*intPoint[0]*intPoint[0] + a[5]*intPoint[1]*intPoint[1] +
                 a[6]*intPoint[2]*intPoint[2] + a[7]*intPoint[0]*intPoint[1] +
                 a[8]*intPoint[0]*intPoint[2] + a[9]*intPoint[1]*intPoint[2];
        
    
    }
    // =========    3D case is done ====================
    
    
    
    
    
    return result;
}

