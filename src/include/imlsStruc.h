//
//  imlsStruc.h
//  GMAT
//
//  Created by PETER ZAJAC on 10/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

struct imlsStruc
{
    MatDoub xmat;
    VecDoub ytem;
    Int numstr; 
    Int ndim;
    Doub cut;
    
    functionRot(MatDoub inte, VecDoub yt, Int nums, Int di, Doub cu):xmat(inte), ytem(yt), numstr(nums), dim(di), cut(cu) { }
    // note ee[0] ee[1] ee[2] hold the vector of rotation
    // and ee[3] is the angle 0 to 2*pi
    
    void operator()(VecDoub intPoint, VecDoub derivative)
    {
    
        Doub result = 0; // here the final result is stored
        Doub distance = 0; // the distance between the interpolated point and the point on the grid
        Doub sfunc = 0, dc; // the scalling function for weight function
        Doub eps = 1.e-12; // small parameter for weighting function
        int i,j,k;
        
        
        // 2D case scenario => our xmat has 2 columns
        if (ndim == 2) 
        {
            Int rowb = 8;
            MatDoub B(numstr,rowb);
            VecDoub weights(numstr);
            
            // assemble the B matrix representing the function:
            //  a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + a6*x^3 + a7*y^3
            for (i=0; i<numstr; i++) 
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
            for (i=0; i<numstr; i++) 
            {
                distance = sqrt ( 
                                 (intPoint[0] - xmat[i][0])*(intPoint[0] - xmat[i][0]) + 
                                 (intPoint[1] - xmat[i][1])*(intPoint[1] - xmat[i][1])
                                 );
                dc = distance/cut; // dc is just a dummy variable to simplify equations
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
            
            MatDoub W(numstr,numstr);
            for (i=0; i<numstr; i++) 
            {
                W[i][i] = weights[i];
            }
            // assembly of W matrix finished
            
            // solve the following system: B' * W * B * a' = B' * W * f'
            // B' * W first: B' is 1 x numstr and W is numstr numstr
            MatDoub Btrans(rowb,numstr);
            MatDoub A1(rowb,numstr);
            MatDoub A(rowb,rowb);
            
            transpose(B,Btrans,numstr,rowb);
            matmult(Btrans,W,A1,rowb,numstr,numstr,numstr);
            matmult(A1,B,A,rowb,numstr,numstr,rowb);
            // so the left hand side is stored in A, now we need rhs:
            
            VecDoub y(rowb);
            MatDoub y1(rowb,numstr);
            
            matmult(Btrans,W,y1,rowb,numstr,numstr,numstr);
            // vector matrix multiplication
            for (k=0; k<rowb; k++) 
            {
                y[k] = 0;
                for (i=0; i<numstr; i++) 
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
            /*     result = 0;
             result = a[0] + a[1]*intPoint[0] + a[2]*intPoint[1] + a[3]*intPoint[0]*intPoint[0] +
             a[4]*intPoint[1]*intPoint[1] + a[5]*intPoint[0]*intPoint[1] +
             a[6]*intPoint[0]*intPoint[0]*intPoint[0] + a[7]*intPoint[1]*intPoint[1]*intPoint[1];
             */
            
            derivative[0] = a[1] + 2*a[3]*intPoint[0] + a[5]*intPoint[1] + 3*a[6]*intPoint[0]*intPoint[0];
            
            derivative[1] = a[2] + 2*a[4]*intPoint[1] + a[5]*intPoint[0] + 3*a[7]*intPoint[1]*intPoint[1];
            
            
            
        }
        // =========    2D case is done ====================
        
        
        
        // 3D case scenario => our xmat has 3 columns 
        if (ndim == 3) 
        {
            Int rowb=7;
            MatDoub B(numstr,rowb);
            VecDoub weights(numstr);
            
            // assemble the B matrix representing the function:
            for (i=0; i<numstr; i++) 
            {
                B[i][0] = 1.0; // a0
                B[i][1] = xmat[i][0]; // a1*x
                B[i][2] = xmat[i][1]; // a2*y
                B[i][3] = xmat[i][2]; // a3*z
                B[i][4] = xmat[i][0]*xmat[i][0]; // a4*x^2
                B[i][5] = xmat[i][1]*xmat[i][1]; // a5*y^2
                B[i][6] = xmat[i][2]*xmat[i][2]; // a6*z^2
                //   B[i][7] = xmat[i][0]*xmat[i][1]; // a7*x*y
                //   B[i][8] = xmat[i][0]*xmat[i][2]; // a8*x*z
                //   B[i][9] = xmat[i][1]*xmat[i][2]; // a9*y*z
                /*  B[i][10] = xmat[i][0]*xmat[i][0]*xmat[i][0];
                 B[i][11] = xmat[i][1]*xmat[i][1]*xmat[i][1];
                 B[i][12] = xmat[i][2]*xmat[i][2]*xmat[i][2];
                 B[i][13] = xmat[i][0]*xmat[i][0]*xmat[i][1];
                 B[i][14] = xmat[i][0]*xmat[i][0]*xmat[i][2];
                 B[i][15] = xmat[i][1]*xmat[i][1]*xmat[i][0];
                 B[i][16] = xmat[i][1]*xmat[i][1]*xmat[i][2];
                 B[i][17] = xmat[i][2]*xmat[i][2]*xmat[i][0];
                 B[i][18] = xmat[i][2]*xmat[i][2]*xmat[i][1];
                 B[i][19] = xmat[i][0]*xmat[i][1]*xmat[i][2];
                 */
            } // Matrix B assembled
            
             MatDoub W(numstr,numstr);
            // assemble the weight matrix
            for (i=0; i<numstr; i++) 
            {
                distance = sqrt ( 
                                 (intPoint[0] - xmat[i][0])*(intPoint[0] - xmat[i][0]) + 
                                 (intPoint[1] - xmat[i][1])*(intPoint[1] - xmat[i][1]) +
                                 (intPoint[2] - xmat[i][2])*(intPoint[2] - xmat[i][2])
                                 );
                dc = distance/cut; // dc is just a dummy variable to simplify equations
                if (dc > 1) 
                {
                    sfunc = 0;
                }
                else
                {
                    sfunc = (1 - dc*dc*dc*dc*dc*dc);
                }
                W[i][i] = sfunc * exp (-(distance*distance) ) / ( pow(distance,8) + eps);
                
            }
            // assembly of W matrix finished
            
            // solve the following system: B' * W * B * a' = B' * W * f'
            // B' * W first: B' is 1 x numstr and W is numstr numstr
            MatDoub Btrans(rowb,numstr);
            MatDoub A1(rowb,numstr);
            MatDoub A(rowb,rowb);
            
            transpose(B,Btrans,numstr,rowb);
            matmult(Btrans,W,A1,rowb,numstr,numstr,numstr);
            matmult(A1,B,A,rowb,numstr,numstr,rowb);
            // so the left hand side is stored in A, now we need rhs:
            
            VecDoub y(rowb);
            MatDoub y1(rowb,numstr);
            
            matmult(Btrans,W,y1,rowb,numstr,numstr,numstr);
            for (k=0; k<rowb; k++) 
            {
                y[k] = 0;
                for (i=0; i<numstr; i++) 
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
            //    result = 0;
            //   result = a[0] + a[1]*intPoint[0] + a[2]*intPoint[1] + a[3]*intPoint[2] +
            //             a[4]*intPoint[0]*intPoint[0] + a[5]*intPoint[1]*intPoint[1] +
            //             a[6]*intPoint[2]*intPoint[2] + a[7]*intPoint[0]*intPoint[1] +
            //    a[8]*intPoint[0]*intPoint[2] + a[9]*intPoint[1]*intPoint[2];
            /*
             a[10]*intPoint[0]*intPoint[0]*intPoint[0] +   
             a[11]*intPoint[1]*intPoint[1]*intPoint[1] +
             a[12]*intPoint[2]*intPoint[2]*intPoint[2] +
             a[13]*intPoint[0]*intPoint[0]*intPoint[1] +
             a[14]*intPoint[0]*intPoint[0]*intPoint[2] +
             a[15]*intPoint[1]*intPoint[1]*intPoint[0] +
             a[16]*intPoint[1]*intPoint[1]*intPoint[2] +
             a[17]*intPoint[2]*intPoint[2]*intPoint[0] +
             a[18]*intPoint[2]*intPoint[2]*intPoint[1] +
             a[19]*intPoint[0]*intPoint[1]*intPoint[2];
             */
            
            derivative[0] = a[1] + 2*a[4]*intPoint[0]; //+ a[7]*intPoint[1] + a[8]*intPoint[2];
            
            derivative[1] = a[2] + 2*a[5]*intPoint[1]; //+ a[7]*intPoint[0] + a[9]*intPoint[2];
            
            derivative[2] = a[3] + 2*a[6]*intPoint[2]; // + a[8]*intPoint[0] + a[9]*intPoint[1];
    
    
    }
    
};
