//
//  helpf.h
//  ILMS
//
//  Created by PETER ZAJAC on 9/28/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//



// findCut finds the cut off parameter for a particular grid
void findCut(MatDoub &xmat, Int &numpoints, Int &ndim, Doub &cut)
{
	Int i;
	Doub maxX, maxY, maxZ;
	Doub minX, minY, minZ;
	Doub xcord,ycord,zcord;
    Doub grid;
    
  /*  if (ndim == 3)
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
        
        
        cut = sqrt(xcord * xcord + ycord * ycord + zcord * zcord);
        
    }
   */ 
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
        
        
        xcord = (maxX - minX)/3;
        ycord = (maxY - minY)/3;
        
        
        cut = sqrt(xcord * xcord + ycord * ycord);
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
        
       cut = (maxX - minX)/5;
        
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
            cout << matTpPrint[i][j] << " ";
        }
        cout << "\n";
    } 
    
}

void printGmat(MatDoub &gmatrix, MatDoub &internals, Int &dim, Int &numstructs)
{
    int i,j,k;
    
    // print the 2D matrix
    if (dim == 2)
    {
        for (i=0; i<numstructs; i++) 
        {
            cout << internals[i][0] << " " << internals[i][1] << " ";
            for (j=0; j<dim*dim; j++) 
            {
                cout << setw(25) << gmatrix[i][j] << " ";
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
            }
            cout << "\n";
            
        }
    }
}

// analytical gmatrix representation:
void gmatAnalytical(MatDoub &gmat, MatDoub &internals, VecDoub &mas, Int &structures)
{
    // gmat numstr X dim*dim   
    int i;
    
    // through structres
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


