//
//  pointSelector.h
//  GMAT
//
//  Created by PETER ZAJAC on 10/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
void howmany(MatDoub &x, VecDoub &y, VecDoub &point, Int &ndim, Int &numpoints, Int &pocet, MatDoub &pts, VecDoub &ynew)
{
    int i,j,k,pos;
    VecDoub source(numpoints);
    VecInt unused(numpoints);
    
    // Careful here VecDoub ind are the points we will take in
    
    for (i=0; i<numpoints; i++)
    {
        source[i] = 0.0; // set distance to zero
        for (j=0; j<ndim; j++)
        {
            source[i] += (x[i][j]-point[j])*(x[i][j]-point[j]);
        }
        source[i] = sqrt(source[i]); // Catesian distance from point p
        unused[i] = 1;
    }
    
    
    for (i = 0; i < pocet; i++)
    {
        pos = -1;
        for (j = 0; j < numpoints; j++)
        {
            if (pos == -1)
            {
                if (unused[j])
                {
                    pos = j;
                }
            }
            else
            {
                if (unused[j] && (source[j] < source[pos]))
                {
                    pos = j;
                }
            }
        }
        
        for (k=0;k<ndim;k++)
        {
            pts[i][k] = x[pos][k];
        }
        ynew[i] = y[pos];
        unused[pos] = 0;
    }
}


void howmany_noref(MatDoub &x, VecDoub &point, Int &ndim, Int &numpoints, Int &pocet, MatDoub &pts)
{
    int i,j,k,pos;
    VecDoub source(numpoints);
    VecInt unused(numpoints);
    
    // Careful here VecDoub ind are the points we will take in
    
    for (i=0; i<numpoints; i++)
    {
        source[i] = 0.0; // set distance to zero
        for (j=0; j<ndim; j++)
        {
            source[i] += (x[i][j]-point[j])*(x[i][j]-point[j]);
        }
        source[i] = sqrt(source[i]); // Catesian distance from point p
        unused[i] = 1;
    }
    
    
    for (i = 0; i < pocet; i++)
    {
        pos = -1;
        for (j = 0; j < numpoints; j++)
        {
            if (pos == -1)
            {
                if (unused[j])
                {
                    pos = j;
                }
            }
            else
            {
                if (unused[j] && (source[j] < source[pos]))
                {
                    pos = j;
                }
            }
        }
        
        for (k=0;k<ndim;k++)
        {
            pts[i][k] = x[pos][k];
        }
        unused[pos] = 0;
    }
}


void ordering(MatDoub &x, VecInt &y, MatDoub &point, Int &ndim, Int &numpoints, Int &pocet)
{
    int i,j,k,pos;
    VecDoub source(numpoints);
    VecInt unused(numpoints);
    
    // Careful here VecDoub ind are the points we will take in
    
    for (i=0; i<numpoints; i++)
    {
        source[i] = 0.0; // set distance to zero
        for (j=0; j<ndim; j++)
        {
            source[i] += (x[i][j]-point[0][j])*(x[i][j]-point[0][j]);
        }
        source[i] = sqrt(source[i]); // Catesian distance from point p
        unused[i] = 1;
    }
    
    
    for (i = 0; i < pocet; i++)
    {
        pos = -1;
        for (j = 0; j < numpoints; j++)
        {
            if (pos == -1)
            {
                if (unused[j])
                {
                    pos = j;
                }
            }
            else
            {
                if (unused[j] && (source[j] < source[pos]))
                {
                    pos = j;
                }
            }
        }
        

        y[i] = pos;
        unused[pos] = 0;
    }
}



int differ(Doub &a, Doub &b)
{
    int l=0;
    if (a > 0 && b < 0)
    {
        l=1;
    }
    if (a < 0 && b > 0) 
    {
        l=1;
    }
    if ((abs(a)-abs(b)) < 0.001) 
    {
        l=1;
    }
    
    return l;
    
    
}
