
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
// 	cout << source[pos] << " ";
       unused[pos] = 0;
    }
// cout << "\n";
}

void selecty(VecDoub &y, Int &numpoints, VecInt &selected)
{
      int i,j,k,pos;
      VecDoub source(numpoints);			   
      VecInt unused(numpoints);

    for (i=0;i<numpoints; i++)
	{
	source[i] = abs(y[i]);
	unused[i] = 1;
 	}


 
// Careful here VecDoub ind are the points we will take in
    for (i = 0; i < numpoints; i++) 
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

	selected[i] = pos;
        unused[pos] = 0;
    }
}

