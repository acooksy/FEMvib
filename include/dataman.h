/*
 Peter Zajac 28th of December 2010
 
 
	this subroutine reads in the vector of functional values and adjusts them according to 
	their properties: if the numbers are very similar it subtracts the first common digists 
	and multiplies the rest of the numbers hence we obtain a grid extention.
 
	oldvec - the original data
	newvec - the new vector of values
	vecsize - size of the vectors (size of array of functional numbers)
	meanval - mean value of the array
	stdval - standard deviation of the array
	whereweare - will denote the transformation that was used - this number is passed to the program so
					it knows how to adjust the values back to its orgininal format
					1 - 1000,0.1		1
					2 - 1000,1.0		1
					3 - 1000,10.0		1
					4 - 1000,100.0		1	
					5 - 1000,other		0
					
					6 - 100,0.1			1
					7 - 100,1.0			1
					8 - 100,10.0		1
					9 - 100,other		0
					
					10 - 10,0.1			1
					11 - 10,1.0			1
					12 - 10, other		0
			
					13 - 1,0.1			1
					14 - 1,other		0
 
					15 - other,other	0
 
					
 
 
 
 */

void functionvalsmanip(VecDoub &oldvec, VecDoub &newvec, Int &vecsize, Doub &meanval, Doub &stdval, Int &whereweare, Int &valuesubtracted)

{
	Int i;
	valuesubtracted = 0;
	whereweare = 0;
	
// lets look how big the number is:
	Doub criterion;
	criterion = abs(meanval/1000.0);
	
// IF THOUSANDS THEN	
	if (criterion >= 1.0) 
	{
		//VAL
		if (stdval <= 0.1)
		{
			valuesubtracted = (int) (meanval);
			
			
			for (i=0; i<vecsize; i++) 
			{
				newvec[i] = (oldvec[i] - valuesubtracted)*100;
			}
			whereweare = 1;
			
		}
		
		//VAL
		if (stdval > 0.1 && stdval <= 1.0) 
		{
			valuesubtracted = (int) (meanval/10);
			
			for (i=0; i<vecsize; i++) 
			{
				newvec[i] = (oldvec[i] - valuesubtracted)*100;
			}
			whereweare = 1;
		}
		
		//VAL
		if (stdval > 1.0 && stdval <= 10.0) 
		{
			valuesubtracted = (int) (meanval/100);
			
			for (i=0; i<vecsize; i++) 
			{
				newvec[i] = (oldvec[i] - valuesubtracted)*100;
			}
			whereweare = 1;
		}
		
		//VAL
		if (stdval > 10.0 && stdval <= 100.0) 
		{
			valuesubtracted = (int) (meanval/1000);
			
			for (i=0; i<vecsize; i++) 
			{
				newvec[i] = (oldvec[i] - valuesubtracted)*100;
			}
			whereweare = 1;
		}
		
		//FINAL
		if (stdval > 100.0) 
		{
			for (i=0;i<vecsize;i++)
			{
				newvec[i] = oldvec[i];
			}
			whereweare = 0;
		}
		
	}

// IF HUNDREDS THEN:
	if (criterion < 1.0 && criterion >= 0.1) 
	{
		//VAL
		if (stdval <= 0.1)
		{
			valuesubtracted = (int) (meanval);
			
			for (i=0; i<vecsize; i++) 
			{
				newvec[i] = (oldvec[i] - valuesubtracted)*100;
			}
			whereweare = 1;
			
		}
		
		//VAL
		if (stdval > 0.1 && stdval <= 1.0) 
		{
			valuesubtracted = (int) (meanval/10);
			
			for (i=0; i<vecsize; i++) 
			{
				newvec[i] = (oldvec[i] - valuesubtracted)*100;
			}
			whereweare = 1;
		}
		
		//VAL
		if (stdval > 1.0 && stdval <= 10.0) 
		{
			valuesubtracted = (int) (meanval/100);
			
			for (i=0; i<vecsize; i++) 
			{
				newvec[i] = (oldvec[i] - valuesubtracted)*100;
			}
			whereweare = 1;
		}
		
		//FINAL
		if (stdval > 10.0) 
		{
			for (i=0;i<vecsize;i++)
			{
				newvec[i] = oldvec[i];
			}
			whereweare = 0;
		}
	}

	
// IF TENS THEN:
	if (criterion < 0.1 && criterion >= 0.01) 
	{
		//VAL
		if (stdval <= 0.1)
		{
			valuesubtracted = (int) (meanval);
			
			for (i=0; i<vecsize; i++) 
			{
				newvec[i] = (oldvec[i] - valuesubtracted)*100;
			}
			whereweare = 1;
			
		}
		
		//VAL
		if (stdval > 0.1 && stdval <= 1.0) 
		{
			valuesubtracted = (int) (meanval/10);
			
			for (i=0; i<vecsize; i++) 
			{
				newvec[i] = (oldvec[i] - valuesubtracted)*100;
			}
			whereweare = 1;
		}
		
		//FINAL
		if (stdval > 1.0) 
		{
			for (i=0;i<vecsize;i++)
			{
				newvec[i] = oldvec[i];
			}
			whereweare = 0;
		}
	}

	
// IF ONES THEN:
	if (criterion < 0.01 && criterion >= 0.001) 
	{
		//VAL
		if (stdval <= 0.1)
		{
			valuesubtracted = (int) (meanval);
			
			for (i=0; i<vecsize; i++) 
			{
				newvec[i] = (oldvec[i] - valuesubtracted)*100;
			}
			whereweare = 1;
			
		}
		
		//FINAL
		if (stdval > 0.1) 
		{
			for (i=0;i<vecsize;i++)
			{
				newvec[i] = oldvec[i];
			}
			whereweare = 0;
		}
	}

// IF DECIMAL THEN:
	if (criterion < 0.001) 
	{
		for (i=0;i<vecsize;i++)
		{
			newvec[i] = oldvec[i];
		}
		whereweare = 0;
	}

// end of the subroutine

}


