

double calcvalkrig(MatDoub &xmat, VecDoub &yvec, Int &numpoints, Int &point, Doub &rval, Int &nd, Int &npts)
{
	Int i, j=0,k;
	Doub vys=0;
	Int newpoints = numpoints - 1;

	
	VecDoub ytemp(newpoints);
	MatDoub xtemp(newpoints,nd);
	VecDoub xstar(nd);
	
	// Initialize 
	for (i=0;i<numpoints;i++)
	{
		
		if (i==point) 
		{
		   for (k=0;k<nd;k++)
			{
			xstar[k]=xmat[i][k];
			}
			j--;
		}
		
		else 
		{
			ytemp[j]=yvec[i];
			
		   for (k=0;k<nd;k++)
			{
			xtemp[j][k] = xmat[i][k];
			}	
		}
		
		j++;
	}

	// Sort the closest 50 values	
        MatDoub pts(npts,nd);
        VecDoub yvecnew(npts); 
	howmany(xtemp, ytemp, xstar, nd, newpoints, npts, pts, yvecnew);
	// Calculate Kriging Interpolation
	Powvargram vgram(pts,yvecnew,rval);
	Krig<Powvargram> krigp(pts,yvecnew,vgram);
	vys=krigp.interp(xstar);
	
	return vys;
	
}


void krigcall(Int &ndim, Doub &r, Int &testpt)
{

int npts,i;
Doub resultofinterp,err;
ifstream indata;

	 indata.open("results/energyP");
		if (!indata) 
			{
				std::cerr << "Error: file could not be opened\n";
				exit(1);
			}
	
 		indata >> npts;
		MatDoub x(npts,ndim); 
		VecDoub y(npts);
		VecDoub com(npts);
		if (ndim == 3)
		{
		for (i=0;i<npts;i++)
		{
			cout.precision(10);
			cout.setf(ios::fixed,ios::floatfield);
			indata >> x[i][0] >> x[i][1] >> x[i][2] >> y[i];
		}
		}

		if (ndim == 2)
		{
		for (i=0;i<npts;i++)
		{
			cout.precision(10);
			cout.setf(ios::fixed,ios::floatfield);
			indata >> x[i][0] >> x[i][1] >> y[i];
		}
		}

	indata.close(); 
	// close the input file ===============================================
	

	// =========== KRIGING INTERPOLATION ========================================================
			Doub aver = 0.0;
			for (i=0;i<npts;i++)
				{
					resultofinterp = calcvalkrig(x, y, npts, i, r, ndim, testpt);
					com[i] = y[i] - resultofinterp;
					aver+=abs(com[i]);
					cout << com[i] << "\n";
				}
			cout << "\n";
			cout << testpt << " " << r << " " << aver/npts << "\n";
	// ============== END KRIGING INTEPOLATION ====================================================
	
//	VecInt ind(npts);
//	selecty(com,npts,ind);

//	Doub rad = 0.0174532925;
//	Doub deg = 90.0;

//	for (i=npts-100;i<npts;i++)
//	{
//	cout << com[ind[i]] << " " << x[ind[i]][0] << " " << x[ind[i]][1]/rad - deg << " " << x[ind[i]][2]/rad - deg << "\n";
//	}

}


