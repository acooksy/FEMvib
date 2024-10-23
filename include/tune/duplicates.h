
void duplicates(MatDoub &x, Int &npt, Int &dim)
{


int i,j,k;


if (dim == 3)
{
for (i=0;i<npts;i++)
{
  for (j=i+1;j<npts;j++)
  {
    if (x[i][0] == x[j][0] && x[i][1] == x[j][1] && x[i][2] == x[j][2])
	{
	cout << "Warning: Duplicates exist! Interpolatio might not be accurate. \n";
	}
  }
	
}
}



  return 0;

}

