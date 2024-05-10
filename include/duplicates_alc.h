
// checks for duplicates and prints warning if they exist

void dupl(MatDoub &x, Int &npt, Int &dim)
{


int i,j,k;


if (dim == 3)
{
for (i=0;i<npt;i++)
{
  for (j=i+1;j<npt;j++)
  {
    if (x[i][0] == x[j][0] && x[i][1] == x[j][1] && x[i][2] == x[j][2])
	{
	std::cout << "Warning: Duplicates exist! Check your surface. \n";
	
	Doub rad = 0.0174532925;
	Doub deg = 90.0;

	std::cout << x[j][0] << " " << x[j][1]/rad - deg << " " << x[j][2]/rad - deg << "\n";
	std::cout << i << "\n"; 

	}
  }
	
}
}


}

