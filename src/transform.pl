#! /usr/bin/perl

# perl script for extracting the values from the gmat and submitting them to the c-code for interpolation

# written by Peter Zajac during January 2011
# modified by Peter Zajac during October 2011 


$filename = $ARGV[0];
$dimension = $ARGV[1];
$hartree = $ARGV[2];

# print "the name is: $filename \n";

if ($dimension == 3)
  {	 
	
$cou = 0; 
	 open(INFILE25, $filename);
		while (<INFILE25>)
			{
			@linenum = split();
			$points = @linenum; 
			if ($points == 4)
				{
				
				$array = $array . " " . $linenum[3];
				$cou++;
				}	

			}
	 close(INFILE25);
	 # the data are in

		$array = substr($array,1); 
		@values = split(' ',$array);

	# find the minimum on the surface and determine its distance from 0

	$minsurf=$values[0];
	for ($j=1;$j<$cou;$j++)
	{
		if ($values[$j] < $minsurf) {$minsurf = $values[$j]; }
	}

	$difference = 0 - $minsurf;

	# add the difference to the surface so that it stays positive and touching 0			
       for ($j=0;$j<$cou;$j++)
           {
		$values[$j] = $values[$j] + $difference;
  	   }
	
	# if the surface is in hartrees then convert to cm-1
	if ($hartree == 1)
	{
        for ($j=0;$j<$cou;$j++)
           {
		$values[$j] = $values[$j]*219474.63;
  	   }
	}
 	 
	 open(INFILE1, $filename); 
 	 open(OUTFILE4, ">energyP"); 

         $g=0;
  	 while (<INFILE1>)
	       {
		@lnn=split();
		$points=@lnn;

		if ($points == 4)
			{
			print(OUTFILE4 "$lnn[0] $lnn[1] $lnn[2] $values[$g]  \n");
			$g++;
			}
                }
	close(INFILE1); close(OUTFILE4);
}


if ($dimension == 2)
  {	 
	
$cou = 0; 
	 open(INFILE25, $filename);
		while (<INFILE25>)
			{
			@linenum = split();
			$points = @linenum; 
			if ($points == 3)
				{
				
				$array = $array . " " . $linenum[2];
				$cou++;
				}	

			}
	 close(INFILE25);
	 # the data are in

		$array = substr($array,1); 
		@values = split(' ',$array);

	# find the minimum on the surface and determine its distance from 0

	$minsurf=$values[0];
	for ($j=1;$j<$cou;$j++)
	{
		if ($values[$j] < $minsurf) {$minsurf = $values[$j]; }
	}

	$difference = 0 - $minsurf;

	# add the difference to the surface so that it stays positive and touching 0			
       for ($j=0;$j<$cou;$j++)
           {
		$values[$j] = $values[$j] + $difference;
  	   }
	
	# if the surface is in hartrees then convert to cm-1
	if ($hartree == 1)
	{
        for ($j=0;$j<$cou;$j++)
           {
		$values[$j] = $values[$j]*219474.63;
  	   }
	}
 	 
	 open(INFILE1, $filename); 
 	 open(OUTFILE4, ">energyP"); 

         $g=0;
  	 while (<INFILE1>)
	       {
		@lnn=split();
		$points=@lnn;

		if ($points == 3)
			{
			print(OUTFILE4 "$lnn[0] $lnn[1] $values[$g]  \n");
			$g++;
			}
                }
	close(INFILE1); close(OUTFILE4);
}


if ($dimension == 1)
  {	 
	
$cou = 0; 
	 open(INFILE25, $filename);
		while (<INFILE25>)
			{
			@linenum = split();
			$points = @linenum; 
			if ($points == 2)
				{
				
				$array = $array . " " . $linenum[1];
				$cou++;
				}	

			}
	 close(INFILE25);
	 # the data are in

		$array = substr($array,1); 
		@values = split(' ',$array);

	# find the minimum on the surface and determine its distance from 0

	$minsurf=$values[0];
	for ($j=1;$j<$cou;$j++)
	{
		if ($values[$j] < $minsurf) {$minsurf = $values[$j]; }
	}

	$difference = 0 - $minsurf;

	# add the difference to the surface so that it stays positive and touching 0			
       for ($j=0;$j<$cou;$j++)
           {
		$values[$j] = $values[$j] + $difference;
  	   }
	
	# if the surface is in hartrees then convert to cm-1
	if ($hartree == 1)
	{
        for ($j=0;$j<$cou;$j++)
           {
		$values[$j] = $values[$j]*219474.63;
  	   }
	}
 	 
	 open(INFILE1, $filename); 
 	 open(OUTFILE4, ">energyP"); 

         $g=0;
  	 while (<INFILE1>)
	       {
		@lnn=split();
		$points=@lnn;

		if ($points == 2)
			{
			print(OUTFILE4 "$lnn[0]  $values[$g]  \n");
			$g++;
			}
                }
	close(INFILE1); close(OUTFILE4);
}


