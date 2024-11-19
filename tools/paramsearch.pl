#! /usr/bin/perl


$file_in = $ARGV[0];
$file_opp = "results/opp";
$file_optr = "results/optr";
	
$cou = 0; 
	 open(INFILE25, $file_in);
		while (<INFILE25>)
			{
			@linenum = split();
			$points = @linenum; 
			if ($points == 3)
				{
				$cou++;
				}	

			}
	 close(INFILE25);


	 open(INFILE30, $file_in);		
  	 while (<INFILE30>)
	       {
               @ln = split();  # number of elements on a single line
               $number = @ln;
		
               if ($number == 3)
                  {	
                   $arrayp = $arrayp . " " . $ln[0];
                   $arrayr = $arrayr . " " . $ln[1];
                   $arraye = $arraye . " " . $ln[2];
		  }
               }
	close(INFILE30);

		$arrayp = substr($arrayp,1); 
		@valuesp = split(' ',$arrayp);
		
		$arrayr = substr($arrayr,1); 
		@valuesr = split(' ',$arrayr);

		$arraye = substr($arraye,1); 
		@valuese = split(' ',$arraye);
       
	$min = $valuese[0];
	$index = 0;
	for ($j=1;$j<$cou;$j++)
	{
		if ($min > $valuese[$j])
		{
			$min = $valuese[$j];
			$index = $j;
		}	

	}
	# we know the minimum stored in $min we need to find corresponding r
		print("$valuesr[$index] \n$valuesp[$index]");
		open OUTPUT, '>', $file_opp;
		print OUTPUT $valuesp[$index];
		close(OUTPUT);
		open OUTPUT, '>', $file_optr;
		print OUTPUT $valuesr[$index];
		close(OUTPUT);








