#! /usr/bin/perl

# perl script for extracting the values from the gmat and submitting them to the c-code for interpolation

# written by Peter Zajac during January 2011
# modified by Peter Zajac during October 2011 
# modified ALC 3-Oct-2024 to allow path as input argument and to put matrixg files in results/
# modified ALC 19-Nov-2024 to add number of entries to top of file

$filename = $ARGV[0];
$dimension = $ARGV[1];
$path = $ARGV[2];

$name = $path."/results/matrixg.";
# print "the input filename is: $filename \n";

if ($dimension == 3)
  {	 
	 $newname3 = $name."xx"; 
	 $newname4 = $name."xy"; 
	 $newname7 = $name."yy";
         $newname9 = $name."zx"; 
	 $newname10 = $name."zy"; 
	 $newname11 = $name."zz";
	
$cou = 0; 
	 open(INFILE25, $filename);
		while (<INFILE25>)
			{
			@linenum = split();
			$points = @linenum; 
			if ($points == 12)
				{
				$cou++;
				}	

			}
	 close(INFILE25);


	 open(INFILE30, $filename);		
  	 while (<INFILE30>)
	       {
               @ln = split();  # number of elements on a single line
               $number = @ln;
		
               if ($number == 12)
                  {	
                        $array4 = $array4 . " " . $ln[3];
                        $array5 = $array5 . " " . $ln[4];
                        $array8 = $array8 . " " . $ln[7];
                        $array10 = $array10 . " " . $ln[9];
                        $array11 = $array11 . " " . $ln[10];
                        $array12 = $array12 . " " . $ln[11];
		}
               }
	close(INFILE30);

		$array4 = substr($array4,1); 
		$array5 = substr($array5,1); 
		$array8 = substr($array8,1); 
		$array10 = substr($array10,1);
		$array11 = substr($array11,1); 
		$array12 = substr($array12,1); 
		
		@values4 = split(' ',$array4);
		@values5 = split(' ',$array5);
		@values8 = split(' ',$array8);
		@values10 = split(' ',$array10);
		@values11 = split(' ',$array11);
		@values12 = split(' ',$array12);

       for ($j=0;$j<$cou;$j++)
           {
		$values4[$j] = abs($values4[$j]);	
		$values5[$j] = abs($values5[$j]);	
		$values8[$j] = abs($values8[$j]);	
		$values10[$j] = abs($values10[$j]);	
		$values11[$j] = abs($values11[$j]);	
		$values12[$j] = abs($values12[$j]);	

	   }
	
	


	$middle = int $cou/2;
	@sor4 = sort {$a <=> $b} @values4; $final4 = $sor4[$middle];
	@sor5 = sort {$a <=> $b} @values5; $final5 = $sor5[$middle];
	@sor8 = sort {$a <=> $b} @values8; $final8 = $sor8[$middle];
	@sor10 = sort {$a <=> $b} @values10; $final10 = $sor10[$middle];
	@sor11 = sort {$a <=> $b} @values11; $final11 = $sor11[$middle];
	@sor12 = sort {$a <=> $b} @values12; $final12 = $sor12[$middle];


	 open(INFILE1, $filename); 
	 open(OUTFILE4, ">$newname3"); 
	 open(OUTFILE5, ">$newname4");	
         open(OUTFILE8, ">$newname7"); 
	 open(OUTFILE10, ">$newname9"); 
	 open(OUTFILE11, ">$newname10");	
         open(OUTFILE12, ">$newname11"); 

  	 while (<INFILE1>)
	       {
               @lnn = split();  # number of elements on a single line
               $numbern = @lnn;
		
		$m=3;
               if ($numbern == 12)
                  {	
			if (abs($lnn[3]) < ($final4 + $m*$final4) && $lnn[3] != 0)
                        {print(OUTFILE4 "$lnn[0] $lnn[1] $lnn[2] $lnn[3]\n");} 
			if (abs($lnn[4]) < ($final5 + $m*$final5) && $lnn[4] != 0)
			{print(OUTFILE5 "$lnn[0] $lnn[1] $lnn[2] $lnn[4]\n");} 
			if (abs($lnn[7]) < ($final8 + $m*$final8) && $lnn[7] != 0)
			{print(OUTFILE8 "$lnn[0] $lnn[1] $lnn[2] $lnn[7]\n");} 
			if (abs($lnn[9]) < ($final10 + $m*$final10) && $lnn[9] != 0)
			{print(OUTFILE10 "$lnn[0] $lnn[1] $lnn[2] $lnn[9]\n");} 
			if (abs($lnn[10]) < ($final11 + $m*$final11) && $lnn[10] != 0)
			{print(OUTFILE11 "$lnn[0] $lnn[1] $lnn[2] $lnn[10]\n");} 
			if (abs($lnn[11]) < ($final12 + $m*$final12) && $lnn[11] != 0)
			{print(OUTFILE12 "$lnn[0] $lnn[1] $lnn[2] $lnn[11]\n");} 
		}
               }

	close(INFILE1); close(OUTFILE4); close(OUTFILE5);
        close(OUTFILE8); close(OUTFILE10); close(OUTFILE11);
        close(OUTFILE12);


chmod 0644, $name."xx" or die $!;
chmod 0644, $name."xy" or die $!;
chmod 0644, $name."yy" or die $!;
chmod 0644, $name."zx" or die $!;
chmod 0644, $name."zy" or die $!;
chmod 0644, $name."zz" or die $!;

}


if ($dimension == 2)
{

	 $newname1 = $name."xx";
         $newname2 = $name."xy";
         $newname4 = $name."yy";
	
$cou = 0; 
	 open(INFILE20, $filename);
		while (<INFILE20>)
			{
			@linenum = split();
			$points = @linenum; 
			if ($points == 6)
				{
				$cou++;
				}	

			}
	 close(INFILE20);

 
	open (INFILE30, $filename);
	while (<INFILE30>)
	       {
               @ln = split();  # number of elements on a single line
               $number = @ln; 
		
               if ($number == 6)
                  {	
                        $array3 = $array3 . " " . $ln[2];
                        $array4 = $array4 . " " . $ln[3];
                        $array6 = $array6 . " " . $ln[5];
		  }
               }
	close(INFILE30);


		$array3 = substr($array3,1); # $aver3 = $sum3/$cou;
		$array4 = substr($array4,1); # $aver4 = $sum4/$cou;
		$array6 = substr($array6,1); # $aver6 = $sum6/$cou;
		@values3 = split(' ',$array3);
		@values4 = split(' ',$array4);
		@values6 = split(' ',$array6);


       for ($j=0;$j<$cou;$j++)
           {
		$values3[$j] = abs($values3[$j]);	
		$values4[$j] = abs($values4[$j]);	
		$values6[$j] = abs($values6[$j]);	

	   }


	$middle = int $cou/2;
	@sor3 = sort {$a <=> $b} @values3; $final3 = $sor3[$middle];
	@sor4 = sort {$a <=> $b} @values4; $final4 = $sor4[$middle];
	@sor6 = sort {$a <=> $b} @values6; $final6 = $sor6[$middle];



	 open(INFILE1, $filename);		
         open(OUTFILE1, ">$newname1");	
         open(OUTFILE2, ">$newname2");	
         open(OUTFILE4, ">$newname4");	

  	 while (<INFILE1>)
	       {
               @lnn = split();  # number of elements on a single line
               $numbern = @lnn;
		$m=3;		
               if ($numbern == 6)
                  {	
			if (abs($lnn[2]) < ($final3 + $m*$final3) && $lnn[2] != 0)
                        {print(OUTFILE1 "$lnn[0] $lnn[1] $lnn[2]\n");} 
			if (abs($lnn[3]) < ($final4 + $m*$final4) && $lnn[3] != 0)
                        {print(OUTFILE2 "$lnn[0] $lnn[1] $lnn[3]\n");} 
			if (abs($lnn[5]) < ($final6 + $m*$final6) && $lnn[5] != 0)
                        {print(OUTFILE4 "$lnn[0] $lnn[1] $lnn[5]\n");} 
		}
               }

	close(INFILE1);
        close(OUTFILE1);
        close(OUTFILE2);
        close(OUTFILE4);

chmod 0644, $name."xx" or die $!;
chmod 0644, $name."xy" or die $!;
chmod 0644, $name."yy" or die $!;

}

	
if ($dimension == 1)
{

	 $newname1 = $name."xx";
	
$cou = 0; 
	 open(INFILE20, $filename);
		while (<INFILE20>)
			{
			@linenum = split();
			$points = @linenum; 
			if ($points == 2)
				{
				$cou++;
				}	

			}
	 close(INFILE20);

 
	open (INFILE30, $filename);
	while (<INFILE30>)
	       {
               @ln = split();  # number of elements on a single line
               $number = @ln; 
		
               if ($number == 2)
                  {	
# print("@ln \n");                       
 			$array3 = $array3 . " " . $ln[1];
		  }
               }
	close(INFILE30);


		$array3 = substr($array3,1); # $aver3 = $sum3/$cou;
		@values3 = split(' ',$array3);


       for ($j=0;$j<$cou;$j++)
           {
		$values3[$j] = abs($values3[$j]);	

	   }


	$middle = int $cou/2;
	@sor3 = sort {$a <=> $b} @values3; $final3 = $sor3[$middle];


	 open(INFILE1, $filename);		
         open(OUTFILE1, ">$newname1");	

  	 while (<INFILE1>)
	       {
               @lnn = split();  # number of elements on a single line
               $numbern = @lnn;
		$m=3;		
               if ($numbern == 2)
                  {	
			if (abs($lnn[1]) < ($final3 + $m*$final3) && $lnn[1] != 0)
                        {print(OUTFILE1 "$lnn[0] $lnn[1] \n");} 
		  }
               }

chmod 0644, $name."xx" or die $!;

	close(INFILE1);
        close(OUTFILE1);

}	

print "DIMENSION = $dimension\n"; 
# Add number of entries to top of each file
if ($dimension == 1)
	{
		my @files = <$name"xx">;
		foreach my $file (@files) {
			print "ADDING LINES TO $file\n";
		        open(OUTFILE,$file);
				        my ($lines) = 0;
		        while (<OUTFILE>) {
		                $lines++;
		                }
		        close(OUTFILE);
		        $tempfile = "temp";
		        open(TEMP,">$tempfile");
		        open(OUTFILE,"<$file");
		        print TEMP $lines,"\n" ;
		        while( <OUTFILE> ) {
		                print TEMP $_;
		        }
		        close(OUTFILE);
		        close(TEMP);
		        unlink($file);
		        rename($tempfile,$file)
		}
	}
if ($dimension == 2)
	{
		my @files = <$name"xx" $name"xy" $name"yy">;
		foreach my $file (@files) {
			print "ADDING LINES TO $file\n";
		        open(OUTFILE,$file);
				        my ($lines) = 0;
		        while (<OUTFILE>) {
		                $lines++;
		                }
		        close(OUTFILE);
		        $tempfile = "temp";
		        open(TEMP,">$tempfile");
		        open(OUTFILE,"<$file");
		        print TEMP $lines,"\n" ;
		        while( <OUTFILE> ) {
		                print TEMP $_;
		        }
		        close(OUTFILE);
		        close(TEMP);
		        unlink($file);
		        rename($tempfile,$file)
		}
	}
if ($dimension == 3)
	{
		my @files = <$name"xx" $name"xy" $name"yy" $name"zx" $name"zy" $name"zz">;
		foreach my $file (@files) {
			print "ADDING LINES TO $file\n";
		        open(OUTFILE,$file);
				        my ($lines) = 0;
		        while (<OUTFILE>) {
		                $lines++;
		                }
		        close(OUTFILE);
		        $tempfile = "temp";
		        open(TEMP,">$tempfile");
		        open(OUTFILE,"<$file");
		        print TEMP $lines,"\n" ;
		        while( <OUTFILE> ) {
		                print TEMP $_;
		        }
		        close(OUTFILE);
		        close(TEMP);
		        unlink($file);
		        rename($tempfile,$file)
		}
	}

