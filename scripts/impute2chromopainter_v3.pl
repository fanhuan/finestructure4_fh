#!/usr/bin/perl
## CONVERTS PHASED IMPUTE2 OUTPUT TO CHROMOPAINTER-STYLE INPUT FILES

sub help {
print("CONVERTS PHASED SHAPEIT/IMPUTE2 OUTPUT TO CHROMOPAINTER-STYLE INPUT FILES\n");

print("usage:   perl impute2chromopainter.pl <options> output_filename_prefix < impute_output_file.haps\n");
print("    or:  zcat impute_output_file.haps.gz | perl impute2chromopainter.pl <options> output_filename_prefix\n");

print("where:\n");
print("        (i) input is read from STDIN (pipe or redirect a .haps or .haps.gz file)\n");
print("        (ii) output_filename_prefix = filename prefix for chromopainter input file(s). The suffix \".phase\" is added\n\n");
print("The output, by default, is in CHROMOPAINTER v2 input format.\n");

print("<options>:\n");
print("-legend <file>:     Legend file (from bcftools --haplegendsample), plain or .gz. Required when the .hap file has no SNP metadata columns.\n");
print("-hap <file>:        Hap file, plain or .gz. If omitted, reads from STDIN as before.\n");
print("-J:                 Jitter (add 1) snp locations if snps are not strictly ascending. Otherwise an error is produced.\n");

print("<further options>   NOTE: YOU ONLY NEED THESE OPTIONS FOR BACKWARDS COMPATABILITY!\n");
print("-v1:                Produce output compatible with CHROMOPAINTER v1, i.e. include the line of \"S\" for each SNP. \n");
print("-f:                 By default, this script produces PHASE-style output, which differs from \n");
print("			   ChromoPainter input which requires an additional first line.  This option creates the correct\n");
print("		           first line for standard fineSTRUCTURE usage (i.e. the first line is \"0\", all other lines are appended)\n\n");

print("NOTE: TO USE IN CHROMOPAINTER: You also need a recombination map. Create this with the \"convertrecfile.pl\" or \"makeuniformrecfile.pl\" scripts provided.\n\n");
print(" !!! WARNING:  THIS PROGRAM DOES NOT SUFFICIENTLY CHECK FOR MISSPECIFIED FILES. WE ARE NOT ACCOUNTABLE FOR THIS RUNNING INCORRECTLY !!!\n");
die "\n";
}

use Switch;


###############################
## INPUT:

$numindsmaxsize=500;     ## ONLY READ IN THIS NUMBER OF IMPUTE2 INDS AT A TIME (TO MINIMIZE RAM REQUIREMENTS -- MAKES PROGRAM RUN A FACTOR OF (N/$numindsmaxsize) SLOWER THAN IF YOU SET $numindsmaxsize=N, where N is the total number of individuals in your ".haps" IMPUTE2 output file, but also uses a factor of (N/$numindsmaxsize) less RAM

###############################
## ARGUMENT PROCESSING

$v1=0; ## version 1 mode
$fsmode=0; ## finestructure mode (i.e. start with an additional line containing 0)
$jitter=0; ## whether we jitter snp locations

$Mb = 1000000.0;
$IMPUTEinfile="";
$outfilePRE="";
$legendfile="";
$hapfile="";

$argon=0;
for (my $i = 0; $i < scalar(@ARGV); ++$i){
	if(@ARGV[$i] eq "-f"){
	    $fsmode=1;
	}elsif(@ARGV[$i] eq "-v1"){
	    $v1=1;
	}elsif(@ARGV[$i] eq "-J"){
	    $jitter=1;
	}elsif(@ARGV[$i] eq "-legend"){
	    $i++;
	    $legendfile=$ARGV[$i];
	}elsif(@ARGV[$i] eq "-hap"){
	    $i++;
	    $hapfile=$ARGV[$i];
	}else{
	    switch($argon){
		case 0 {$outfilePRE="$ARGV[$i]";}
		else {
		    help();
		}
	    }
	    $argon++;
	}
}

if($outfilePRE eq "" || $argon != 1) {help();}
$outfilePRE =~ s/.phase$//;

##############################
## PROGRAM:


## (II) GET NUMBER OF SITES AND INDS:

## Read hap file from -hap argument or STDIN
if ($hapfile ne "") {
    if ($hapfile =~ /\.gz$/) {
        open(HAPFH, "zcat \Q$hapfile\E |") or die "Cannot open hap file '$hapfile': $!\n";
    } else {
        open(HAPFH, "<", $hapfile) or die "Cannot open hap file '$hapfile': $!\n";
    }
    @INlines = <HAPFH>;
    close(HAPFH);
} else {
    @INlines = <STDIN>;
}

## If a legend file is given, prepend its SNP metadata to each hap line so the
## existing +5 offset logic below works unchanged.  Legend format (bcftools):
##   id position a0 a1   (header on line 1, data from line 2)
## We prepend a dummy field so indices match what the script expects:
##   [0]=X  [1]=id  [2]=position  [3]=a0  [4]=a1  [5..]=haplotypes
if ($legendfile ne "") {
    if ($legendfile =~ /\.gz$/) {
        open(LEGFH, "zcat \Q$legendfile\E |") or die "Cannot open legend file '$legendfile': $!\n";
    } else {
        open(LEGFH, "<", $legendfile) or die "Cannot open legend file '$legendfile': $!\n";
    }
    <LEGFH>; ## skip legend header line
    for my $idx (0..$#INlines) {
        my $legline = <LEGFH>;
        die "Legend file has fewer lines than hap file (at hap line ".($idx+1).")\n" unless defined $legline;
        chomp $legline;
        chomp $INlines[$idx];
        $INlines[$idx] = "X $legline $INlines[$idx]\n";
    }
    close(LEGFH);
}
$line=$INlines[0];
@linearray=split(/\s+/,$line);
$totalINDS=(@linearray-5)/2;
$totalhaps=@linearray-5;
$numsplits=int($totalINDS/$numindsmaxsize);
if (($numsplits*$numindsmaxsize)<$totalINDS)
{
    $numsplits=$numsplits+1;
}
$nsites=scalar(@INlines);
#print "$numsplits $nsitesFULL $totalINDS\n";

              ## (III) READ IN IMPUTE2 HAPLOTYPES AND MAKE CHROMOPAINTER HAPLOTYPE INPUT FILE:
open(OUT,">${outfilePRE}.phase");
if($fsmode==1) {
	print OUT "0\n";
}
for ($a=0; $a < $numsplits; $a+=1)
{
    $startIND=$a*$numindsmaxsize;
    $endIND=($a+1)*$numindsmaxsize;
    if ($endIND > $totalINDS)
    {
	$endIND=$totalINDS;
    }

                              ## read in:
    @rsvec=();
    @posvec=();
    @genomat=();
    $snpcount=0;
    foreach $line (@INlines)
    {
	@linearray=split(/\s+/,$line);
	push(@rsvec,$linearray[1]);
#	print("BEFORE: SNP location $linearray[2] lup $lastuniquesnp\n");
	if(scalar(@posvec)>0 &&($linearray[2] <= $posvec[-1]) && $linearray[2]>=0){
	    if(!$jitter){
		die("ERROR: SNPs are not strictly ascending, exiting. Rerun with -J to jitter the SNP locations.\n");
	    }
#	    print("Duplication found: setting $linearray[2] to $posvec[-1]+1\n");
	    $linearray[2]=$posvec[-1]+1;
	}
#	print("AFTER:  SNP location $linearray[2] lup $lastuniquesnp\n");
	if( $linearray[2] == $posvec[-1] ){
	    die("Strange error due to jittering?\n");
	}
	push(@posvec,$linearray[2]);
	# shift(@linearray);
	# shift(@linearray);
	# shift(@linearray);
	# shift(@linearray);
	# shift(@linearray);

	for ($i=$startIND; $i < $endIND; $i+=1)
	{
	    $genomat[(($i-$startIND)*2)] .= $linearray[(2*$i)+5];
	    $genomat[(($i-$startIND)*2+1)] .= $linearray[(2*$i+1)+5];
	}

	$snpcount=$snpcount+1;
    }

                                ## print out:	
    if ($a==0)
    {
	if($v1){
	    print OUT "$totalINDS\n";
	}else {
	    print OUT "$totalhaps\n";
	}
	print OUT "$nsites\n";
	print OUT "P @posvec\n";
	if($v1){
	    for ($j=0; $j < $nsites; $j+=1)
	    {
		print OUT "S";
	    }
	    print OUT "\n";
	}
    }
    for ($i=0; $i < (2*($endIND-$startIND)); $i+=1)
    {
	print OUT "$genomat[$i]\n";
    }
}

