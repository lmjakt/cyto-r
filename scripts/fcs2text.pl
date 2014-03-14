#!/usr/bin/perl -w

## reads the .xml file produced by Aria in order to define a mapping
## note this .xml file has just two lines. 
## run the 

($infile) = @ARGV;
open(IN, $infile) || die "unable to open $infile $!\n";

$line = <IN>;
$line = <IN>;

@items = split /tube name=\"/, $line;

for $i(1..$#items){
    if($items[$i] =~ /([^"]+)\"\>/){
	$tube = $1;
	if($items[$i] =~ /<data_filename>([^<]+)</){
	    $file = $1;
	    system("read_fcs $file $tube")
	}
	print "$tube\t$file\n";
    }
}

