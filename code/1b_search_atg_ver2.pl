#!/usr/bin/perl
use strict;

# Usage:
# perl ./search_atg_ver2.pl [path to 5' fasta file] [output file]

# Change log:
# Version 2 - de Ruijter:
#   Changed the code so that it also runs on files where sequences are not seperated by an empty line.
# Version 1 - de Ruijter: 
#   Version 1 doesn't exist (I went from 0 to 2 to avoid confusion).
# Version 0 - Oda:
#   Initial version

############################
# find ATG and STOP codons #
############################
#my $file = "fbp1.fasta";
my $file = $ARGV[0];
my $out_file = $ARGV[1];
open(IN, "<$file") or die;
 open(OUT, ">$out_file") or die;

my($x,$y,$len,$frame,$id,$num)=("","",0,0,"",0);
my @atg=();my @stop=();my @base=();
#my ($count,$count2) =(0,0);
### OUTPUT ###
print OUT "GeneName\t5UTR_Length\tFrame\tATG\tATG-FRAME(j)\tATG-pos\tSTOP\tSTOP-FRAME(j)\tSTOP-pos\tuORF_len\n";

#Count rows.###############
my    $rows = 0;
open(FILE, $file) or die "Can't open `$file': $!";
while (my $a =<FILE>) {
  $rows ++;
}
close FILE;
# print "$rows\n";
my $row_n=0;
#Yomikomi###################
#do {my $line=<IN>;
while(my $line=<IN>){
  chomp $line;
  if($line =~/^>/ ){
    if($row_n > 0){
		@base=split//,$x;
		$len=length($x);
		for(my$k=0;$k<$len;$k++){
			if($base[$k]=~/[^A|T|G|C|a|t|g|c]/){ #CHECK FOR NON NUCLEOTIDE.
				die "BASE_ERROR";
			}
		}
		#  print "$len,\t";
		for(my$i=0;$i<3;$i++){ # "$i" represents the frame to be checked.
			for(my$j=0;$j>-$len/3;$j--){ # the "$j"_th frame
				#check codon
				#	my $codon=$base[($len-3*$j-$i-2)].$base[($len-3*$j-$i-1)].$base[($len-3*$j-$i)];
				my $codon=$base[($len+3*$j-$i)].$base[($len+3*$j-$i+1)].$base[($len+3*$j-$i+2)];
				#	print "$i,$j,";
				#	print "$codon\n";
				#
				if($codon =~/ATG/){
					push (@atg, $j, $codon);
				}
				elsif($codon=~/TAA/ or $codon=~/TAG/ or $codon=~/TGA/){
					push (@stop, $j, $codon);
				  #
				}
			}  
			#Compare ATG-pos and STOP-pos
			#      $"="\t";
			#	  print "i=","$i\t","atg(J/I/codon)", "@atg","\n";
		  	my$atg_count=  @atg;
			#	print "$atg_count,@atg\n";
			my$stop_count=  @stop;
			for(my$n=0;$n<$atg_count/2;$n++){
				my $atg_pos=3*$atg[2*$n]-$i;
				for(my$m=0;$m<$stop_count/2;$m++){
					if(@stop[2*$m]>@atg[2*$n]){
						my $stop_pos=2+3*$stop[2*$m]-$i;
						my $uorf=(-$atg_pos+$stop_pos-2)/3;
						##out codon
						print OUT "$id\t$len\t-$i\t$atg[2*$n+1]\t$atg[2*$n]\t$atg_pos\t$stop[2*$m+1]\t$stop[2*$m]\t$stop_pos\t$uorf\n";
						#	  print "i=","$i\t","stop(J/I/codon)", "@stop","\n";
					}
				}
			}
			@atg=();@stop=();
			$x="";
		}
    }
    $id=$line;
    $row_n++;
#    print "$line,";
  }elsif($line=~/A|T|G|C|a|t|g|c/){
    $x=$x.$line;
    $row_n++;
  }elsif($line=~/^$/){
    $row_n++;
  }else{
    print "INPUT_FILE_FORMAT_ERROR\n";
  }
#  print "$row_n\n";
}#while($row_n <= $rows);
# print "$.\n";
############################
