use warnings;
use strict;


open(OF, ">../Data/HeatProc_Learn.txt");
print OF "patRIL","\t","chamber","\t","Pre","\t","Learning","\t","Memory","\t","file","\n";

open(AC, ">../Data/HeatProc_Activity.txt");
print AC "patRIL","\t","chamber","\t","Pre","\t","Learning","\t","Memory","\n";


my $directory = '/home/pwilliams/DSPR/RawData/Main_Heat_Calc_2';

    opendir (DIR, $directory) or die ("Can't open!");
    

while (my $file = readdir(DIR)) {

next if ($file =~ m/^\./);

print "$directory/$file\n"; 
my @rilset = split(/_/, $file);
my $rils = $rilset[1];
my $ril = substr $rils, 0, 5;
#print "$ril\n";
    
open (HC, "$directory/$file")  || die ("Can't open!");

my $counter = 0;

while (my $line = <HC>){
	chomp $line;
	if(($line =~ /^[0]/ && ($counter== 0 || $counter==1))|| $line =~ /why/){
	$counter = $counter+1;
	}
	
	if($line =~ /^[0-9]/ && $counter==1)
	{
	my @lineset = split /\s+/,$line;
	my $len = scalar @lineset;
	my $mem = $lineset[($len-1)];
	my $learn = $lineset[($len-2)];
	my $pre = $lineset[($len-3)];
	my $chamb = $lineset[($len-4)];
	my $gn = $lineset[($len-5)];
	#my($nn,$ril,$gp,$gn,$chamb,$pre,$learn,$mem) = split /\s+/,$line;
	#my($groupnum, $blank) = split/\./,$gn;
	print $gn, "\n";
	print OF $ril, "\t",$chamb, "\t",$pre, "\t",$learn, "\t",$mem,"\t",$gn, "\n";
	}
	
	if($line =~ /^[0-9]/ && $counter==2)
	{
  my @lineset = split /\s+/,$line;
	my $len = scalar @lineset;
	my $mem = $lineset[($len-1)];
	my $learn = $lineset[($len-2)];
	my $pre = $lineset[($len-3)];
	my $chamb = $lineset[($len-4)];
		
	#my($nn,$ril,$gp,$gn,$chamb,$pre,$learn,$mem) = split /\s+/,$line;
	#my($groupnum, $blank) = split/\./,$gn;
	print AC $ril, "\t",$chamb, "\t",$pre, "\t",$learn, "\t",$mem, "\n";
	#print $line, "\n";
	}
	
}

close (HC);
}

closedir(DIR);
close(OF);
close(AC);
