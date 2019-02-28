#!/usr/bin/perl

$file = "qc2.ene";
$etol = 0.01; #kcal/mol

@lll = qx ! ls -d */ !;
for ($i=0;$i<=$#lll;$i++) { $lll[$i] =~ s/\/\n// }
@lll = sort {$a <=> $b} @lll;

foreach $i (@lll) {
	$ls = qx ! ls $i/$file !;
	if (!(-e "$i/$file") && $ls =~ /CONF/) { die('$file not available yet\n') }
	if (-e "$i/$file") {
	$tmp = qx ! cat $i/$file !;
	@dump = split/ +/,$tmp;
	$ee[$i] = "$dump[0]";
	if ($ee[$i] < $emin) {$emin = $ee[$i]; $imin = $i};
	} else {
	$ee[$i] = 0.;
	}
}

open(OUT,">MIN$imin");
close(OUT);


