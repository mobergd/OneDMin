#!/usr/bin/perl

$file = "qc1.ene";
$etol = 0.01; #kcal/mol

@lll = qx ! ls -d */ !;
for ($i=0;$i<=$#lll;$i++) { $lll[$i] =~ s/\/\n// }
@lll = sort {$a <=> $b} @lll;

foreach $i (@lll) {
	if (-e "$i/$file") {
	$tmp = qx ! cat $i/$file !;
	@dump = split/ +/,$tmp;
	$ee[$i] = "$dump[0]";
	if ($ee[$i] < $emin) {$emin = $ee[$i]};
	} else {
	$ee[$i] = 0.;
	}
}
foreach $i (@lll) {
	qx ! /bin/rm $i/FAIL $i/CONF* $i/DUPE* >& /dev/null !;
	$dupecheck = "";
	if ($ee[$i] == 0.) {
		$dupecheck = "FAIL";
	} else {
		$e=($ee[$i]-$emin)*627.509;
		foreach $j (@unique) { $tmp = abs($e-$eu[$j]); if ( abs($e-$eu[$j]) < $etol ) { $dupecheck = "DUPE$j"; last; } }
		if ( $dupecheck eq "") {
			@unique = (@unique, $i);
			$eu[$i]=$e;
			$dupecheck = "CONF".($#unique+1);
		}
	}
	open(OUT,">$i/$dupecheck");
	close(OUT);
}
