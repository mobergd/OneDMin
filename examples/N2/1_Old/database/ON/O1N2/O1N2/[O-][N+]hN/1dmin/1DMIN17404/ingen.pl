#!/usr/bin/perl

$suf="a";

$cpupernode=20;

$r1 = $ARGV[0];
$r2 = $ARGV[1];
$"="";
$pre = "sub";
$exe = "auto1dmin.x";

for ($i=$r1;$i<=$r2;$i++) { @dirs=(@dirs,"$suf$i"); }

open(SUB,">$pre-$r1.x");
print SUB "cd $dirs[0]\n";
foreach $dd (@dirs) {
	$count++;
	if ($count == 1) {
		$job++;
		open(SUB,">$pre-$job.x");
		print SUB "cd $dirs[0]\n";
	}
	if (-e $dd) {die("directory $dd already exists!\n")};
	mkdir $dd;
	qx ! cp r0/$exe $dd/. !;
	qx ! cp r0/smallest.geo $dd/. !;
	qx ! cp r0/qc.mol $dd/. !;
	qx ! cp r0/m.x $dd/. !;
	open(IN,"<r0/input");
	open(OUT,">$dd/input");
	$ran = int(rand(1000000000));
	foreach $line (<IN>) {
		$line =~ s/RANSEED/$ran/;
		print OUT $line;
	}
	close(IN);
	close(OUT);
	print SUB "cd ../$dd\n";
	print SUB "time ./$exe < input >  output &\n";
	if ($count == $cpupernode) {
		$count=0;
		print SUB "wait\n";
		close(SUB);
		qx ! chmod u+x $pre-$job.x !;
	}
}
qx ! chmod u+x $pre-$job.x !;
