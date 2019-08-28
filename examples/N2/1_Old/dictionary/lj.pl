#!/usr/bin/perl

open(IN,"<$ARGV[0]");

$vcutlo = -20000;
$vcuthi = 200000;

$rhi=-1000;
$shi=-1000;
$vhi=-1000;
$rlo=1000;
$slo=1000;
$vlo=1000;

foreach $line (<IN>) {
chomp($line);
@dat=split/ +/,$line;

if ($dat[8] != 0 && $dat[1] =~ /[0-9]/) {

#print "$dat[6]\n";

$nn++;
if ($dat[6] > $vcutlo && $dat[6] < $vcuthi) {
$r+=$dat[4];
$v+=$dat[6];
$s+=$dat[7];
$r2+=($dat[4])**2;
$v2+=($dat[6])**2;
$s2+=($dat[7])**2;
$n++;

if ($dat[4] < $rlo) {$rlo = $dat[4]};
if ($dat[4] > $rhi) {$rhi = $dat[4]};
if ($dat[6] < $vlo) {$vlo = $dat[6]};
if ($dat[6] > $vhi) {$vhi = $dat[6]};
if ($dat[7] < $slo) {$slo = $dat[7]};
if ($dat[7] > $shi) {$shi = $dat[7]};

}
}

}

$r/=$n;
$v/=$n;
$s/=$n;
$r2/=$n;
$v2/=$n;
$s2/=$n;
$sr=sqrt($r2-$r**2);
$sv=sqrt($v2-$v**2);
$ss=sqrt($s2-$s**2);

#print "$n hits of $nn total\n";
#print "\n";
#print "       lo      hi\n";
#print " v   $vlo    $vhi\n";
#print " r   $rlo    $rhi\n";
#print " s   $slo    $shi\n";
#print "\n";
#printf("% 15.7f",$r);
#printf("% 15.7f",$sr);
#printf("% 15.7f",$v);
#printf("% 15.7f",$sv);
#printf("% 15.7f",($r/2**(1/6)));
#printf("% 15.7f",($sr/2**(1/6)));
#print "\n";
#printf("% 15.7f",0);
#printf("% 15.7f",0);
#printf("% 15.7f",$v);
#printf("% 15.7f",$sv);
#printf("% 15.7f",$s);
#printf("% 15.7f",$ss);
#print "\n";

#print "\n";
#print " $r $v \n";
#if ($s != $r) {#print " $s 0. \n";}
$v=-$v;
print " CollisionFrequency\n   LennardJones\n     Epsilons[1/cm]        $v   $v \n     Sigmas[angstrom]         $s  $s  \n     Masses[amu]             0. 0. \n   End\n";

