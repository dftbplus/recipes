#!/usr/bin/env perl

use Math::Trig;

# DFTB+ supercell and helical gen format repeat tool
# repeat a periodic .gen structure a specified number of times along
# lattice vectors. Input consists of gen file with an
# additional last line consisting of the number of additional repeats
# in each lattice vector direction, for example:
# 2 F
# Si
# 1 1 -0.125 -0.125 -0.125
# 2 1 0.125 0.125 0.125
# 0.0 0.0 0.0
# 2.71 2.71 0.0
# 2.71 0.0 2.71
# 0.0 2.71 2.71
# 2 3 4
# <- so repeat cell to give a 2x3x4 supercell
#
# or
# 2  H
# C
# 1 1    0.2756230044E+01    0.2849950460E+01    0.1794011798E+01
# 2 1    0.2656226397E+01    0.2949964389E+01    0.3569110265E+00
# 0 0 0
# 0.2140932670E+01 18.0 10
# 2 10
# <- so repeat cell to give a 2x length and full C_10 ring

if (!@ARGV) {@ARGV='-';
             $_=<>;
         }
else {$_=<ARGV>}

#get first line of gen file
s/^\s*(\d+)\s*([SsCcFfHh])//;
$atoms=$1;
if (!$atoms) {die"not a recognised gen file\n"}

#check if it's a fractional coords supercell file
if ($2 =~ /[Ff]/) {$frac=1} else {$frac=0}
if ($2 =~ /[Hh]/) {$helix=1} else {$helix=0}
# stop if a cluster file
if ($2 =~ /[Cc]/) {die "cluster file -- no point in repeating!\n"}

#get atom types from the comment file
($_=" ".<>)=~s/[^A-Za-z]+/ /g;
$types = $_;

#get $atoms data lines
for ($n=0;$n<$atoms;$n++) {
    $_=<>;
    s/^\s*\d+\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)//;
    if ($4 eq "") {die"there are only ",$n," atoms in the data, not $atoms!\n $1\n $2\n $3\n"}
    $atom[$n]=$1;$x[$n]=$2;$y[$n]=$3;$z[$n]=$4;
    s/.*//;
}

if (!$helix) {
    $_=<>; # origin line
    for ($i=0;$i<3;$i++) {
        $_=<>;
        s/^\s*(\S+)\s+(\S+)\s+(\S+)\s*[\n\r\f]// ||
            next;
        push(@vectsx,$1);push(@vectsy,$2);push(@vectsz,$3);
    }
    if (@vectsx!=3) {die "lattice vector data corrupted\n"}
# get repeat values
    $_=<>;
    s/^\s*(\d+)\s+(\d+)\s+(\d+)\s*[\n\r\f]//;
    $rep[0] = $1;$rep[1] = $2;$rep[2] = $3;
    for ($i=0;$i<3;$i++) {
        if ($rep[$i]<1 && int($rep[$i])!=$rep[$i]) {die "nonsensical repeat value $i: $rep[$i]\n"}
    }
    $rep[0]--;$rep[1]--;$rep[2]--;

# repeat each coord the specified number of times
    if ($frac) {
        for ($arep=0;$arep<=$rep[0];$arep++) {
            for ($brep=0;$brep<=$rep[1];$brep++) {
                for ($crep=0;$crep<=$rep[2];$crep++) {
                    $shiftAt[0] = $arep;
                    $shiftAt[1] = $brep;
                    $shiftAt[2] = $crep;
                    if ($arep==0 && $brep==0 && $crep==0) {next}
                    for ($n=0;$n<$atoms;$n++) {
                        push(@atom,$atom[$n]);
                        push(@x,$x[$n]+$shiftAt[0]);
                        push(@y,$y[$n]+$shiftAt[1]);
                        push(@z,$z[$n]+$shiftAt[2]);
                    }
                }
            }
        }
# increase number of atoms to new total
        $atoms=@atom;
# re-normalize the coordinates to fit inside new cell in fractional coords
# between 0 and 1
        for ($n=0;$n<$atoms;$n++) {
            $x[$n] /= (1+$rep[0]);
            $y[$n] /= (1+$rep[1]);
            $z[$n] /= (1+$rep[2]);
        }
    } else {
        for ($arep=0;$arep<=$rep[0];$arep++) {
            for ($brep=0;$brep<=$rep[1];$brep++) {
                for ($crep=0;$crep<=$rep[2];$crep++) {
                    $shiftAt[0] = $vectsx[0]*$arep + $vectsx[1]*$brep + $vectsx[2]*$crep;
                    $shiftAt[1] = $vectsy[0]*$arep + $vectsy[1]*$brep + $vectsy[2]*$crep;
                    $shiftAt[2] = $vectsz[0]*$arep + $vectsz[1]*$brep + $vectsz[2]*$crep;
                    if ($arep==0 && $brep==0 && $crep==0) {next}
                    for ($n=0;$n<$atoms;$n++) {
                        push(@atom,$atom[$n]);
                        push(@x,$x[$n]+$shiftAt[0]);
                        push(@y,$y[$n]+$shiftAt[1]);
                        push(@z,$z[$n]+$shiftAt[2]);
                    }
                }
            }
        }
# increase number of atoms to new total
        $atoms=@atom;
    }

# scale cell vectors
    for ($i=0;$i<3;$i++) {
        $vectsx[$i] *= ($rep[$i]+1);
        $vectsy[$i] *= ($rep[$i]+1);
        $vectsz[$i] *= ($rep[$i]+1);
    }

    if ($frac) {
        print "$atoms F\n";
    } else {
        print "$atoms S\n";
    }
} elsif ($helix) {
    $_=<>; # origin line
    $_=<>;
    /^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/ || die "Malformed helical lattice line!\n";
    $transA = $1; $angA = $2; $angRadA =  $angA * pi / 180.0;
    $transB = 0.0; $angB = $3; $angRadB =  2.0*pi / $angB;
    $_=<>;
    /^\s*(\d+)\s+(\d+)\s*[\n\r\f]/ || die "Malformed repeat line!\n";
    # translation, then rotation
    $rep[0] = $1;$rep[1] = $2;
    $rep[0]--;$rep[1]--;

    if ($rep[0]<0 && int($rep[0])!=$rep[0]) {die "nonsensical first repeat value!\n"}
    if ($rep[1]<0 && int($rep[1])!=$rep[1]) {die "nonsensical second repeat value!\n"}
    if ($angB % (1+$rep[1]) != 0) {die "Non commensurate 2nd repeat!\n"}
    if ((1+$rep[1]) > $angB) {die "Too many 2nd repeats!\n"}
    for ($arep=1;$arep<=$rep[0];$arep++) {
        $shiftAt[0] = 0.0;
        $shiftAt[1] = 0.0;
        $shiftAt[2] = $arep * $transA;
        $atmp = $angRadA * $arep;
        for ($n=0;$n<$atoms;$n++) {
            push(@atom,$atom[$n]);
            push(@x,($x[$n]+$shiftAt[0])*cos($atmp)-($y[$n]+$shiftAt[1])*sin($atmp));
            push(@y,($y[$n]+$shiftAt[1])*cos($atmp)+($x[$n]+$shiftAt[0])*sin($atmp));
            push(@z,$z[$n]+$shiftAt[2]);
        }
    }
    $atoms=@atom;
    for ($arep=1;$arep<=$rep[1];$arep++) {
        $atmp = $angRadB * $arep;
        for ($n=0;$n<$atoms;$n++) {
            if (($x[$n]**2 + $y[$n]**2) > 1.0E-10) { # avoid repeating
                # axial atoms, tollerance used is tolSameDist2 in
                # DFTB+ main code
                push(@atom,$atom[$n]);
                push(@x,($x[$n])*cos($atmp)-($y[$n])*sin($atmp));
                push(@y,($y[$n])*cos($atmp)+($x[$n])*sin($atmp));
                push(@z,$z[$n]);
            }
        }
    }
    $atoms=@atom;
    $angA *= ($rep[0]+1.0); $angA %= 360;
    $transA *= ($rep[0]+1.0);
    $angB /= ($rep[1]+1.0);
    print "$atoms H\n";
}


print "$types\n";
for ($n=0;$n<$atoms;$n++) {
    $i = $n+1;
    printf("%i %i %12.8f %12.8f %12.8f\n",$i,$atom[$n],$x[$n],$y[$n],$z[$n]);
}
if (!$helix) {
    printf("%12.8f %12.8f %12.8f\n",0,0,0);
    for ($i=0;$i<3;$i++) {
        printf("%12.8f %12.8f %12.8f\n",$vectsx[$i],$vectsy[$i],$vectsz[$i]);
    }
} else {
    printf("%12.8f %12.8f %12.8f\n",$transA,$angA,$angB);
}
