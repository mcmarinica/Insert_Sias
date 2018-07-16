#!/usr/bin/perl
#The following will read one file and will read first line
$file=$ARGV[0];
$aOLD=$ARGV[1];
$aNEW=$ARGV[2];

open(INFILE,"<$file") || die "cant open \@ $file";
$i=0;
while($line=<INFILE>){
    if ( $i == 0 ) {
    @seqs = split / /, $line;
    @nv = &remove_space(@seqs);
    $n1=$nv[0];
    $n2=$nv[1];
    $n3=$nv[2];
    } 
    if ( $i == 1 ) {
    @seqs = split / /, $line;
    @nv = &remove_space(@seqs);
    $a11=$nv[0]*$aNEW/$aOLD;
    $a12=$nv[1]*$aNEW/$aOLD;
    $a13=$nv[2]*$aNEW/$aOLD;
    }; 
    if ( $i == 2 ) {
    @seqs = split / /, $line;
    @nv = &remove_space(@seqs);
    $a21=$nv[0]*$aNEW/$aOLD;
    $a22=$nv[1]*$aNEW/$aOLD;
    $a23=$nv[2]*$aNEW/$aOLD;
    };
    if ( $i == 3 ) {
    @seqs = split / /, $line;
    @nv = &remove_space(@seqs);
    $a31=$nv[0]*$aNEW/$aOLD;
    $a32=$nv[1]*$aNEW/$aOLD;
    $a33=$nv[2]*$aNEW/$aOLD;
    };     
    if ( $i == 4 ) {
    @seqs = split / /, $line;
    @nv = &remove_space(@seqs);
    $nat=$nv[0];
    } 
$i++;
};
close(INFILE); 


 $n1=$n1+0;
 $n2=$n2+0;
 $n3=$n3+0;
 $nat=$nat+0;
 printf "  %7i          %7i        %7i\n",$n1,$n2,$n3;
 printf"    %18.15f     %18.15f        %18.15f\n",$a11,$a12,$a13;
 printf"    %18.15f     %18.15f        %18.15f\n",$a21,$a22,$a23;
 printf"    %18.15f     %18.15f        %18.15f\n",$a31,$a32,$a33;
 printf "$nat\n"; 

open(INFILE,"<$file") || die "cant open \@ $file";
$i=0;
while($line=<INFILE>){
    if ( $i >= 5 ) {
    @seqs = split / /, $line;
    @nv = &remove_space(@seqs); 
   
    $x=$nv[0];
    $y=$nv[1];
    $z=$nv[2];
    $nn=$nv[3];
    printf "  %15.9f  %15.9f  %15.9f   %i\n",$x,$y,$z,$nn; 
    }
$i++;
};
close(INFILE); 

sub remove_space

{
   local($ii, $lo, @newvec,@tmpvec,@vect);

    @vect=@_ ;
    $lo = @vect;
    $ii=0;
    $ij=0;
     foreach (@vect) {
           $tmpvec[$ii]='v'.$vect[$ii].'v'; 
           if ($tmpvec[$ii] ne 'vv'){    
	    @newvec[$ij]=$vect[$ii];
	    $ij++;
	   } 
	   $ii++;
     }
     return @newvec;
}   
