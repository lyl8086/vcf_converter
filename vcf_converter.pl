#!/usr/bin/perl
use strict;
use warnings;
use List::MoreUtils qw{uniq};
use Data::Dumper;
die "Usage:$0 [vcf] [popmap] [format] [prefix]\n" unless @ARGV >= 3;
my $vcf_in = $ARGV[0];
my $pop    = $ARGV[1];
my $format = $ARGV[2];
my $name   = $ARGV[3];

my %genpop_num = (
    'A'=>100,
    'C'=>110,
    'G'=>120,
    'T'=>130,
    'N'=>000,
); 
my %str_num = (
    $genpop_num{'A'}=>1,
    $genpop_num{'C'}=>2,
    $genpop_num{'G'}=>3,
    $genpop_num{'T'}=>4,
    '000'=>'-9',
);
my %ref_alt = (
    '0'=>'refb',
    '1'=>'altb',
    'altb'=>'1',
    'refb'=>'0',
);

my %FMT = (
    'bayescan'=>1,
    'genepop'=>1,
    'str'=>1,
    'arp'=>1,
    'all'=>1,
    'bayenv'=>1,
    'bayesass'=>1,
);

die "  Unknown format '$format', accept ".join(",",keys %FMT).".\n" unless defined $FMT{$format};

my ($popmap,$order) = parse_popmap($pop);
my ($vcf, $tot)     = vcf_parser($vcf_in,$popmap,$order);

# print Dumper($vcf);

if ($format eq 'bayescan') {
    vcf_2_bayescan($vcf, $order);
} elsif ($format eq 'genepop') {
    vcf_2_genepop($vcf, $popmap, $order);
} elsif ($format eq 'str') {
    vcf_2_structure($vcf, $popmap, $order);
} elsif ($format eq 'arp') {
    vcf_2_arlequin($vcf, $popmap, $order);
} elsif ($format eq 'bayenv') {
    vcf_2_bayenv($vcf, $popmap, $order);
} elsif ($format eq 'bayesass') {
    vcf_2_bayesass($vcf, $popmap, $order);
} elsif ($format eq 'all') {
    
    open STDOUT, ">$name.bayescan";
    vcf_2_bayescan($vcf, $order);
    open STDOUT, ">$name.gen";
    vcf_2_genepop($vcf, $popmap, $order);
    open STDOUT, ">$name.str";
    vcf_2_structure($vcf, $popmap, $order);
    open STDOUT, ">$name.arp";
    vcf_2_arlequin($vcf, $popmap, $order);
    open STDOUT, ">$name.inp";
    vcf_2_bayesass($vcf, $popmap, $order);
} else {
    print STDERR "Unkonw foramt $format, accept 'bayescan' 'genepop' 'str' 'arp' 'all'\n";
}

sub vcf_parser {
    #
    # vcf = [ { 
    #           'id' =>{} 
    #           'ind'=> {
    #                           ind=>{
    #                           'geno'=> '',
    #                           'A'=>'', 'C'=>'', 'G'='', 'T'=>'',
    #                           'ref'=>'', 'alt'=>'',}
    #                         }
    #           'pop'=> { pop=>{'A'=>'', 'C'=>'', 'G'='', 'T'=>'',
    #                           'ref'=>'', 'alt'=>'',}
    #                         } 
    #          }
    #        ]
    my $vcf_f = shift;
    my $pops  = shift;
    my $order = shift;
    my $head  = 0;
    my $tot   = 0;
    my (@header, $fmt, %samples, $vcf);
    open(my $in_fh, $vcf_f) or die $!;
    while(<$in_fh>) {
        next if /^##|^$/;
        next if $head > 0 && /^#/; 
        chomp;
        if ($head == 0 && /^#/) {
            push @header, split(/\t/);
            for (my $i=9; $i<=$#header; $i++) {
                $samples{$header[$i]} = $i;
            }
            $head = 1;
            next;
        }
        my $vcfline = $_;
        $tot++;
        my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genos) = split(/\t/, $vcfline);
        my @formats = split(/:/, $format); 
        map { $fmt->{$formats[$_]} = $_;} 0..$#formats; # Geno => Order.
        my $dat->{'id'} = $chrom . '_' . $pos;
        foreach my $pop (@$order) {
            # each population.
            foreach my $ind (@{$pops->{$pop}}) {
                # each individual.
                my $rank = $samples{$ind};
                my $geno0 = $genos[$rank-9];
                my @geno = split(/:/, $geno0);
                my $GT   = $geno[$fmt->{'GT'}];
                my $cnt  = calc_stat($GT, $ref,$alt);
                $dat->{'ind'}->{$ind} = {
                            'gt'   => $cnt->{'gt'},
                            'geno' => $GT,
                            'ref'  => $cnt->{'ref'},
                            'alt'  => $cnt->{'alt'},
                            'A'    => $cnt->{'A'},
                            'C'    => $cnt->{'C'},
                            'G'    => $cnt->{'G'},
                            'T'    => $cnt->{'T'},
                            'refb' => $ref,
                            'altb' => $alt,
                        };
                $dat->{'pop'}->{$pop}->{'geno'} .= $GT;
                $dat->{'pop'}->{$pop}->{'ref'}  += $cnt->{'ref'};
                $dat->{'pop'}->{$pop}->{'alt'}  += $cnt->{'alt'};
                $dat->{'pop'}->{$pop}->{'A'}    += $cnt->{'A'};
                $dat->{'pop'}->{$pop}->{'C'}    += $cnt->{'C'};
                $dat->{'pop'}->{$pop}->{'G'}    += $cnt->{'G'};
                $dat->{'pop'}->{$pop}->{'T'}    += $cnt->{'T'};
            }    
        }
        push @$vcf, $dat;
    }
    return $vcf, $tot;
}

sub parse_popmap {
    my $pop = shift;
    my ($pops, @order);
    open(my $in_fh, $pop) or die "$!";
    while (<$in_fh>) {
            next if /^#|^$/;
            $_ =~ s/[\r\n]|^\s+|\s+$//g;
            my @part = split;
            push @order, $part[1];
            push @{$pops->{$part[1]}}, $part[0]; 
            
        }
    close $in_fh;
    @order = uniq @order;
    return $pops, \@order;
}

sub calc_stat {
    my $gt  = shift;
    my $ref = shift;
    my $alt = shift;
    my $cnt;
    # initialize.
    $cnt->{'ref'} = 0;
    $cnt->{'alt'} = 0;
    $cnt->{'A'}   = 0;
    $cnt->{'C'}   = 0;
    $cnt->{'G'}   = 0;
    $cnt->{'T'}   = 0;
    if ($gt eq '0/1' || $gt eq '1/0' || $gt eq '1|0' || $gt eq '0|1') {
        $cnt->{'ref'} = 1;
        $cnt->{'alt'} = 1;
        $cnt->{$ref}  = 1;
        $cnt->{$alt}  = 1;
    } elsif ($gt eq '0/0' || $gt eq '0|0') {
        $cnt->{'ref'} = 2;
        $cnt->{$ref}  = 2;
    } elsif ($gt eq '1/1' || $gt eq '1|1') {
        $cnt->{'alt'} = 2;
        $cnt->{$alt}  = 2;
    }
    $cnt->{'gt'} = gencode($ref, $alt, $cnt->{'ref'}, $cnt->{'alt'});
    return $cnt;
}

sub gencode {
    my $ref  = shift;
    my $alt  = shift;
    my $cnt1 = shift;
    my $cnt2 = shift;
    my $code;
    if (($cnt1+$cnt2) == 0) {
        $code = '000000';
    } elsif ($cnt1==1&&$cnt2==1) {
        $code = $genpop_num{$ref}.$genpop_num{$alt};
    } elsif ($cnt1==2) {
        $code = $genpop_num{$ref}.$genpop_num{$ref};
    } elsif ($cnt2==2) {
        $code = $genpop_num{$alt}.$genpop_num{$alt};
    } else {
        $code = '000000';
    }
    return $code;
}


sub vcf_2_structure {
    #
    my $vcf    = shift;
    my $pops   = shift;
    my $order  = shift;
    my $pop_id = 0;
    my @id;
    my $i = 0;
    
    map {++$i;push @id,"$i"} @$vcf;
    
    print "\t\t", join(" ", @id), "\n";
    
    foreach my $pop (@$order){
        $pop_id++;
        foreach my $ind (@{$pops->{$pop}}) {
            my @ref_gts = ();
            my @alt_gts = ();
            map {
                my $gt = $_->{'ind'}->{$ind}->{'gt'};
                push @ref_gts, $str_num{substr($gt,0,3)};
                push @alt_gts, $str_num{substr($gt,3,3)};
                } @$vcf;
            print "$ind\t$pop_id\t", join(" ", @ref_gts), "\n";
            print "$ind\t$pop_id\t", join(" ", @alt_gts), "\n";
        }
    }
 

}

sub vcf_2_bayescan {
    my $vcf    = shift;
    my $order  = shift;
    my $n_loci = scalar(@$vcf);
    my $n_pops = scalar(@$order);
    print "[loci]=$n_loci\n\n[populations]=$n_pops\n\n";
    my $i = 0;
    foreach my $pop (@$order) {
        my $j = 0;
        $i++;
        print "[pop]=$i\n";
        foreach my $id (@$vcf) {
            my $ref = $id->{'pop'}->{$pop}->{'ref'};
            my $alt = $id->{'pop'}->{$pop}->{'alt'};
            print join("\t", ++$j, $ref+$alt, 2, $ref, $alt), "\n";
        }
    }

}

sub vcf_2_genepop {
    #
    my $vcf    = shift;
    my $pops   = shift;
    my $order  = shift;
    my @id;
    my $i = 0;
    map {++$i;push @id,"$i"} @$vcf;
    
    print "vcf2genpop\n",
    join("\n", @id), "\n";
    foreach my $pop (@$order){
        print "Pop\n";
        foreach my $ind (@{$pops->{$pop}}) {
            my @gts = ();
            map {push @gts, $_->{'ind'}->{$ind}->{'gt'}} @$vcf;
            print "$ind ,  ", join(" ", @gts), "\n"; 
        }
    }
 
}

sub vcf_2_arlequin {
    my $vcf    = shift;
    my $pops   = shift;
    my $order  = shift;
    my $n_loci = scalar(@$vcf);
    my $n_pops = scalar(@$order);
    my $profile= "[Profile]
  Title = \"\"
  NbSamples = $n_pops
  DataType = DNA
  GenotypicData = 1
  LocusSeparator = WHITESPACE
  MissingData = \"?\"
  GameticPhase = 0
  RecessiveData = 0
  
[Data]\n";
    print $profile;
    
    my $max_len = 0;
    foreach my $pop (@$order) {
        foreach my $ind (@{$pops->{$pop}}) {
            my $len  = length($ind);
            $max_len = $len if $len > $max_len;
        }
    }
    my $len = $max_len + 9;
    
    foreach my $pop (@$order){
         
         my ($sample_dat,$n);
         foreach my $ind (@{$pops->{$pop}}) {
            $n++;
            my ($g1, $g2);
            map {
                my $g  = $_->{'ind'}->{$ind}->{'geno'};
                my @g  = split(/\//, $g);
                if ($g[0] eq '.') {
                    $g1 .= ' ?';
                    $g2 .= ' ?';
                } else {
                    $g1 .= ' ' . $_->{'ind'}->{$ind}->{$ref_alt{$g[0]}};
                    $g2 .= ' ' . $_->{'ind'}->{$ind}->{$ref_alt{$g[1]}};
                }
                
            } @$vcf;
            $sample_dat .= '    '.sprintf("%-${max_len}s",$ind).'  1  '.$g1."\n";
            $sample_dat .= ' 'x$len.$g2."\n";
        }
        print "  [[Samples]]\n";
        print "    SampleName = \"$pop\"\n";
        print "    SampleSize = $n\n";
        print "    SampleData = {\n$sample_dat";
        print "    }\n";
    }

}

sub vcf_2_bayenv {
    my $vcf    = shift;
    my $pops   = shift;
    my $order  = shift;
    my $n_loci = scalar(@$vcf);
    my $n_pops = scalar(@$order);
    foreach my $snp (@$vcf) {
        my @out1 = ();
        my @out2 = ();
        foreach my $pop (@$order){
            push @out1, $snp->{'pop'}->{$pop}->{'ref'};
            push @out2, $snp->{'pop'}->{$pop}->{'alt'};
        }
        print join("\t", @out1), "\n";
        print join("\t", @out2), "\n";
    }
    
}

sub vcf_2_bayesass {
    my $vcf    = shift;
    my $pops   = shift;
    my $order  = shift;
    my $n_loci = scalar(@$vcf);
    my $n_pops = scalar(@$order);
    my $i      = 0;
    foreach my $snp (@$vcf) {
        $i++;
        foreach my $pop (@$order){
            foreach my $ind (@{$pops->{$pop}}) {
                my $gt = $snp->{'ind'}->{$ind}->{'gt'};
                my ($g1,$g2) = (0,0);
                if ($gt ne '000000') {
                   ($g1,$g2) = (substr($gt,0,3),substr($gt,3,3));
                }
                print join("\t", $ind, $pop, $i, $g1, $g2), "\n";
            }
        }
    }
}
