#!perl 

#v2 imposes a mising data threshold -- max = 10% missing data
#PA

# count_derived_mutations.pl

# catalogs eight site categories
# array of site categories  my @site_categories =(0,0,0,0,0,0,0,0);
# 1 = fixed Dsim
# 2 = fixed Dmau
# 3 = fixed Dsim poly Dmau
# 4 = fixed Dmau poly Dsim
# 5 = poly Dmau and poly Dsim
# 6 = poly Dmau
# 7 = poly Dsim
# 8 = fixed Dsim and fixed Dmau

if (@ARGV < 1){print "\nusage: perl count_derived_mutations_ELS infile.fasta \nNote: input is a fasta multialignment:\n\n>Dmel_outgroup\nGACTGACTGACTGACTGACTGACT\n>Dsim\nGACTGACTGACTGACTGACTGACT\n>Dsim\nGACTGACTGACTGACTGACTGACT\n>Dmau\nGACTGACTGACTGACTGACTGACT\n>Dmau\nGACTGACTGACTGACTGACTGACT\n\n Note -- this is currently hard-coded to look for Dsim and Dmau sequences. ELS version produces an ELS_results file, which is a table of frequencies at variant sites in the two species\n\n";exit;}

my $infile = shift(@ARGV);
open IN, "$infile" or die "cannot open input fasta file -- sure you got the name/path right?\n";

system("rm ELS_results");  # remove previous results file

my %mau_sites_freqs =();      # Dmau sites and their derived allele frequences
my %sim_sites_freqs =();      # Dsim sites and their derived allele frequences

my @site_categories =(0,0,0,0,0,0,0,0);
# print "\n";

#************************ READ IN DATA HERE **************************************
my %seq_names = ();
my %sequences = ();
my $nsam_mau=0; my $nsam_sim=0;
while (my $line = <IN>) {
        chomp $line; if ($line =~ /mau/){$nsam_mau++;} else{if ($line =~ /sim/){$nsam_sim++;}}
        my $name = $line;
        $line = <IN>; chomp $line;
        $sequences{$name}=$line;
        #}
# close IN;  # changed to allow multiple simulated replicates
# print $name, "\n";
if (($name =~ /outgroup/) || ($name =~ /mel/)){
    
    # print "\nnumber of sequences: ", scalar(keys %sequences), "\n";
    # print "number Dmau: $nsam_mau\n";
    # print "number Dsim: $nsam_sim\n";

    my $seqlen = length((values %sequences)[0]);
    # print "seq length: $seqlen\n";

    #open LIST, ">sequence_list"; foreach $key (sort keys %sequences) {print LIST "$key\n"; # print "$key\n"; # print LIST "%sequences{$key}\n";
    #}
    # print "\n";
    
    open ELS, ">>ELS_results";
    print ELS "position\tfreq_Dmau\tfreq_Dsim\n";
    
    for ($i=0; $i < $seqlen; $i++){  # for each site
        
      my $nsam_mau_mod = $nsam_mau; my $nsam_sim_mod = $nsam_sim; # reset samples sizes
        
      # check outgroup and ignore site if outgroup = N
      my $ignore=0;
      foreach $key (sort keys %sequences) {
          if (($key =~ /outgroup/) || ($key =~ /mel/)){
              if (substr($sequences{$key},$i,1) eq "N"){$ignore=1;}
              } # if outgroup set to ignore
           } # foreach sequence
      if ($ignore ==0){
        my %mau_states_freqs=();
        my %sim_states_freqs=();
        my %anc_states =();
          
        # go through all sequences and only process alleles that are not N
        foreach $key (sort keys %sequences) {
            # print "$key\t", substr($sequences{$key},$i,1), "\n";
            if (substr($sequences{$key},$i,1) eq "N"){ # if state == N substract 1 from sample sizes
                    if ($key =~ /mau/){$nsam_mau_mod=$nsam_mau_mod-1;}
                    if ($key =~ /sim/){$nsam_sim_mod=$nsam_sim_mod-1;}
                    # print "ignoring site ", $i+1," in $key\n";
            }
            else{
                if (($key =~ /mau/) && (substr($sequences{$key},$i,1) ne "N")){$mau_states_freqs{substr($sequences{$key},$i,1)}++;}
                if (($key =~ /sim/) && (substr($sequences{$key},$i,1) ne "N")){$sim_states_freqs{substr($sequences{$key},$i,1)}++;}
                if (($key =~ /outgroup/) || ($key =~ /mel/)){$anc_states{substr($sequences{$key},$i,1)}=substr($sequences{$key},$i,1);}
                } # if state != N
            }# foreach sequence
    
        # check all is good with counting alleles
        # print "\n";
        # foreach $key (sort keys %mau_states_freqs){print "Dmau states: $key\t", $mau_states_freqs{$key}, "\n";}
        # foreach $key (sort keys %sim_states_freqs){print "Dsim states: $key\t", $sim_states_freqs{$key}, "\n";}
        # foreach $key (sort keys %anc_states){print "outgroup state: $key\t", $anc_states{$key}, "\n";}
        # print "\n";
          
    if (($nsam_mau_mod/$nsam_mau >0.9) && (($nsam_sim_mod/$nsam_sim >0.9))){ # 10% missing data max
          
        # hashes of derived allele frequencies (DAF)
        $mau_sites_freqs{$i}=0;
        foreach $key (sort keys %mau_states_freqs){
        if ($key ne $anc_states{$key}){
            $mau_sites_freqs{$i} = $mau_states_freqs{$key}/$nsam_mau_mod;
            } # if derived
        } # for each key
    
        $sim_sites_freqs{$i}=0;
        foreach $key (sort keys %sim_states_freqs){
            if ($key ne $anc_states{$key}){
                $sim_sites_freqs{$i} = $sim_states_freqs{$key}/$nsam_sim_mod;
                } # if derived
            } # for each key
        
        # if site is variant, print to results file
        if (($mau_sites_freqs{$i}>0) || ($sim_sites_freqs{$i}>0)){
        print ELS $i+1,"\t", $mau_sites_freqs{$i},"\t", $sim_sites_freqs{$i},"\n";
        }
    
        # count up categories of sites
        # print "site ", $i+1,"\tsim:",$sim_sites_freqs{$i}, "\tmau:",$mau_sites_freqs{$i}, "\n";
        if (($mau_sites_freqs{$i}==1) && ($sim_sites_freqs{$i}==1)){$site_categories[7]++;} # if both fixed ignore
        else{
            # print "not fixed in both species\n";
            if ($sim_sites_freqs{$i}==1) {                                  # sim is fixed
                # print "site ", $i+1," sim fixed\t";
                if ($mau_sites_freqs{$i}==0){$site_categories[0]++;
                    # print "type 1\n";                                           # site type 1
                    }
                else {$site_categories[2]++; #print "type 3\n";                   # site type 3
                    }
                }
            elsif ($mau_sites_freqs{$i}==1) {                                # mau is fixed
                # print "site ", $i+1," mau fixed\t";
                if ($sim_sites_freqs{$i}==0){$site_categories[1]++;
                    # print "type 2\n";                                            # site type 2
                }
                else {$site_categories[3]++;  # print "type 4\n";                  # site type 4
                    }
                }
            else{                                              # neither mau or sim are fixed
                # print "site $i not fixed in either species\n";
                if (($mau_sites_freqs{$i}>0) and ($sim_sites_freqs{$i}>0)){         # site type 5
                    $site_categories[4]++;
                    # print "site ", $i+1," poly both\tprint type 5\n";
                    } #both polymorphic site 5
                elsif (($mau_sites_freqs{$i}>0) and ($sim_sites_freqs{$i}==0)){     # site type 6
                    $site_categories[5]++;
                    # print "site ", $i+1," private mau poly\tprint type 6\n";
                    } #mau polymorphic site 6
                elsif (($mau_sites_freqs{$i}==0) and ($sim_sites_freqs{$i}>0)){     # site type 7
                    $site_categories[6]++;
                    # print "site ", $i+1," private sim poly\tprint type 7\n";
                    } #sim polymorphic site 7
                } # neither mau or sim are fixed
            } # if not fixed in both species
        # print "\n";
      } # if ignore ==0;
      # else {print "\nignoring site ", $i+1," because ref = N\n";}
      } # if both species have > threshold % sample sizes after removing Ns -- hardcoded as 90% for now
    } # for each site

    for (my $j=0; $j<(scalar @site_categories); $j++){
        print "$site_categories[$j]\t";
        }
    print "\n";
    
    # re-initialize hashes
    %mau_sites_freqs =();      # Dmau sites and their derived allele frequences
    %sim_sites_freqs =();      # Dsim sites and their derived allele frequences
    %seq_names = ();
    %sequences = ();
    $nsam_mau=0; $nsam_sim=0;
    @site_categories =(0,0,0,0,0,0,0,0);
    } # if name contains outgroup

} # while there are lines in the file
close IN;
# print "\n";
    
