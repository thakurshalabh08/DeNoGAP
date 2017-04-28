##### Module to Cluster Ortholog Pairs #######
##### Author: Shalabh Thakur #################
##### Date: 28-SEP-2013 ######################

#!/usr/bin/perl
package OrthologCluster;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Cluster;
use SQLiteDB;
use Hash::Merge qw(merge );
use File::Path qw(remove_tree);
use Tie::File;

my $cluster_id=(shift);
my($ortholog_pair_dir)=(shift);
my($inparalog_pair_dir)=(shift);
my($cluster_dir)=(shift);
my($db_dir)=(shift);
my($mcl_inflation_value)=(shift);

        print "$cluster_id\n";
   
        #### CREATE TEMP DATABASE ######

        chdir($db_dir);   
    
        my $dbh=DBI->connect("dbi:SQLite:dbname=$cluster_id","","") or die $DBI::errstr;
        my $stmt=$dbh->prepare("SELECT SQLITE_VERSION()");
        $stmt->execute();   
        $stmt->finish();

        $dbh->do("CREATE TABLE IF NOT EXISTS $cluster_id(taxaA TEXT, idA TEXT , taxaB TEXT, idB TEXT, divergence REAL, Cluster TEXT)");

        chdir($Bin);

        ##### READ ORTHOLOG PAIR ###
        my %list_protein_id=();
        
        open(ORTHOLOG_FILE,"$ortholog_pair_dir/$cluster_id.txt");
        open(INPARALOG_FILE,"$inparalog_pair_dir/$cluster_id.txt");
        
        my @ortholog_pair=map{s/\|/	/g; $_;}<ORTHOLOG_FILE>;
        my @inparalog_pair=map{s/\|/	/g; $_;}<INPARALOG_FILE>;

        my @inparalog_pair_2=();

        foreach my $inpara_line(@inparalog_pair){

           my @inpara=split(" ",$inpara_line);

           splice(@inpara,5,1);

           push(@inparalog_pair_2,join("\t",@inpara));
        }

        my $sql_stmt="INSERT into $cluster_id (taxaA,idA,taxaB,idB,divergence,Cluster) VALUES (?,?,?,?,?,?)";

        SQLiteDB::load_from_array($db_dir,$cluster_id,\@ortholog_pair,$sql_stmt);
        SQLiteDB::load_from_array($db_dir,$cluster_id,\@inparalog_pair_2,$sql_stmt);
   
        my $pair_stmt="Select distinct(taxaA),idA,taxaB,idB,divergence from $cluster_id";
      
        my($get_all_pairs)=SQLiteDB::get_record($db_dir,$cluster_id,$pair_stmt);

        open(MCL_WEIGHT_TABLE,">$db_dir/$cluster_id.abc");
        
       if(scalar(@{$get_all_pairs})>=1){

            foreach my $row(@{$get_all_pairs}){

                my $taxaA=shift(@{$row});
                my $idA=shift(@{$row});
                my $taxaB=shift(@{$row});
                my $idB=shift(@{$row});
                my $divergence=shift(@{$row});

                my $lableA=$taxaA."|".$idA;
                my $lableB=$taxaB."|".$idB;
                my $weight=1-$divergence;
                   
                print MCL_WEIGHT_TABLE "$lableA\t$lableB\t$weight\n";
            }
       }  

       close MCL_WEIGHT_TABLE;

      system("mcxload -abc $db_dir/$cluster_id.abc --stream-mirror -ri max --stream-neg-log10 -o $db_dir/$cluster_id.mci -write-tab $db_dir/$cluster_id.tab --write-binary");

      system("mcl $db_dir/$cluster_id.mci -I $mcl_inflation_value -use-tab $db_dir/$cluster_id.tab -o $cluster_dir/orthoCluster_$cluster_id.txt");
    
      my @cluster=();

      tie @cluster, 'Tie::File', "$cluster_dir/orthoCluster_$cluster_id.txt";

      my $clust_count=1;

      foreach my $cluster_line(@cluster){

          my $ortho_id=$cluster_id.".".$clust_count;

          $cluster_line=$ortho_id."\t".$cluster_line;

          $clust_count++;       
      } 

      system("rm $db_dir/$cluster_id");
      system("rm $db_dir/$cluster_id.abc");
      system("rm $db_dir/$cluster_id.mci");
      system("rm $db_dir/$cluster_id.tab");

