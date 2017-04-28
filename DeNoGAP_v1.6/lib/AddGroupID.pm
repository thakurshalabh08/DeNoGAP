##### Module to Perform MCL CLUSTERING#####
##### Author: Shalabh Thakur ################
##### Date: 6-AUG-2013 #####################

#!/usr/bin/perl
package AddGroupID;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(add_id);


sub add_id{

    my($mcl_output)=(shift);
    my($cluster_num)=(shift);
    my($hmm_program)=(shift);
    my(%homologue_group)=%{(shift)};
    my($cluster_out_dir)=(shift);
    my($mcl_dir)=(shift);

    my(%group_for_query)=();

    my $group_file="$cluster_out_dir/Cluster_$cluster_num.txt";

    open(MCL_OUT,"$mcl_dir/$mcl_output");

    open(GROUP_OUT,">$group_file");
   
    if($hmm_program eq "phmmer"){

        my @timedata = localtime();
        $timedata[5] += 1900;
        my $id_counter=100000;
        #my $id_counter=$timedata[0].$timedata[1].$timedata[2].$timedata[3].$timedata[4].$timedata[5];
        
        my $prefix="Group";

        foreach my $group_line(<MCL_OUT>){
               chomp($group_line);
               my $group_id=$prefix.$id_counter;
               $homologue_group{$group_id}=$group_line;
               $group_for_query{$group_id}=$group_line;

               $id_counter++;
        }

    }elsif($hmm_program eq "hmmscan"){

        my $prefix="Group";
        my @timedata = localtime();
        $timedata[5] += 1900;
        my $id_counter=100000;
        #my $id_counter=$timedata[0].$timedata[1].$timedata[2].$timedata[3].$timedata[4].$timedata[5];

        foreach my $group_line(<MCL_OUT>){

                chomp($group_line);

                if($group_line=~/(Group\d+)/){                  
                    $group_line=~/(Group\d+)/;
                    my $group_id=$1;                    
                    $group_line=~s/(Group\d+)(\t*)//g;
                    $group_line=~s/\t$//g;
                    
                    if(!defined($homologue_group{$group_id})){
                       print "Cannot find protein id's for cluster $group_id, Error in AddGroupID.pm script at line 70\n";
                       exit;
                    }
                    $group_line=$homologue_group{$group_id}."\t".$group_line;
                    $homologue_group{$group_id}=$group_line;
                    $group_for_query{$group_id}=$group_line;
                 }
                 else{   

                    GET_NEW_GROUP_ID:

                    my $group_id=$prefix.$id_counter; 

                    if($homologue_group{$group_id}){
                       $id_counter++;
                       goto GET_NEW_GROUP_ID;
                    }

                    $group_line=~s/\t$//g;
                    $homologue_group{$group_id}=$group_line;
                    $group_for_query{$group_id}=$group_line;
 
                    $id_counter++;
                 }
        }      

   }

   ####### Print Cluster Group #######

   foreach(sort{$a cmp $b}keys %homologue_group){
      print GROUP_OUT "$_:\t$homologue_group{$_}\n";
   }
   
   close GROUP_OUT;
   close MCL_OUT;  
   undef %homologue_group;
   
   return($group_file,\%group_for_query);
}
