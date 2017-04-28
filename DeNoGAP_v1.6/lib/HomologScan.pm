##### Module to Perform Homolog Prediction Analysis#####
##### Author: Shalabh Thakur ################
##### Date: 6-AUG-2013 #####################

#####Version 2.1: Removed CreateModel Module, the feature to create hmm model is merged into SelectModelGene ####

#!/usr/bin/perl
package HomologScan;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Hmmer;
use FilterPair;
use Getopt::Long;
use List::MoreUtils qw(uniq);
use List::Util qw(sum);
use File::Basename;
use File::Copy;
use File::Path qw(remove_tree);
use Parallel::ForkManager;
use Hash::Merge qw( merge );
use Tie::File;
use MclClustering;
use AddGroupID;
use Initialize;
use Cluster;
use PartialMap;
use SQLiteDB;

use vars qw(@ISA @EXPORT @EXPORT_OK $db_name $db_dir);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(run_compare_reference run_homologue_scan);


sub run_compare_reference {

    $db_dir=(shift);
    $db_name=(shift);
    my %config_param=%{(shift)};
    my %out_dir=%{(shift)};
    my $output_dir=(shift);
    my $tmp_dir=(shift);
    

    if($config_param{ACTIVATE_ANALYSIS}->{COMPARE_REFERENCE}=~/YES/i){          
 
          my %homologue_group=();
          my %group_gene=();
          my $start_reference_scan = time;
          my $hmm_program="phmmer";
          my $cpu_core=$config_param{PARAMETERS}->{PARALLEL_CPU_CORE};
          my $hmmer_opt=$config_param{HMMER_PARAMETERS}->{HMMER_OPT};
          my $mcl_inflation_value=$config_param{PARAMETERS}->{MCL_INFLATION_VALUE};
          my $similarity_file="SimilarityPair.txt";
          my $decoy_cluster_file="$tmp_dir/decoy_file.txt";

          open(DECOY_FILE,">$decoy_cluster_file");
          close DECOY_FILE;
          
          #### Fetch names of reference genomes from the database #####
          
          my $reference_genome_name_sql="Select distinct(abbreviation) from OrganismInfo where genome_type='reference'";
          
          my $list_reference_genome=SQLiteDB::get_record($db_dir,$db_name,$reference_genome_name_sql);
          
          my @reference_genome=();
          
          if(scalar(@{$list_reference_genome})>0){
          
              foreach my $row(@{$list_reference_genome}){
              
                   my $genome_abbreviation=shift(@{$row});
                   
                   push(@reference_genome,$genome_abbreviation);
              }         
          }else{        
             die "Cannot find any defined reference genome in the database\n";
          }
          
          ############# Prepare Sequence Database for comparing reference sequences #######
          
          my $seq_database_file="$out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQ_DB}";
          
          open(SEQ_DB,">$seq_database_file");
          
          my $search_genome=join("','",@reference_genome);
          
          my $get_all_sequence_sql="Select * from ProteinSequence where genome_abbreviation IN ('$search_genome')";
          
          my $reference_protein_sequences=SQLiteDB::get_record($db_dir,$db_name,$get_all_sequence_sql);
          
          if(scalar(@{$reference_protein_sequences})>0){
          
                foreach my $row(@{$reference_protein_sequences}){
                
                        my $protein_index=shift(@{$row});
                        my $protein_id=shift(@{$row});
                        my $genome_abbreviation=shift(@{$row});
                        my $seq_type=shift(@{$row});
                        my $seq_length=shift(@{$row});
                        my $seq=shift(@{$row});
                        
                        print SEQ_DB ">$genome_abbreviation|$protein_id\n$seq\n";              
                }
          }
          close SEQ_DB;
          
          
          ######### RUN PHMMER SCAN FOR SEQUENCES FROM EACH REFERENCE GENOME ######
          
          foreach my $genome_name(@reference_genome){
          
              system("perl RunHMMER.pl $genome_name $hmm_program $db_dir $db_name $seq_database_file $tmp_dir 5000 1000 $hmmer_opt $cpu_core");
              
              system("cp $tmp_dir/dom_phmmer_$genome_name.out  $out_dir{hmm_domout_dir}");
              
              system("cp $tmp_dir/full_phmmer_$genome_name.out $out_dir{hmm_fullout_dir}");
     
              print STDOUT "Parsing $genome_name PHMMER scan result\n";

              my($read_phmmer_output)=Hmmer::Read_domain_table($hmm_program,$genome_name,$out_dir{hmm_domout_dir});

              my($sequence_alignment)=Hmmer::Read_aligned_sequence($hmm_program,$genome_name,$out_dir{hmm_fullout_dir},$read_phmmer_output); 

              FilterPair::getHitForQuery($genome_name,$config_param{PARSE_HMMER},$out_dir{similarity_dir},$out_dir{all_similarity_dir},$out_dir{chimeric_similarity_dir},$read_phmmer_output,$sequence_alignment); 
 
              print "Loading Similarity data in database\n"; 
              
              my $load_data_sql_stmt="INSERT or REPLACE into Similarity(query_id, subject_id, query_length, subject_length, total_domain, high_scoring_domain, qstart, qend, sstart, send, evalue, bitscore, percent_identity, percent_similarity, query_coverage, subject_coverage, pair_relation, note) VALUES
                                      (?, ?, ?, ? ,? ,? ,?, ? ,? ,? ,? ,? ,? ,? ,? ,? ,?, ?)";

              SQLiteDB::load_data($db_dir,$db_name,"$out_dir{all_similarity_dir}/allhit_$genome_name.txt",$load_data_sql_stmt);
              
          }

          system("cat $out_dir{similarity_dir}/* > $out_dir{mcl_dir}/$similarity_file");
          
          ##### MCL CLUSTERING ######          
          my $mcl_output=MclClustering::run_mcl($similarity_file,$out_dir{mcl_dir},$mcl_inflation_value);

       	  my $cluster_num = localtime();
             $cluster_num=~s/\:/\_/g;
             $cluster_num=~s/\s+/\_/g;      
          
          my ($cluster_file,$group_for_query)=AddGroupID::add_id($mcl_output,$cluster_num,$hmm_program,\%homologue_group,$out_dir{cluster_dir},$out_dir{mcl_dir});

          my $count_groups=keys(%{$group_for_query});
          
          system("cp $cluster_file $out_dir{result_dir}/reference_model_cluster_$cluster_num.txt");
          
          $cluster_file="$out_dir{result_dir}/reference_model_cluster_$cluster_num.txt";
          
          print "Creating HMM Model for gene cluster\n";
          
          my($ref_group_gene)=GroupModel($db_dir,$db_name,$cluster_file,$count_groups,\%config_param,\%homologue_group,$group_for_query,\%group_gene,\%out_dir,$cpu_core,$tmp_dir);    
         
          system("cp $out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB} -t $out_dir{result_dir}");
      
          system("cp $out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQ_DB} -t $out_dir{result_dir}");
          
          my $end_reference_scan = time;
          
          my $run_time_reference=$end_reference_scan-$start_reference_scan;
          
          print "\n\nREFERENCE GENOME COMPARISION IS COMPLETED\n\n";
          
          #print LOG_TIME "$cluster_file\t$count_groups\tReference Scan\t$run_time_reference\t",int($run_time_reference/ 3600),"h:",int(($run_time_reference % 3600) / 60),"m:",int($run_time_reference % 60),"s\n";            
    } 
}


##### PHASE II: ##### DENOVO PROTEIN FAMILY PREDICTION USING SEED / REFERENCE PROTEIN FAMILY ######

sub run_homolog_scan {

    $db_dir=(shift);
    $db_name=(shift);
    my %config_param=%{(shift)};
    my %out_dir=%{(shift)};
    my $output_dir=(shift);
    my $tmp_dir=(shift);
    
    my %genome=();
    
if($config_param{ACTIVATE_ANALYSIS}->{PREDICT_HMM}=~/YES/i){

    ##### DECLARE GLOBAL VARAIBLES ####
    my $cpu_core=$config_param{PARAMETERS}->{PARALLEL_CPU_CORE};
    my $hmmer_opt=$config_param{HMMER_PARAMETERS}->{HMMER_OPT};
    my $similarity_file="SimilarityPair.txt";
    my $count_query=0;
     
    ######## FETCH Abbreviation for Query Genome from the database #######
    
    my $query_genome_name_sql="Select distinct(abbreviation) from OrganismInfo where genome_type='query'";
          
    my $list_query_genome=SQLiteDB::get_record($db_dir,$db_name,$query_genome_name_sql);
          
    my @query_genome=();
          
    if(scalar(@{$list_query_genome})>0){
          
        foreach my $row(@{$list_query_genome}){
              
            my $genome_abbreviation=shift(@{$row});
                   
            push(@query_genome,$genome_abbreviation);
        }          
    }else{        
        die "Cannot find any defined reference genome in the database\n";
    } 


    open(LOG_TIME,">$out_dir{tmp_log}/execution_time.log");
    
    ###### READ information from user defined cluster file #######
    
        print "CHECKING INPUT PARAMTERS.....\n";
        
        my($group_file,$ref_group_gene,$homologue_group)=Initialize::check_input(\%config_param,\%out_dir); 
        
        my %group_gene=%{$ref_group_gene};   
               
        my %homologue_group=%{$homologue_group};
        
        my $group_for_query=$homologue_group;
        
        my $cluster_file=$group_file;                                      
       
        my $count_groups=keys(%{$group_for_query});
        
    ############ PREPARE HMM MODEL DATABASE AND SINGLETON SEQUENCE DATABASE ###########    

     if(($config_param{GROUP}->{MODEL_DB_FILE}) and ($config_param{GROUP}->{SINGLETON_DB_FILE})){

          if(-s "$config_param{GROUP}->{MODEL_DB_FILE}" and -s "$config_param{GROUP}->{SINGLETON_DB_FILE}"){    
           
               print "Copying HMM Model Database into the project directory\n";
                            
                  unless(-s "$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}"){

                      system("cp $config_param{GROUP}->{MODEL_DB_FILE} -t $out_dir{hmm_db_dir}");

                      my $hmm_db_file=basename("$config_param{GROUP}->{MODEL_DB_FILE}");

                      system("mv $out_dir{hmm_db_dir}/$hmm_db_file $out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}");
                 
                  }else{
                    
                      print "\nHMM model database with same name already exist in directory $out_dir{hmm_db_dir}\n\n";

                      print "Do you want to continue copying the file [Y or N]?\n\n";

                      my $continue_copy=<STDIN>;

                      if($continue_copy=~/(Y)/i){

                          system("cp $config_param{GROUP}->{MODEL_DB_FILE} -t $out_dir{hmm_db_dir}");
                          my $hmm_db_file=basename("$config_param{GROUP}->{MODEL_DB_FILE}");
                          system("mv $out_dir{hmm_db_dir}/$hmm_db_file $out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}");

                      }else{
                         print "SKIP COPYING HMM MODEL DATABASE FILE \n\n";
                      }
                  }
           
               print "Copying Singleton Sequence Database into the project directory\n";

                  unless(-s "$out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQ_DB}"){

                       system("cp $config_param{GROUP}->{SINGLETON_DB_FILE} -t $out_dir{hmm_db_dir}");

                       my $singleton_db_file=basename("$config_param{GROUP}->{SINGLETON_DB_FILE}");

                       system("mv $out_dir{hmm_db_dir}/$singleton_db_file $out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQ_DB}");             
                  
                  }else{
                       
                       print "\nSingleton Sequence database with same name already exist in directory $out_dir{hmm_db_dir}\n\n";

                       print "Do you want to continue copying the file [Y or N]?\n\n";

                       my $continue_copy=<STDIN>;

                       if($continue_copy=~/(Y)/i){

                          system("cp $config_param{GROUP}->{SINGLETON_DB_FILE} -t $out_dir{hmm_db_dir}");
                          my $singleton_db_file=basename("$config_param{GROUP}->{SINGLETON_DB_FILE}");
                          system("mv $out_dir{hmm_db_dir}/$singleton_db_file $out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQ_DB}");

                       }else{
                         print "SKIP COPYING SINGLETON DATABASE FILE\n\n";
                       }
                  }
               ##### Printing out HMM Files for each group #####

            print "Do you want to prepare hmm files for each group [Y or N]?\n\n";

            my $print_hmm_file=<STDIN>;

            if($print_hmm_file=~/(Y)/i){

               print "Printing HMM Files\n";
               my @hmm_db=();
               tie @hmm_db, 'Tie::File', "$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}";           
              
               my $hmm_block='';
               my $group_id='';

              for(my $l=0;$l<=scalar(@hmm_db);$l++){

                   if($hmm_db[$l] and $hmm_db[$l] eq "//"){

                      $hmm_block=$hmm_block.$hmm_db[$l];
                      
                      if($homologue_group{$group_id}){  
                          print "$group_id\n";                    
                          open(HMM_BLOCK,">$out_dir{hmm_file_dir}/$group_id.hmm");
                          print HMM_BLOCK "$hmm_block\n";
                          close HMM_BLOCK;
                      }   

                      $hmm_block='';
                      $group_id='';

                   }elsif($hmm_db[$l]){

                      $hmm_block=$hmm_block.$hmm_db[$l]."\n"; 
  
                      if($hmm_db[$l]=~/^(NAME)(\s+)(\w+)/){                          
                          $hmm_db[$l]=~/^(NAME)(\s+)(\w+)/;
                          $group_id=$3; 
                       }                 
                   }
               }
            }
           
            print "Do you want to prepare singleton files for each group [Y or N]?\n\n";

            my $print_singleton_file=<STDIN>;

            if($print_singleton_file=~/(Y)/i){

               #####Printing Singleton Files for each group ####
               print "Printing Singleton Files\n"; 
               my @singleton_db=();              

               tie @singleton_db, 'Tie::File', "$out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQ_DB}";
      
               my $single_db=join("\n",@singleton_db);
               my @singleton_group=split(">",$single_db);
               undef($single_db);

               foreach my $singleton_seq(@singleton_group){
                   if($singleton_seq=~/^(Group)(\d+)/){
                       $singleton_seq=~/^(Group)(\d+)/;
                       my $group_id=$1.$2;  
                        
                       open(SINGLETON,">$out_dir{singleton_group}/$group_id.fasta");
                       print SINGLETON ">$singleton_seq\n";
                       print ">$singleton_seq\n";
                       close SINGLETON;                                        
                   }
               } 
            }

          }
     }


     unless($config_param{GROUP}->{MODEL_DB_FILE} and $config_param{GROUP}->{SINGLETON_DB_FILE}){

       %group_gene=();  

       ########## Call Group Model Function ######
       print "Creating HMM Model for gene cluster\n";

       my($ref_group_gene)=GroupModel($db_dir,$db_name,$cluster_file,$count_groups,\%config_param,\%homologue_group,$group_for_query,\%group_gene,\%out_dir,$cpu_core,$tmp_dir);    

       %group_gene=%{$ref_group_gene};                 
     }
       
    my $jump_to_step=4;
    
    
 ######## ITERATIVE ANALYSIS FOR EACH ADDITIONAL GENOME ##################
       
      foreach my $genome_name(@query_genome){

             ###############

             my $start_query_scan = time;                       
             
             undef %homologue_group;

             my %genes_in_group=();
             
             print "READ CLUSTER $cluster_file\n"; 
                       
             my($homologue_group,$query_status)=Cluster::read_cluster($cluster_file,$genome_name);
             %homologue_group=%{$homologue_group};

             my($genes_in_group)=Cluster::readGroupForGene($cluster_file);
             %genes_in_group=%{$genes_in_group};

             #### count number of clusters for each analyzed genome #####
             my $count_groups=keys(%{$group_for_query});   
             print "Number of groups: $count_groups\n";

             ### If any sequence from query file is already present in cluster, do not process that query further #####
             if($query_status==1){ 
                  print "Warning: Cluster file contains sequences from $genome_name, skip this genome\n";
                  $count_query++;
                  next;
             } 
 
             if($jump_to_step==4){
                 $jump_to_step=0;
                 goto(BUILD_DB);                
             }elsif($jump_to_step==5){
                 $jump_to_step=0;
                 goto(BUILD_SINGLETON_ONLY); 
             }

             ##### CALL Group Model Function #######         
             my($ref_group_gene)=GroupModel($db_dir,$db_name,$cluster_file,$count_groups,\%config_param,\%homologue_group,$group_for_query,\%group_gene,\%out_dir,$cpu_core,$tmp_dir);    
         
             #### UPDATE GLOBAL VARIABLES ##### 
             undef %group_gene;          

             %group_gene=%{$ref_group_gene};             
             
             ###### JUMP HERE IF JUMP_STEP IS EQUAL TO 4 ########
             BUILD_DB:
             
             ######### MODEL DATABASE ###########################
             
             print "BUILD HMM MODEL DATABASE\n";
             
             if(-e "$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}.h3i"){
                   unlink("$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}.h3i");
                   unlink("$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}.h3f");
                   unlink("$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}.h3m");
                   unlink("$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}.h3p");
             }

             system("hmmpress $out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}")==0 or die "$?\n";
             
             ######## SINGLETON DATABASE #########################
             BUILD_SINGLETON_ONLY:
             print "BUILD HMM SINGLETON DATABASE\n";
              
             system("cp $out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQ_DB} $out_dir{hmm_db_dir}/SINGLETON_SEQ_DB");

             my $unanalyzed_gene=0;

             open(SINGLETON_DB,">>$out_dir{hmm_db_dir}/SINGLETON_SEQ_DB");

             my $get_sequence_sql="Select * from ProteinSequence where genome_abbreviation='$genome_name'";
          
             my $query_protein_sequences=SQLiteDB::get_record($db_dir,$db_name,$get_sequence_sql);
          
                    if(scalar(@{$query_protein_sequences})>0){
          
                         foreach my $row(@{$query_protein_sequences}){
                
                                 my $protein_index=shift(@{$row});
                                 my $protein_id=shift(@{$row});
                                 my $genome_abbreviation=shift(@{$row});
                                 my $seq_type=shift(@{$row});
                                 my $seq_length=shift(@{$row});
                                 my $seq=shift(@{$row});
                                 my $gene_id=$genome_abbreviation."|".$protein_id;
                                 
                                 if($genes_in_group{$gene_id}){
                                    next;
                                 }

                                 print SINGLETON_DB ">$genome_abbreviation|$protein_id\n$seq\n";

                                 $unanalyzed_gene++;
                         }
                    }
                    
                    close SINGLETON_DB;

             if($unanalyzed_gene==0){
                $jump_to_step=5;
                goto SKIP_FURTHER_PROCESS;
             }

             ######## HMMSCAN ANALYSIS ##########################
             print "RUN HMMSCAN for $genome_name\n";
             
             my $start_hmmscan = time;
             
             my $program="hmmscan";

             system("perl RunHMMER.pl $genome_name $program $db_dir $db_name $out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB} $tmp_dir 5000 1000 $hmmer_opt $cpu_core");
             
             #if(-z "$out_dir{tmp_log}/dom_hmmscan_$genome_name.out" or -z "$out_dir{tmp_log}/full_hmmscan_$genome_name.out"){
             
             #  print STDERR "HMMSCAN process failed to generate output for $genome_name\n";
             #  exit;
             #}

             print "Compiling HMM Search Output\n";
             
             system("cat $out_dir{tmp_log}/dom_hmmscan_$genome_name.out > $out_dir{hmm_domout_dir}/dom_hmmscan_$genome_name.out");
             
             system("cat $out_dir{tmp_log}/full_hmmscan_$genome_name.out > $out_dir{hmm_fullout_dir}/full_hmmscan_$genome_name.out");
             
             unlink("$out_dir{tmp_log}/dom_hmmscan_$genome_name.out");
             
             unlink("$out_dir{tmp_log}/full_hmmscan_$genome_name.out"); 

             my $end_hmmscan = time;
             
             my $run_hmmscan=$end_hmmscan-$start_hmmscan;
             
             print LOG_TIME "$cluster_file\thmmscan\t$run_hmmscan\t",int($run_hmmscan/ 3600),"h:",int(($run_hmmscan % 3600) / 60),"m:",int($run_hmmscan % 60),"s\n"; 
           
             ######## PHMMER ANALYSIS FOR SINGLETON FAMILY ##########
             print "Run Phmmer for Singlton groups\n";
             
             my $start_phmmer = time;     
                           
             $program="phmmer";

             system("perl RunHMMER.pl $genome_name $program $db_dir $db_name $out_dir{hmm_db_dir}/SINGLETON_SEQ_DB $tmp_dir 5000 1000 $hmmer_opt $cpu_core");

             print "Compiling PHMMER Search Output\n";
 
             #if(-z "$out_dir{tmp_log}/dom_phmmer_$genome_name.out" or -z "$out_dir{tmp_log}/full_phmmer_$genome_name.out"){
             #  print STDERR "PHMMER process failed to generate output for $genome_name\n";
             #  exit;
             #}

              system("cat $out_dir{tmp_log}/dom_phmmer_$genome_name.out > $out_dir{hmm_domout_dir}/dom_phmmer_$genome_name.out");
              
              system("cat $out_dir{tmp_log}/full_phmmer_$genome_name.out > $out_dir{hmm_fullout_dir}/full_phmmer_$genome_name.out");
              
              unlink("$out_dir{tmp_log}/dom_phmmer_$genome_name.out");
              
              unlink("$out_dir{tmp_log}/full_phmmer_$genome_name.out");

              my $end_phmmer= time;
              my $run_phmmer=$end_phmmer-$start_phmmer;
              
              print LOG_TIME "$cluster_file\tphmmer\t$run_phmmer\t",int($run_phmmer/ 3600),"h:",int(($run_phmmer % 3600) / 60),"m:",int($run_phmmer % 60),"s\n";                         
                

              ####### PARSE OUTPUT FILE FROM HMMSCAN AND PHMMER SEARCH #########             
              
              print STDOUT "PARSING OUTPUT FOR $genome_name\n"; 

              my %domain_table=();
              my %sequence_alignment=();
              my $generate_output_flag=0;

       if(-s "$out_dir{tmp_log}/dom_hmmscan_$genome_name.out" and -s "$out_dir{tmp_log}/full_hmmscan_$genome_name.out"){

              print STDOUT "READING HMM SCAN DOMAIN TABLE\n";
              my($read_hmmscan_output)=Hmmer::Read_domain_table("hmmscan",$genome_name,$out_dir{hmm_domout_dir});

              print STDOUT "READING HMMSCAN ALIGNMENT OUTPUT\n";
              my($hmmscan_sequence_alignment)=Hmmer::Read_aligned_sequence("hmmscan",$genome_name,$out_dir{hmm_fullout_dir},$read_hmmscan_output); 

              %domain_table=%{merge(\%domain_table, $read_hmmscan_output)};
              
              my $count_keys1=keys(%domain_table);
              
              print "HMMSCAN LINES:\t", $count_keys1,"\n";

              %sequence_alignment=%{merge(\%sequence_alignment, $hmmscan_sequence_alignment)};
              
              my $count_keys3=keys(%sequence_alignment);
              
              print "HMMSCAN LINES:\t", $count_keys3,"\n";

              if($count_keys1>0){
                $generate_output_flag=1;
              }

        }else{
             print STDERR "No hits in HMM model database for $genome_name\n";
        }

        if(-s "$out_dir{tmp_log}/dom_phmmer_$genome_name.out" and -s "$out_dir{tmp_log}/full_phmmer_$genome_name.out"){    
  
              print STDOUT " READING PHMMER SCAN DOMAIN TABLE\n";
              my($read_phmmer_output)=Hmmer::Read_domain_table("phmmer",$genome_name,$out_dir{hmm_domout_dir});

              print STDOUT "READING PHMMER ALIGNMENT OUTPUT\n";
              my($phmmer_sequence_alignment)=Hmmer::Read_aligned_sequence("phmmer",$genome_name,$out_dir{hmm_fullout_dir},$read_phmmer_output); 

              %domain_table=%{merge(\%domain_table, $read_phmmer_output)};
              
               my $count_keys2=keys(%domain_table);
               
              print "PHMMER LINES:\t", $count_keys2,"\n";
 
              %sequence_alignment=%{merge(\%sequence_alignment, $phmmer_sequence_alignment)};
              
              my $count_keys4=keys(%sequence_alignment);
              
              print "PHMMER LINES:\t", $count_keys4,"\n";

              if($count_keys2>0){
                $generate_output_flag=1;
              }

         }else{

            print STDERR "No hits in HMM Singleton database for $genome_name\n";
         }

             if($generate_output_flag==0){
                $jump_to_step=5;
                goto SKIP_FURTHER_PROCESS;
             }

 
              print STDOUT "FILTERING PAIRS\n";
              
              FilterPair::getHitForQuery($genome_name,$config_param{PARSE_HMMER},$out_dir{similarity_dir},$out_dir{all_similarity_dir},$out_dir{chimeric_similarity_dir},\%domain_table,\%sequence_alignment);       
              
              ###### PREPARE BEST MATCH SIMILAIRTY FILE #########   
              
              system("cat $out_dir{similarity_dir}/besthit_$genome_name.txt > $out_dir{mcl_dir}/$similarity_file");

              ###### LOAD ALL MATCH SIMILARITY FILE IN DATABASE #####
              
              print "Loading Similarity data in database\n";
     
              my @full_similarity=();   
                         
              tie @full_similarity, 'Tie::File', "$out_dir{all_similarity_dir}/allhit_$genome_name.txt";  

              my $load_data_sql_stmt="INSERT or REPLACE into Similarity(query_id, subject_id, query_length, subject_length, total_domain, high_scoring_domain, qstart, qend, sstart, send, evalue, bitscore, percent_identity, percent_similarity, query_coverage, subject_coverage, pair_relation, note) VALUES
                            (?, ?, ?, ? ,? ,? ,?, ? ,? ,? ,? ,? ,? ,? ,? ,? ,?, ?)";
              SQLiteDB::load_from_array($db_dir,$db_name,\@full_similarity,$load_data_sql_stmt);

              $group_for_query={};
                
              my $hmm_program="hmmscan"; 

              my $mcl_inflation_value=$config_param{PARAMETERS}->{MCL_INFLATION_VALUE};
              
              my $mcl_output=MclClustering::run_mcl($similarity_file,$out_dir{mcl_dir},$mcl_inflation_value);

       	      my $cluster_num = localtime();
                 $cluster_num=~s/\:/\_/g;
                 $cluster_num=~s/\s+/\_/g;

              ($cluster_file,$group_for_query)=AddGroupID::add_id($mcl_output,$cluster_num,$hmm_program,\%homologue_group,$out_dir{cluster_dir},$out_dir{mcl_dir});
              
              $count_query++;
              
              $genome{$genome_name}=$genome_name;  

              $count_groups=keys(%{$group_for_query});

              print "Number of groups in query $genome_name: $count_groups\n";

              ##### point here only if there are no genes in query genome for analysis ####
              SKIP_FURTHER_PROCESS:

              unlink("$out_dir{tmp_log}/dom_hmmscan_$genome_name.out");
             
              unlink("$out_dir{tmp_log}/full_hmmscan_$genome_name.out"); 

              unlink("$out_dir{tmp_log}/dom_phmmer_$genome_name.out");
              
              unlink("$out_dir{tmp_log}/full_phmmer_$genome_name.out");
 
              my $end_query_scan = time;  
              
              my $run_time_query=$end_query_scan-$start_query_scan;
              
              print LOG_TIME "$cluster_file\t$count_groups\tQuery Scan\t$run_time_query\t",int($run_time_query/ 3600),"h:",int(($run_time_query % 3600) / 60),"m:",int($run_time_query % 60),"s\n";                      
      }

      ###### MODEL BUILDING AFTER PROCESSING LAST QUERY IN THE QUEUE #####
      print "Computing Cluster and Distance for analysis\n";  
       
       	my $name_value = localtime();
           $name_value=~s/\:/\_/g;
           $name_value=~s/\s+/\_/g; 

       system("cp $cluster_file $out_dir{result_dir}/model_cluster_$name_value.txt");

       $cluster_file="$out_dir{result_dir}/model_cluster_$name_value.txt";
      
      ($homologue_group,my $query_status)=Cluster::read_cluster($cluster_file,'xxx');
      
      %homologue_group=%{$homologue_group};
      
      $count_groups=keys(%{$group_for_query}); 

      ##### CALL GROUP MODEL FUNCTION ########   
      ($ref_group_gene)=GroupModel($db_dir,$db_name,$cluster_file,$count_groups,\%config_param,\%homologue_group,$group_for_query,\%group_gene,\%out_dir,$cpu_core,$tmp_dir);    
            
      undef %group_gene;

      my $final_model_db=$config_param{DATABASE}->{MODEL_DB}."_".$name_value;
      
      my $final_singleton_db=$config_param{DATABASE}->{SEQ_DB}."_".$name_value;

      system("cp $out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB} -t $out_dir{result_dir}");
      
      system("cp $out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQ_DB} -t $out_dir{result_dir}"); 

      system("mv $out_dir{result_dir}/$config_param{DATABASE}->{MODEL_DB} $out_dir{result_dir}/$final_model_db"); 

      system("mv $out_dir{result_dir}/$config_param{DATABASE}->{SEQ_DB} $out_dir{result_dir}/$final_singleton_db"); 

      print "\n\nITERATIVE PROTEIN COMPARISION IS COMPLETED\n\n";

      print "[RESULT FILES]\n\n";

      print "CLUSTER FILE: $cluster_file\n";
      print "HMM MODEL DATABASE: $out_dir{result_dir}/$final_model_db\n";
      print "SINGLETON SEQUENCE DATABASE: $out_dir{result_dir}/$final_singleton_db\n\n";
  }
  
} #### END OF BLOCK HOMOLOG SCAN  


####### PARTIAL GENE MAPPING and SUPER-HOMOLOG PREDICTION ######

sub run_super_homolog_prediction {

    $db_dir=(shift);
    $db_name=(shift);
    my $homolog_database_name=(shift);
    my %config_param=%{(shift)};
    my %out_dir=%{(shift)};
    my $output_dir=(shift);
    my $tmp_dir=(shift);
    my $name_value=(shift);
       
	if($config_param{ACTIVATE_ANALYSIS}->{PREDICT_SUPER_HOMOLOG}=~/YES/i){
	
	 ##### DECLARE GLOBAL VARAIBLES ####
    
         my $cluster_file=$config_param{GROUP}->{HMM_CLUSTER_FILE};
         
         my $parameters=$config_param{PARAMETER};          
	
	 print "CLUSTER HMM FAMILY INTO HOMOLOG FAMILY \n";  

         my($homolog_cluster_file,$GenetoSuperfamily_file)=PartialMap::mapPartialSequence($cluster_file,$parameters,\%out_dir,$db_dir,$db_name,$homolog_database_name,$name_value);

         print "\n\nPARTIAL PROTEIN MAPPING IS COMPLETED\n\n";

         print "[RESULT FILES]\n\n";

         print "HOMOLOG CLUSTER FILE: $homolog_cluster_file\n";
         print "GENE To SUPERFAMILY MAP FILE: $GenetoSuperfamily_file\n\n";
	}
}

##################### ORTHOLOG PREDICTION PHASE ###################################################

sub run_ortholog_prediction{
    
    $db_dir=(shift);
    $db_name=(shift);
    my($old_ortholog_database)=(shift);
    my($new_ortholog_database)=(shift);
    my %config_param=%{(shift)};
    my %out_dir=%{(shift)};
    my $output_dir=(shift);
    my $tmp_dir=(shift);
    my $name_value=(shift);
    
    
    if($config_param{ACTIVATE_ANALYSIS}->{PREDICT_ORTHOLOG}=~/YES/i or $config_param{ACTIVATE_ANALYSIS}->{CLUSTER_ORTHOLOG}=~/^YES$/i){
    
         my $cluster_file=$config_param{GROUP}->{HMM_CLUSTER_FILE};
         my $homolog_cluster_file=$config_param{GROUP}->{HOMOLOG_CLUSTER_FILE};
         my $cpu_core=$config_param{PARAMETERS}->{PARALLEL_CPU_CORE};
         my $ortholog_divergence_threshold=$config_param{PARAMETERS}->{ORTHOLOG_DIVERGENCE_THRESHOLD};
         my $inparalog_divergence_threshold=$config_param{PARAMETERS}->{INPARALOG_DIVERGENCE_THRESHOLD};
         my $alignment_overlap_threshold=$config_param{PARAMETERS}->{GLOBAL_ALIGNMENT_OVERLAP_THRESHOLD};
         my $pair_identity_threshold=$config_param{PARAMETERS}->{PAIRWISE_SEQUENCE_IDENTITY_THRESHOLD};
         my $pair_similarity_threshold=$config_param{PARAMETERS}->{PAIRWISE_SEQUENCE_SIMILARITY_THRESHOLD};
     
         print "$cluster_file\n$homolog_cluster_file\n";
 
         unless(-s "$cluster_file" and -s "$homolog_cluster_file"){     
             die "ERROR: Cannot read group files\n";
         }

            my($homolog_cluster,$query_status)=Cluster::read_cluster($homolog_cluster_file,'xxx');

            my(%homolog_cluster)=%{($homolog_cluster)};   
         
            my @query_file=();

            foreach my $cluster_id(keys %homolog_cluster){

                my $cluster_gene=$homolog_cluster{$cluster_id};

                my @cluster_gene=split(/\t/,$cluster_gene);

                my $list_gene_id=join("\t",@cluster_gene);

                my $size_cluster=scalar(@cluster_gene);

                if($size_cluster > 1){      
                    push(@query_file,$cluster_id);
                }
            }
        
         
         if($config_param{ACTIVATE_ANALYSIS}->{PREDICT_ORTHOLOG}=~/YES/i){

            print "Predict Ortholog Pairs\n";  
         
            my $fork=Parallel::ForkManager->new($cpu_core);

            foreach my $cluster_id(@query_file){

                $fork->start and next;
     
                system("perl PredictOrtholog.pl $cluster_id $db_dir $db_name $old_ortholog_database $homolog_cluster_file $tmp_dir $out_dir{pair_distance} $out_dir{pair_ortholog} $out_dir{pair_inparalog} $out_dir{homolog_alignment} $out_dir{pair_alignment} $ortholog_divergence_threshold $inparalog_divergence_threshold $alignment_overlap_threshold $pair_identity_threshold $pair_similarity_threshold");
           
                $fork->finish();
            } 
            $fork->wait_all_children; 


              print "Loading Distance Pairs\n";
         
              opendir(DIST,"$out_dir{pair_distance}"); 
              my @dist_file=readdir(DIST);
         
              my $dist_pair_insert_stmt="INSERT INTO DistancePair(taxonA,idA,taxonB,idB,divergence,homolog_cluster_id) VALUES (?, ?, ?, ?, ?, ?)";
       
              foreach my $dist_file(@dist_file){
         
                  SQLiteDB::load_data($db_dir,$new_ortholog_database,"$out_dir{pair_distance}/$dist_file",$dist_pair_insert_stmt);
              }
         
         
              print "Loading Ortholog Pairs\n";
         
              opendir(ORTHO,"$out_dir{pair_ortholog}"); 
              my @ortho_file=readdir(ORTHO);
         
              my $ortho_insert_stmt="INSERT INTO Ortholog(taxonA,idA,taxonB,idB,divergence,homolog_cluster_id) VALUES (?, ?, ?, ?, ?, ?)";
       
              foreach my $ortho_file(@ortho_file){
         
                  SQLiteDB::load_data($db_dir,$new_ortholog_database,"$out_dir{pair_ortholog}/$ortho_file",$ortho_insert_stmt);
              }
                
               print "Loading InParalog Pairs\n";
         
               opendir(INPARA,"$out_dir{pair_inparalog}"); 
               my @inpara_file=readdir(INPARA);
         
               my $inpara_insert_stmt="INSERT INTO Inparalog(taxonA,idA,taxonB,idB,divergence,min_ortholog_divergence,homolog_cluster_id) VALUES (?, ?, ?, ?, ?, ?, ?)";
       
               foreach my $inpara_file(@inpara_file){
         
                  SQLiteDB::load_data($db_dir,$new_ortholog_database,"$out_dir{pair_inparalog}/$inpara_file",$inpara_insert_stmt);
               }
                 
               print "Loading Pairwise Alignment\n";
         
               opendir(PAIR_ALN,"$out_dir{pair_alignment}"); 
               my @pair_aln_file=readdir(PAIR_ALN);
         
               my $pair_alignment_insert_stmt="INSERT INTO PairwiseAlignment(idA, idB, sequence_identity, sequence_similarity, homolog_cluster_id) VALUES (?, ?, ? , ?, ?)";
       
               foreach my $pair_aln_file(@pair_aln_file){
         
                  SQLiteDB::load_data($db_dir,$new_ortholog_database,"$out_dir{pair_alignment}/$pair_aln_file",$pair_alignment_insert_stmt);
               }
                    
               print "Loading Homolog Alignment\n";
         
               opendir(ALN,"$out_dir{homolog_alignment}"); 
               my @aln_file=readdir(ALN);
         
               my $load_alignment_stmt="INSERT INTO Alignment (feature_id, genome_abbreviation, seq_type, alignment_id, alignment_length, alignment_sequence) VALUES(?,?,?,?,?,?)";
       
               foreach my $aln_file(@aln_file){
                  SQLiteDB::load_data($db_dir,$new_ortholog_database,"$out_dir{homolog_alignment}/$aln_file",$load_alignment_stmt);      
               } 

                print "\n\nORTHOLOG PAIR PREDICTION IS COMPLETED\n\n";
         }  
          

         if($config_param{ACTIVATE_ANALYSIS}->{CLUSTER_ORTHOLOG}=~/^YES$/i){   
    
              print "ORTHOLOG CLUSTERING\n";

              my($homolog_cluster,$query_status)=Cluster::read_cluster($homolog_cluster_file,'xxx');

              my(%homolog_cluster)=%{($homolog_cluster)};
 
              my $mcl_inflation_value=$config_param{PARAMETERS}->{MCL_INFLATION_VALUE};

              foreach my $cluster_id(keys %homolog_cluster){

                   my $cluster_gene=$homolog_cluster{$cluster_id};
 
                   my @cluster_gene=split(' ',$cluster_gene);

                   my $size_cluster=scalar(@cluster_gene);
                  
                   if($size_cluster>1){ 
       
                       system("perl OrthologClustering_MCL.pl $cluster_id $out_dir{pair_ortholog} $out_dir{pair_inparalog} $out_dir{ortholog_cluster} $tmp_dir $mcl_inflation_value"); 

                       #system("perl OrthologClustering.pl $homolog_cluster_file $out_dir{ortholog_cluster} $db_dir $db_name $tmp_dir $cpu_core $transitivity_threshold"); 
                   }else{                     
                       open(CLUSTER,">$out_dir{ortholog_cluster}/orthoCluster_$cluster_id.txt");    
                       print CLUSTER "$cluster_id.1:\t$cluster_gene\n";
                       print "$cluster_id.1:\t$cluster_gene\n";
                       close CLUSTER;
                   }
              }    
 
              system("find $out_dir{ortholog_cluster} -name '*' -type f | xargs cat > $out_dir{result_dir}/ortholog_cluster_$name_value.txt");

              print "\n\nORTHOLOG FAMILY CLUSTERING IS COMPLETED\n\n";

              print "[RESULT FILES]\n\n";

              print "ORTHOLOG CLUSTER FILE: $out_dir{result_dir}/ortholog_cluster_$name_value.txt\n";
         }
         
    }
             
} ##### END OF ORTHOLOG PREDICTION BLOCK
  


###### GROUP MODEL FUNCTION ######
sub GroupModel {
    
    my($db_dir)=(shift);
    my($db_name)=(shift);
    my($cluster_file)=(shift);
    my($count_groups)=(shift);
    my(%config_param)=%{(shift)};
    my($homologue_group)=(shift);
    my($group_for_query)=(shift);
    my($group_gene)=(shift);    
    my(%out_dir)=%{(shift)};  
    my($process)=(shift);
    my($tmp)=(shift);

    my $hmm_file_dir=$out_dir{hmm_file_dir};
    
    my $singleton_group_dir=$out_dir{singleton_group};
    
    my $hmm_database="$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}"; 
    
    my $singleton_database="$out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQ_DB}";
    
    my $group_for_query_file="$tmp/group_for_query.txt";
    
    my $group_org="$tmp/group_org.txt";
    
    print "Select Proteins for HMM MODEL\n";
                       
     my $start_protein_selection = time;

     ####### print out group for query genome in tmp file ####
     open(QUERY_GROUP,">$group_for_query_file");
     foreach my $group_id(keys %{$group_for_query}){
       print QUERY_GROUP "$group_id\n";
     }
     close QUERY_GROUP;

     ######  print list of previously analyzed genome in tmp file #### 
     open(GENOME_LIST,">$group_org");
     foreach my $genome_name(keys %{$group_gene}){
       print GENOME_LIST "$genome_name\n";
     }
     close GENOME_LIST;

     system("perl BuildHMMModel.pl $cluster_file $group_for_query_file $group_org $db_dir $db_name $hmm_file_dir $singleton_group_dir $hmm_database $singleton_database $tmp $process");
     
     my $end_protein_selection = time;
     
     my $run_time_divergence=$end_protein_selection-$start_protein_selection;  
                    
    open(MODEL_DB,">$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}") or die "$?";
     
    close MODEL_DB;

    system("find $out_dir{hmm_file_dir} -name '*' -type f | xargs cat > $out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}");        
    
    system("find $out_dir{singleton_group} -name '*' -type f | xargs cat > $out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQ_DB}");  

    open(GENOME_LIST,"$tmp/group_org.txt") or die "$?";  
    my @genome_list=<GENOME_LIST>;
    close GENOME_LIST;
    my %group_gene=map{$_=>$_}@genome_list;         
    my $ref_group_gene=\%group_gene;

    return($ref_group_gene);
}

#### END #####
