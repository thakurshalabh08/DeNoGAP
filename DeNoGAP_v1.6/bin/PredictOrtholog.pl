#### Module to Predict Ortholog Genes #######
##### Author: Shalabh Thakur #################
##### Date: 28-SEP-2013 ######################

#!/usr/bin/perl
package Ortholog;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Bio::Tools::Phylo::Phylip::ProtDist;
use Bio::AlignIO;
use Bio::SeqIO;
use Cluster;
use SQLiteDB;
use Hash::Merge qw( merge );
use File::Path qw(remove_tree);
use File::Basename;

my($cluster_id)=(shift);
my($db_dir)=(shift);
my($db_name)=(shift);
my($old_ortholog_database)=(shift);
my($homolog_cluster_file)=(shift);
my($tmp_out_dir)=(shift);
my($distance_pair_dir)=(shift);
my($ortholog_pair_dir)=(shift);
my($inparalog_pair_dir)=(shift);
my($alignment_dir)=(shift);
my($pairwise_alignment_dir)=(shift);
my($divergence_threshold)=(shift);
my($inparalog_divergence_threshold)=(shift);
my($alignment_overlap_threshold)=(shift);
my($pair_identity_threshold)=(shift);
my($pair_similarity_threshold)=(shift); 

my %outgroup=();
my @process_cluster=();
my $distance_pair=undef;
my %db_distance_pair=();
my %pairwise_alignment=();


  print "$cluster_id\n";

##### GET LIST OF OUT-GROUP GENOMES ####

my($get_outgroup)="SELECT DISTINCT(abbreviation) from OrganismInfo WHERE outgroup='Yes'";
my ($outgroup_record)=SQLiteDB::get_record($db_dir,$db_name,$get_outgroup);

if(scalar($outgroup_record)>1){           
  foreach my $row(@{$outgroup_record}){           
     my $out_taxa=shift(@{$row});
     chomp($out_taxa);
     $outgroup{$out_taxa}=$out_taxa;
  }
}

#### READ FROM HOMOLOG CLUSTER FILE ####
my($homolog_cluster,$query_status)=Cluster::read_cluster($homolog_cluster_file,'xxx');
#my($groupforgene)=Cluster::readGroupForGene($group_file);
my(%homolog_cluster)=%{($homolog_cluster)};
#my(%groupforgene)=%{($groupforgene)};
my $count_cluster=keys(%homolog_cluster);

        
        ###### GET PROTEIN ID FOR CLUSTER ######
        my $cluster_gene=$homolog_cluster{$cluster_id};
        
        my @cluster_gene=split(' ',$cluster_gene);

        my $list_gene='';

        my $list_gene_2='';

         foreach my $gene_id(@cluster_gene){ 

             $gene_id=~/(\w+)(\|)(.+)/;

             $list_gene=$list_gene.",'".$3."'";

             $list_gene_2=$list_gene_2.",'".$gene_id."'";
         }
         
         $list_gene=~s/^\,//g;

         $list_gene_2=~s/^\,//g;
           
        
        ##### CREATE TEMPORARY DIRECTORY FOR THE PROCESS #####
        my $tmp_process_dir="$tmp_out_dir/$cluster_id";
        mkdir($tmp_process_dir);
        
        ####### GET HOMOLOG SEQUENCE MSA FROM DATABASE FOR CLUSTER ID #####
        if($old_ortholog_database eq "NULL"){
          goto SKIP_OLD_DB;
        }

        unless(-s "$old_ortholog_database"){

            die "Error $old_ortholog_database cannot be found\n";
        }

           my $old_ortholog_database_name=basename($old_ortholog_database);
           my $old_ortholog_database_dir=dirname($old_ortholog_database);

           my $get_aln_cluster="Select * from Alignment where alignment_id='$cluster_id'";
           my ($cluster_aln_record)=SQLiteDB::get_record($old_ortholog_database_dir,$old_ortholog_database_name,$get_aln_cluster);

           my $get_distance_stmt="Select * from  DistancePair where idA IN ($list_gene) or idB IN ($list_gene)";
           my ($distance_record)=SQLiteDB::get_record($old_ortholog_database_dir,$old_ortholog_database_name,$get_distance_stmt);

           my $get_ortholog_stmt="Select * from  Ortholog where homolog_cluster_id='$cluster_id'";
           my ($ortholog_record)=SQLiteDB::get_record($old_ortholog_database_dir,$old_ortholog_database_name,$get_ortholog_stmt);

           my $get_inparalog_stmt="Select * from  Inparalog where homolog_cluster_id='$cluster_id'";
           my ($inparalog_record)=SQLiteDB::get_record($old_ortholog_database_dir,$old_ortholog_database_name,$get_inparalog_stmt);

           my $get_pair_aln_stmt="Select * from  PairwiseAlignment where idA IN ($list_gene_2) or idB IN ($list_gene_2)";
           my ($pair_aln_record)=SQLiteDB::get_record($old_ortholog_database_dir,$old_ortholog_database_name,$get_pair_aln_stmt);

        

        #### get pairwise alignment from database for all genes in cluster ####
    
       if(scalar(@{$distance_record})>0){
    
          foreach my $row(@{$distance_record}){           
            my @data_column=@{$row};   
            my $idA=$data_column[0]."|".$data_column[1];
            my $idB=$data_column[2]."|".$data_column[3];
            my $dist_value=$data_column[4];
            
            $db_distance_pair{$idA}->{$idB}=$dist_value;
            $db_distance_pair{$idB}->{$idA}=$dist_value;
          }    
       }

       if(scalar(@{$pair_aln_record})>0){
          
             foreach my $row(@{$pair_aln_record}){  
          
                 my @column=@{$row};
              
                 $pairwise_alignment{$column[0]}->{$column[1]}->{identity}=$column[2]; 
                 $pairwise_alignment{$column[1]}->{$column[0]}->{identity}=$column[2];  

                 $pairwise_alignment{$column[0]}->{$column[1]}->{similarity}=$column[3]; 
                 $pairwise_alignment{$column[1]}->{$column[0]}->{similarity}=$column[3]; 
           }
       }

        SKIP_OLD_DB:  
        #### (1) Build Multiple Sequence Alignment for Homolog Group #######
        print "Making Alignment\n";

        my($family_sequence,$alignment_file)=create_msa($db_dir,$db_name,$cluster_id,$cluster_aln_record,$cluster_gene,$tmp_process_dir);

        if(!defined($alignment_file)){

           $alignment_file="$tmp_process_dir/$cluster_id.afa";

           open(ALN_SEQ,">$alignment_file"); 
           
           ##### print alignment file ###

           my $found_distance=1;

           foreach my $row(@{$cluster_aln_record}){            
              my @column=@{$row};
              my $seq_id=$column[0];
              my $aln_seq=$column[5];
              print ALN_SEQ ">$seq_id\n$aln_seq\n";   
 
              if(!defined($db_distance_pair{$seq_id})){
                $found_distance=0;
              }              
           } 
           close ALN_SEQ;

           $distance_pair=\%db_distance_pair;

          if($found_distance==1){
            goto STEP_4;
          }
        }
        
        #### (2) Convert format of alignment to phylip #####
        print "Making Phylip File\n";
        my($phylip_alignment_file,$id_map)=convertFastaToPhylip($alignment_file,"$tmp_process_dir/tmp_aln.afa");
        
        #### (3) Calculate Pairwise Genetic Distance for Homolog Sequences ####    
        print "Calculating Distance\n";
        $distance_pair=calculate_pairwise_distance($db_dir,$db_name,$cluster_id,$phylip_alignment_file,$id_map,$tmp_process_dir,\%db_distance_pair,$distance_pair_dir);      
        
        #######(4) Filter Non-overlapping aligned sequences ######
        STEP_4:
        print "Check overlap between sequence pairs\n";
        my($alignment_type)=check_overlap($cluster_id,$alignment_file,$alignment_dir,\%pairwise_alignment,$distance_pair,$pairwise_alignment_dir,$tmp_process_dir,$alignment_overlap_threshold,$pair_identity_threshold,$pair_similarity_threshold);
        
        #### (5) Identify reciprocal smallest distance pair #####
        print "Predicting Reciprocal Match\n";
        my($reciprocal_pair,$duplicate_pair)=get_reciprocal_smallest_distance_pair($cluster_id,$distance_pair,$alignment_type,$tmp_process_dir);
        
        #### (6) Get distance threshold #####
        print "Getting Distance Threshold\n";
        my($distance_threshold)=get_distance_threshold($cluster_id,$reciprocal_pair,\%outgroup,$divergence_threshold);
        
        #### (7) Identify Ortholog Pairs ####
        print "Predicting Ortholog\n";
        my($ortholog_pair)=predict_ortholog_pair($cluster_id,$reciprocal_pair,$distance_threshold,\%outgroup);
        
        #### (8) Identify Inparalog Pairs #######
        print "Predicting InParalog\n";
        my($inparalog_pair,$duplicate_minus_inparalog_pair)=predict_inparalog_pair($cluster_id,$ortholog_pair,$duplicate_pair,$inparalog_divergence_threshold);
           
        REMOVE_DIR:
        remove_tree("$tmp_process_dir"); 
  
        
######### SUBROUTINES #######

######## CREATE MULTIPLE ALIGNMENT FOR HOMOLOG FAMILY ########

sub create_msa{

  my($db_dir)=(shift);
  my($db_name)=(shift);
  my($cluster_id)=(shift);
  my($cluster_aln_record)=(shift);
  my($cluster_gene)=(shift);
  my($tmp_process_dir)=(shift);
  
  my $list_gene_id='';
  my $list_genome_name='';
  my $count_new_seq=0;
  my $name_old_aln_file=undef;
  my $name_new_aln_file=undef;
  my $aln_file=undef;
  my %family_sequence=();
  my %aligned_sequence=();
  
  my @cluster_gene=split(' ',$cluster_gene);
  
  my %list_gene_in_cluster=();
  
  foreach my $gene_id(@cluster_gene){
  
      $gene_id=~/(\w+)(\|)(.+)/;
      my $id=$3;
      my $genome_name=$1;
      $list_gene_id=$list_gene_id.",'".$id."'";   
      $list_genome_name=$list_genome_name.",'".$genome_name."'";    
      
      $list_gene_in_cluster{$gene_id}=$gene_id;
               
  }
  $list_gene_id=~s/^\,//g;
  $list_genome_name=~s/^\,//g;
  
   ###### GET ALIGNMENT SEQUENCE FOR HOMOLOG FAMILY IF PRESENT IN DATABASE ####
   if($cluster_aln_record and scalar(@{$cluster_aln_record})>1){  
   
        open(OLD_ALN,">$tmp_process_dir/old_aln_$cluster_id.afa");    
        
        $name_old_aln_file="$tmp_process_dir/old_aln_$cluster_id.afa";
        
        foreach my $row(@{$cluster_aln_record}){            
            my @column=@{$row};
            my $gene_id=$column[0];
            chomp($gene_id);
            my $aln_seq=$column[5]; 
            my $aln_group_id=$column[3];
            chomp($aln_group_id);
            
            if(!defined($list_gene_in_cluster{$gene_id})){
               next;
            }

            
            
            $aligned_sequence{$aln_group_id}->{$gene_id}=$aln_seq;    

            $aln_seq=~s/[UJO]/X/gi;   
            
            print OLD_ALN ">$gene_id\n$aln_seq\n";      
        }        
        close OLD_ALN;
   }
  
   ######## GET SEQUENCES FOR HOMOLOG FAMILY FROM PROTEIN SEQUENCE DATABASE######
        my($get_sequence)="SELECT * from ProteinSequence WHERE protein_id IN ($list_gene_id) and genome_abbreviation IN ($list_genome_name)";

        my($family_sequence)=SQLiteDB::get_record($db_dir,$db_name,$get_sequence);
                
        if(scalar($family_sequence)>1){  
        
             open(NEW_FASTA_SEQ,">$tmp_process_dir/$cluster_id.fasta");              
                  
             foreach my $row(@{$family_sequence}){           
                     my @data_column=@{$row};     
                     my $genome_name=$data_column[2];
                     my $gene_id=$genome_name."|".$data_column[1];
                     my $seq=$data_column[5];                          
                     $family_sequence{$genome_name}->{$gene_id}=$seq;

                     if($genome_name=~/\~/){next;}
                     
                     if(!defined($list_gene_in_cluster{$gene_id})){
                        next;
                     }
                     
                     if(!defined($aligned_sequence{$cluster_id}->{$gene_id})){ 


                        $seq=~s/[UJO]/X/gi;  ###### replace bad amino acids from the sequence with X ####
                   
                        print NEW_FASTA_SEQ ">$gene_id\n$seq\n";                      
                        $count_new_seq++;
                     }               
             }
            close NEW_FASTA_SEQ;

        } 
        
        ########## PERFORM MULTIPLE SEQUENCE ALIGNMENT ########
        if((-s "$tmp_process_dir/$cluster_id.fasta") and ($count_new_seq>1)){
        
           system("kalign -in  $tmp_process_dir/$cluster_id.fasta -output $tmp_process_dir/new_aln_$cluster_id.afa -f fasta -q");
          
           $name_new_aln_file="$tmp_process_dir/new_aln_$cluster_id.afa";
           
        }elsif((-s "$tmp_process_dir/$cluster_id.fasta") and ($count_new_seq==1)){ 
             
           $name_new_aln_file="$tmp_process_dir/$cluster_id.fasta";
           
        }
        
        ####### PROFILE ALIGNMENT WITH EXISTING HOMOLOG ALIGNMENT #########
        if($cluster_aln_record and scalar(@{$cluster_aln_record})>=1 and $name_new_aln_file  and -s "$name_new_aln_file"){   
           
            system("muscle -profile -in1 $name_old_aln_file -in2 $name_new_aln_file -out $tmp_process_dir/$cluster_id.afa");   
             
            $aln_file="$tmp_process_dir/$cluster_id.afa";
            
        }else{    
           
             $aln_file=$name_new_aln_file;           
        }
        
   return(\%family_sequence,$aln_file);
}

##### PARSE FASTA TO PHYLIP FORMAT ######

sub convertFastaToPhylip {

    my($fasta_aln_file)=(shift);
    my($tmp_aln_file)=(shift);
    my %id_map=();

    my $in_fasta_aln=Bio::SeqIO->new(-file =>"$fasta_aln_file",-format => "fasta");

    open(TMP_ALN,">$tmp_aln_file");
    my $seq_counter="0001";
  
    while(my $aln_seq=$in_fasta_aln->next_seq()){
        my $seq_id=$aln_seq->id;           
        my $temp_id="SEQ_".$seq_counter;
        print TMP_ALN ">$temp_id\n",$aln_seq->seq(),"\n";
        $id_map{$temp_id}=$seq_id;
        $seq_counter++; 
    }
    close TMP_ALN;

    my $phylip_file=$tmp_aln_file;
       $phylip_file=~s/\.afa/\.phy/g;

    system("perl Fasta2Phylip.pl $tmp_aln_file $phylip_file");

    return($phylip_file,\%id_map);
}

##### CALCULATE PAIRWISE GENETIC DISTANCE USING PHYLIP ####

sub calculate_pairwise_distance{

    my($db_dir)=(shift);
    my($db_name)=(shift);
    my($cluster_id)=(shift);
    my($phylip_alignment)=(shift);
    my($id_map)=(shift);
    my($tmp_process_dir)=(shift);
    my($db_distance_pair)=(shift);
    my($distance_pair_dir)=(shift);
    
    ##### PHYLIP PARAM FILE ####
    
    open(PARAMS,">$tmp_process_dir/param");
    print PARAMS "tmp_aln.phy\n";
    print PARAMS "Y";
    close PARAMS;
    
    #### EXECUTE PROTEIN DISTANCE #####
    chdir("$tmp_process_dir"); 
    
    system("protdist < param > protdist.log");
    
    chdir($Bin);
    
    my($distance_pair)=ParseDistanceMatrix("$tmp_process_dir/outfile",$id_map,$cluster_id,$db_distance_pair,$distance_pair_dir);
    
    return($distance_pair);
}

####### Parse Distance Matrix #######################

sub ParseDistanceMatrix {

    my($outfile)=(shift);    
    my(%temp_id_map)=%{(shift)};
    my($cluster_id)=(shift);   
    my(%db_distance_pair)=%{(shift)};
    my($distance_pair_dir)=(shift);
    
    #my @pair=();    
    
    my %distance_pair=();

    my $dist = Bio::Tools::Phylo::Phylip::ProtDist->new(-file=> $outfile,                  
	                                                    -format=> 'phylip',
                                                        -program=>"ProtDist");  
    my $matrix=$dist->next_matrix;
    
    open(DIST,">$distance_pair_dir/distance_$cluster_id.txt");

    foreach my $temp_idA(keys %temp_id_map){

           my $seq_idA=$temp_id_map{$temp_idA}; 
           $seq_idA=~/(\w+)(\|)(.+)/;
           my $genomeA=$1;
           my $gene_idA=$3;                 

        foreach my $temp_idB(keys %temp_id_map){

              my $seq_idB=$temp_id_map{$temp_idB};
              $seq_idB=~/(\w+)(\|)(.+)/;
              my $genomeB=$1;
              my $gene_idB=$3;   

                   if($distance_pair{$seq_idA}->{$seq_idB} and $distance_pair{$seq_idB}->{$seq_idA}){
                      next;
                   }

                   if($db_distance_pair{$seq_idA}->{$seq_idB} or $db_distance_pair{$seq_idB}->{$seq_idA}){
                      
                          my $distance_value=$db_distance_pair{$seq_idA}->{$seq_idB};

                          my $distance_line_F="$genomeA\t$gene_idA\t$genomeB\t$gene_idB\t$distance_value\t$cluster_id"; 
              
                          print DIST "$distance_line_F\n";

                          $distance_pair{$seq_idA}->{$seq_idB}=$distance_value;
                          $distance_pair{$seq_idB}->{$seq_idA}=$distance_value; 

                          next;
                   }           
                              
              my $distance_value = $matrix->get_entry($temp_idA,$temp_idB);

              if($distance_value<0){
                next;
              }
  
              my $distance_line_F="$genomeA\t$gene_idA\t$genomeB\t$gene_idB\t$distance_value\t$cluster_id"; 
              
              print DIST "$distance_line_F\n";
                                      
              $distance_pair{$seq_idA}->{$seq_idB}=$distance_value;
              $distance_pair{$seq_idB}->{$seq_idA}=$distance_value;                                                          
       }            
   } 
   
   close DIST;
   
   return(\%distance_pair);
}

####### Filter non-overlapping or little overlapping sequence ######

sub check_overlap{

    my($cluster_id)=(shift);
    my($alignment_file)=(shift);
    my($alignment_dir)=(shift);
    my(%pairwise_alignment)=%{(shift)};
    my(%distance_pair)=%{(shift)};
    my($pairwise_alignment_dir)=(shift);
    my($tmp_dir)=(shift);
    my($alignment_overlap_threshold)=(shift);
    my($pair_identity_threshold)=(shift);
    my($pair_similarity_threshold)=(shift);
    
    my %alignment_type=();
    my $seq_id_1='';
    my $seq_id_2='';
    my $identity=0.0;
    my $similarity=0.0;
    my $aln_region=0.0;
    my $evalue=1e-03;
    my $coverage=0.0;

    open(ALN_SEQ,">$alignment_dir/$cluster_id.aln.txt");

    open(PAIR_ALIGNMENT,">$pairwise_alignment_dir/$cluster_id.txt");
    
    my $in_fasta_aln=Bio::SeqIO->new(-file =>"$alignment_file",-format => "fasta");
    
    open(SEQ1,">$tmp_dir/QUERY.fasta");
    
    open(SEQ2,">$tmp_dir/SUBJECT.fasta");         
    
    while(my $aln_seq=$in_fasta_aln->next_seq()){
    
        $seq_id_1=$aln_seq->id;           
        my $seq_1=$aln_seq->seq;
        my $raw_seq1=$seq_1;
           $raw_seq1=~s/\-//g;
       
        $seq_id_1=~/(\w+)(\|)(.+)/;
        my $genome=$1;

        my $aln_len=length($seq_1);

        print ALN_SEQ "$seq_id_1\t$genome\tPROTEIN\t$cluster_id\t$aln_len\t$seq_1\n";      
        
        if(!defined($pairwise_alignment{$seq_id_1})){     
        
             print SEQ1 ">$seq_id_1\n$raw_seq1\n";
             
             print SEQ2 ">$seq_id_1\n$raw_seq1\n";
                      
        }else{
             print SEQ2 ">$seq_id_1\n$raw_seq1\n";
        }
        
    }    
    close SEQ1;
    close SEQ2;
    close ALN_SEQ;
    
    system("blastp -query $tmp_dir/QUERY.fasta -subject $tmp_dir/SUBJECT.fasta -evalue $evalue -out $tmp_dir/BLAST_OUTPUT.txt -outfmt '6 qseqid qlen sseqid slen pident ppos qcovs bitscore evalue qstart qend sstart send'");
    
    open(BLAST_RESULT,"$tmp_dir/BLAST_OUTPUT.txt");
    my @blast_result=<BLAST_RESULT>;

    foreach my $line(@blast_result){

         chomp($line);
         
         my @column=split(/\t/,$line);
         
         $seq_id_1=$column[0];
         
         $seq_id_2=$column[2];
         
         $identity=$column[4];
         
         $similarity=$column[5];
         
         $coverage=$column[6];

         if(!defined($pairwise_alignment{$seq_id_1}->{$seq_id_2}->{identity}) and !defined($pairwise_alignment{$seq_id_2}->{$seq_id_1}->{identity})){
           print PAIR_ALIGNMENT "$seq_id_1\t$seq_id_2\t$identity\t$similarity\t$cluster_id\n"; 
         }
         
         $pairwise_alignment{$seq_id_1}->{$seq_id_2}->{identity}=$identity;
         $pairwise_alignment{$seq_id_2}->{$seq_id_1}->{identity}=$identity;

         $pairwise_alignment{$seq_id_1}->{$seq_id_2}->{similarity}=$similarity;
         $pairwise_alignment{$seq_id_2}->{$seq_id_1}->{similarity}=$similarity;
         
    }
    close PAIR_ALIGNMENT;

    my $in_fasta_aln_1=Bio::SeqIO->new(-file =>"$alignment_file",-format => "fasta");

    while(my $aln_seq_1=$in_fasta_aln_1->next_seq()){
    
        my $id_1=$aln_seq_1->id;   

        my $in_fasta_aln_2=Bio::SeqIO->new(-file =>"$alignment_file",-format => "fasta");

        while(my $aln_seq_2=$in_fasta_aln_2->next_seq()){
    
             my $id_2=$aln_seq_2->id;   

             my $identity=0;
             my $similarity=0;

             unless(!defined($pairwise_alignment{$id_1}->{$id_2}->{identity}) or !defined($pairwise_alignment{$id_1}->{$id_2}->{similarity})){

                $identity=$pairwise_alignment{$id_1}->{$id_2}->{identity};
                $similarity=$pairwise_alignment{$id_1}->{$id_2}->{similarity};
             }

             if($identity>=$pair_identity_threshold and $similarity>=$pair_similarity_threshold){
                  
                    $alignment_type{$id_1}->{$id_2}=1;
                    $alignment_type{$id_2}->{$id_1}=1;
                    
             }else{
                    $alignment_type{$id_1}->{$id_2}=0;
                    $alignment_type{$id_2}->{$id_1}=0;
             }  
        }
    }        
    
    #foreach my $id_1(keys %pairwise_alignment){
    
    #     my %list_id_2=%{$pairwise_alignment{$id_1}};
         
    #     foreach my $id_2(keys %list_id_2){
          
    #              my $identity=$pairwise_alignment{$id_1}->{$id_2}->{identity};
    #              my $similarity=$pairwise_alignment{$id_1}->{$id_2}->{similarity};
              
    #              if($identity>=$pair_identity_threshold and $similarity>=$pair_similarity_threshold){
                  
    #                $alignment_type{$id_1}->{$id_2}=1;
    #                $alignment_type{$id_2}->{$id_1}=1;
                    
    #              }else{
    #                $alignment_type{$id_1}->{$id_2}=0;
    #                $alignment_type{$id_2}->{$id_1}=0;
    #              }  
    #     }
    #}

    return(\%alignment_type); 
}
 
###### get reciprocal smallest distance pair #####

sub get_reciprocal_smallest_distance_pair{

    my($cluster_id)=(shift);
    my(%distance_pair)=%{(shift)};
    my(%alignment_type)=%{(shift)};
    my($tmp_process_dir)=(shift);
    
    my %minimum_dist_pair=();
    
    my %duplicate_pair=();
    
    my %reciprocal_pair=();
    
    foreach my $idA(keys %distance_pair){
    
          my %protein_idB=%{$distance_pair{$idA}};
          
          $idA=~/(\w+)(\|)(.+)/;
          
          my $taxaA=$1;
          
          my %smallest_distance=();
          
          foreach my $idB(sort{$protein_idB{$a}<=>$protein_idB{$b}} keys %protein_idB){
          
               $idB=~/(\w+)(\|)(.+)/;
          
               my $taxaB=$1;  
               
               my $distance=$protein_idB{$idB};  

               if($alignment_type{$idA}->{$idB}==0 or !defined($alignment_type{$idA}->{$idB})){
                  next;
               }
               
               if($distance<0){
                  next;
               }
               
               #### duplicate protein pair ####
               if($taxaA eq $taxaB){               
                  $duplicate_pair{$idA}->{$idB}=$distance;
                  next;
               }
               
               #### non-overlapping pair ###
               if($distance<0){
                 next;
               }
               
               if(!defined($smallest_distance{$taxaB}) or $smallest_distance{$taxaB}==$distance){
                   $minimum_dist_pair{$idA}->{$idB}=$distance;
                   
                   $smallest_distance{$taxaB}=$distance;
               }
          }
    } 
    
    ######## get reciprocal closest pair #####
    
    foreach my $idA(keys %minimum_dist_pair){
    
           my %min_idB=%{$minimum_dist_pair{$idA}};
           
           foreach my $idB(keys %min_idB){
           
               if($minimum_dist_pair{$idB}->{$idA}){
               
                   my $reciprocal_distance=$min_idB{$idB};
               
                   $reciprocal_pair{$idA}->{$idB}=$reciprocal_distance;
                   
                   $reciprocal_pair{$idB}->{$idA}=$reciprocal_distance;
               }
           }      
    }
    
   return(\%reciprocal_pair,\%duplicate_pair); 
} 

####### GET DISTANCE THRESHOLD #####

sub get_distance_threshold{

    my($cluster_id)=(shift);
    my(%reciprocal_pair)=%{(shift)};
    my(%outgroup)=%{(shift)};
    my($divergence_threshold)=(shift);
    
    my %distance_threshold=();
    
    foreach my $idA(keys %reciprocal_pair){
    
         my %idB_list=%{$reciprocal_pair{$idA}};
         
        foreach my $idB(sort{$idB_list{$a}<=>$idB_list{$b}} keys %idB_list){
         
               $idB=~/(\w+)(\|)(.+)/;              
               my $taxaB=$1;
               my $distance_min=$idB_list{$idB};
               
               if($outgroup{$taxaB} and !defined($distance_threshold{$idA}->{minimum_threshold})){ 
                            
                   $distance_threshold{$idA}->{minimum_threshold}=$distance_min;
               }          
         }
         
         foreach my $idB(sort{$idB_list{$b}<=>$idB_list{$a}} keys %idB_list){
         
               $idB=~/(\w+)(\|)(.+)/;              
               my $taxaB=$1;
               my $distance_max=$idB_list{$idB};
                
               if($outgroup{$taxaB} and !defined($distance_threshold{$idA}->{maximum_threshold})){  
                          
                   $distance_threshold{$idA}->{maximum_threshold}=$distance_max;
               }          
         }
         
         #### IF NO OUT-GROUP PRESENT SET THRESHOLD TO USER_DEFINED VALUE ######
         if(!defined($distance_threshold{$idA}->{minimum_threshold}) and !defined($distance_threshold{$idA}->{maximum_threshold})){      
              $distance_threshold{$idA}->{minimum_threshold}=$divergence_threshold;
              $distance_threshold{$idA}->{maximum_threshold}=$divergence_threshold;          
         }
         ######## Check if maximum and minimum out-group distance is higher than user threshold, if yes set value to user- threshold ####
         if($distance_threshold{$idA}->{minimum_threshold}>$divergence_threshold){ 
              $distance_threshold{$idA}->{minimum_threshold}=$divergence_threshold;      
         }
         if($distance_threshold{$idA}->{maximum_threshold}>$divergence_threshold){
              $distance_threshold{$idA}->{maximum_threshold}=$divergence_threshold;
         }
    }

   return(\%distance_threshold);
}

###### PREDICT ORTHOLOG #######

sub predict_ortholog_pair{

    my($cluster_id)=(shift);
    my(%reciprocal_pair)=%{(shift)};
    my(%distance_threshold)=%{(shift)};
    my(%outgroup)=%{(shift)};

    my %ortholog_pair=();
    
    open(ORTHO,">$ortholog_pair_dir/$cluster_id.txt");
    
    foreach my $idA(keys %reciprocal_pair){
    
         $idA=~/(\w+)(\|)(.+)/;
          
         my $taxaA=$1;
         my $seq_idA=$3;
    
         my %idB_list=%{$reciprocal_pair{$idA}};
         
         my $min_cut_off=$distance_threshold{$idA}->{minimum_threshold};
         my $max_cut_off=$distance_threshold{$idA}->{maximum_threshold};
         
         if($outgroup{$taxaA}){
            $min_cut_off=$max_cut_off;
         }
         
         foreach my $idB(keys %idB_list){
         
            $idB=~/(\w+)(\|)(.+)/;              
            my $taxaB=$1;
            my $seq_idB=$3;
            
            my $distance=$idB_list{$idB};
            
            if(!defined($outgroup{$taxaB})){
                if($distance<=$min_cut_off){               
                   $ortholog_pair{$idA}->{$idB}=$distance;
                   print ORTHO "$taxaA\t$seq_idA\t$taxaB\t$seq_idB\t$distance\t$cluster_id\n";
                }
            }else{
               if($distance<=$max_cut_off){               
                 $ortholog_pair{$idA}->{$idB}=$distance;
                  print ORTHO "$taxaA\t$seq_idA\t$taxaB\t$seq_idB\t$distance\t$cluster_id\n";
               }
            }      
         }
    } 
    close ORTHO;
    
    return(\%ortholog_pair);    
}

#### Predict InParalog Pairs #####

sub predict_inparalog_pair{

    my($cluster_id)=(shift);
    my(%ortholog_pair)=%{(shift)};
    my(%duplicate_pair)=%{(shift)};
    my($inparalog_divergence_threshold)=(shift);

    my %inparalog_pair=();
    
    my %duplicate_minus_inparalog_pair=();
    
    open(INPARA,">$inparalog_pair_dir/$cluster_id.txt");
    
    foreach my $idA(keys %duplicate_pair){
    
        my $min_ortho_distA='';

          $idA=~/(\w+)(\|)(.+)/;
          
         my $taxaA=$1;
         my $seq_idA=$3;
        
        #### Minimum ortholog distance for Protein_A ####
        if($ortholog_pair{$idA}){
          my %ortholog_idB=%{$ortholog_pair{$idA}};
        
          foreach my $ortho_idB(sort{$ortholog_idB{$a}<=>$ortholog_idB{$b}} keys %ortholog_idB){  
            $min_ortho_distA=$ortholog_idB{$ortho_idB};
            last;
          }
        }
        
        my %dup_idB_list=%{$duplicate_pair{$idA}};
        
        foreach my $dup_idB(keys %dup_idB_list){
        
            #if($dup_idB eq $idA){next;}

            $dup_idB=~/(\w+)(\|)(.+)/;
          
            my $taxaB=$1;
            my $seq_idB=$3;
        
            #### Minimum ortholog distance for Protein_B ####
            my $min_ortho_distB='';
        
            if($ortholog_pair{$dup_idB}){
               my %ortholog_idB=%{$ortholog_pair{$dup_idB}};
        
                foreach my $ortho_idB(sort{$ortholog_idB{$a}<=>$ortholog_idB{$b}} keys %ortholog_idB){  
                   $min_ortho_distB=$ortholog_idB{$ortho_idB};
                   last;
                }
            }  
            
            if($min_ortho_distA eq '' and $min_ortho_distB eq ''){
               $min_ortho_distA=$inparalog_divergence_threshold;
               $min_ortho_distB=$inparalog_divergence_threshold;
               
            }elsif($min_ortho_distA ne '' and $min_ortho_distB eq ''){
            
                 $min_ortho_distB=$min_ortho_distA;
                 
            }elsif($min_ortho_distA eq '' and $min_ortho_distB ne ''){ 
            
                 $min_ortho_distA=$min_ortho_distB;   
            }  
            
            my $dup_distance=$dup_idB_list{$dup_idB};
            
            if($dup_distance<=$min_ortho_distA and $dup_distance<=$min_ortho_distB){
               $inparalog_pair{$idA}->{$dup_idB}=$dup_distance;               
               print INPARA "$taxaA\t$seq_idA\t$taxaB\t$seq_idB\t$dup_distance\t$min_ortho_distA\t$cluster_id\n";        
            
            }else{
                $duplicate_minus_inparalog_pair{$idA}->{$dup_idB}=$dup_distance;
            }       
        }
    }
    
 close INPARA;
 
 return(\%inparalog_pair,\%duplicate_minus_inparalog_pair);
}

