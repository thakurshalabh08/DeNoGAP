##### Module to Parse GeneBank File#####
##### Author: Shalabh Thakur ################
##### Date: 6-AUG-2013 #####################

#!/usr/bin/perl
package ParseGBK;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocationI;
use Bio::Range;
use File::Basename;
use vars qw(@ISA @EXPORT @EXPORT_OK $db_dir $db_name);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(gbk_genome print_genome_seq_file print_coding_seq_file print_feature_file);


sub gbk_genome {

    my($gbk_file)=(shift);
    my($genome_name)=(shift);

    my %genome_seq=();
    my %coding_seq=();
    my %protein_seq=();
    my %feature=();

    my $cds_id_count="00001";

    my $organism_name='';
    my $genome_length=0;

    my  $infile = new Bio::SeqIO(-file =>"$gbk_file", '-format' => 'genbank');

    print "$gbk_file\n\n";

    while(my $seq=$infile->next_seq){

         my $genome_id=$seq->display_id;    #### contig or chromosome ID ####
         
         print "$genome_id\n";

         my $cds_index=1;

          for my $feature_obj($seq->get_SeqFeatures){

             if(!defined($feature_obj->entire_seq->seq)){
               next;
             }

             $genome_seq{$genome_name}->{$genome_id}=$feature_obj->entire_seq->seq;  #### contig or chromosome sequence #####

            if(!defined( $genome_seq{$genome_name}->{$genome_id})){
               print "$genome_id do not have contig sequence\n";
               exit;
            }
      
            if($feature_obj->primary_tag eq "region" or $feature_obj->primary_tag eq "source"){

                 $genome_length=$feature_obj->location->end;   #### genome length ####                
                 ######## get all tags  for source ####### 
                 for my $tag ($feature_obj->get_all_tags){ 
                     if($tag eq "organism"){ #### organism name ####
		         for my $value ($feature_obj->get_tag_values($tag)){
                             $organism_name=$value;                            
                         }
		     }
                 }

            }elsif($feature_obj->primary_tag eq "gene"){

                for my $tag ($feature_obj->get_all_tags){

                    if($tag eq "locus_tag"){ 
                         for my $value ($feature_obj->get_tag_values($tag)){
                            $feature{$genome_id}->{$cds_index}->{locus_tag}=$value;                            
                         }
                    }
                }

            }elsif($feature_obj->primary_tag eq "CDS"){ ##### Coding Sequence features ####

                $feature{$genome_id}->{$cds_index}->{organism_name}=$organism_name;
                $feature{$genome_id}->{$cds_index}->{genome_length}=$genome_length;
                $feature{$genome_id}->{$cds_index}->{genome_id}=$genome_id;
                $feature{$genome_id}->{$cds_index}->{cds_start}=$feature_obj->location->start; 
                $feature{$genome_id}->{$cds_index}->{cds_end}=$feature_obj->location->end; 
                $feature{$genome_id}->{$cds_index}->{strand}=$feature_obj->location->strand; 
                $feature{$genome_id}->{$cds_index}->{frame}=1;    
                $feature{$genome_id}->{$cds_index}->{db_xref}='';       
                                       
                ##### get all tags for CDS ############
                for my $tag ($feature_obj->get_all_tags){

                    if($tag eq "locus_tag"){ #### get locus_tag id ####
                         for my $value ($feature_obj->get_tag_values($tag)){
                            $feature{$genome_id}->{$cds_index}->{locus_tag}=$value;
                         }
                    }
              
                    if($tag eq "protein_id"){ #### get protein id #####
                         for my $value ($feature_obj->get_tag_values($tag)){
                            $feature{$genome_id}->{$cds_index}->{protein_id}=$value;
                         }
                    }

                    if($tag eq "product"){  #### get product description ####
                         for my $value ($feature_obj->get_tag_values($tag)){
                            $feature{$genome_id}->{$cds_index}->{product_description}=$value;
                         }
                    }
                  
                    if($tag eq "codon_start"){ #### get frame of translation ####
                        for my $value ($feature_obj->get_tag_values($tag)){
                            $feature{$genome_id}->{$cds_index}->{frame}=$value;
                        }
                    }

                    if($tag eq "protein_length"){ #### get frame of translation ####
                        for my $value ($feature_obj->get_tag_values($tag)){
                            $feature{$genome_id}->{$cds_index}->{protein_length}=$value;
                        }
                    }

                    if($tag eq "pseudo"){ ### get information if gene is pseudo ###
                        $feature{$genome_id}->{$cds_index}->{product_description}= "pseudogene";
                    }

                    if($tag eq "db_xref"){
                        for my $value ($feature_obj->get_tag_values($tag)){
                            $feature{$genome_id}->{$cds_index}->{db_xref}=$feature{$genome_id}->{$cds_index}->{db_xref}."|".$value;
                        }
                    }if($tag eq "translation"){
                        for my $value ($feature_obj->get_tag_values($tag)){
                            my $locus_id='';
                            if($feature{$genome_id}->{$cds_index}->{locus_tag}){

                               $locus_id=$feature{$genome_id}->{$cds_index}->{locus_tag};

                            }elsif($feature{$genome_id}->{$cds_index}->{protein_id}){

                               $locus_id=$feature{$genome_id}->{$cds_index}->{protein_id};

                            }elsif($feature{$genome_id}->{$cds_index}->{db_xref}){

                              $locus_id=$feature{$genome_id}->{$cds_index}->{db_xref};

                            }
                            
                            $protein_seq{$genome_name}->{$genome_id}->{$locus_id}=$value;
                            $feature{$genome_id}->{$cds_index}->{protein_sequence}=$value;
                            $feature{$genome_id}->{$cds_index}->{protein_length}= length($protein_seq{$genome_name}->{$genome_id}->{$locus_id});
                        }
                    }
                } #### End of get all tags

                ##### check protein product annotation ##########

                if(!defined($feature{$genome_id}->{$cds_index}->{product_description})){
                   $feature{$genome_id}->{$cds_index}->{product_description}="undefined protein sequence";
                }

                ##### check if locus tag or protein id found ####

                if($feature{$genome_id}->{$cds_index}->{protein_id} and !defined($feature{$genome_id}->{$cds_index}->{locus_tag}) and !defined($feature{$genome_id}->{$cds_index}->{db_xref})){
 
                   $feature{$genome_id}->{$cds_index}->{locus_tag}=$feature{$genome_id}->{$cds_index}->{protein_id};
                   $feature{$genome_id}->{$cds_index}->{db_xref}=$feature{$genome_id}->{$cds_index}->{protein_id};

                }elsif($feature{$genome_id}->{$cds_index}->{locus_tag} and !defined($feature{$genome_id}->{$cds_index}->{protein_id}) and !defined($feature{$genome_id}->{$cds_index}->{db_xref})){

                   $feature{$genome_id}->{$cds_index}->{protein_id}=$feature{$genome_id}->{$cds_index}->{locus_tag};
                   $feature{$genome_id}->{$cds_index}->{db_xref}=$feature{$genome_id}->{$cds_index}->{locus_tag};

                }elsif($feature{$genome_id}->{$cds_index}->{db_xref} and !defined($feature{$genome_id}->{$cds_index}->{protein_id}) and !defined($feature{$genome_id}->{$cds_index}->{locus_tag})){

                   $feature{$genome_id}->{$cds_index}->{protein_id}=$feature{$genome_id}->{$cds_index}->{db_xref};
                   $feature{$genome_id}->{$cds_index}->{locus_tag}=$feature{$genome_id}->{$cds_index}->{db_xref};

                }elsif(!defined($feature{$genome_id}->{$cds_index}->{locus_tag}) and !defined($feature{$genome_id}->{$cds_index}->{protein_id}) and !defined($feature{$genome_id}->{$cds_index}->{db_xref})){

                   $feature{$genome_id}->{$cds_index}->{locus_tag}=$genome_name."_".$cds_id_count;
                   $feature{$genome_id}->{$cds_index}->{protein_id}=$genome_name."_".$cds_id_count; 
      
                   $cds_id_count++;
                }

                my $locus_id=$feature{$genome_id}->{$cds_index}->{locus_tag};                
                
                $feature{$genome_id}->{$cds_index}->{nuc_length}=($feature{$genome_id}->{$cds_index}->{cds_end}-$feature{$genome_id}->{$cds_index}->{cds_start})+1;

                if($feature{$genome_id}->{$cds_index}->{product_description} eq "pseudogene"){
                
                   $feature{$genome_id}->{$cds_index}->{protein_length}= 0;
                   
                }else{  

                   $feature{$genome_id}->{$cds_index}->{nuc_length}=($feature{$genome_id}->{$cds_index}->{cds_end}-$feature{$genome_id}->{$cds_index}->{cds_start}) + 1;  
                  
                        my $start_pos=$feature{$genome_id}->{$cds_index}->{cds_start};
                        my $end_pos=$feature{$genome_id}->{$cds_index}->{cds_end};

                        my $len=($end_pos-$start_pos)+1;
                        
                        $coding_seq{$genome_name}->{$genome_id}->{$locus_id}=substr($genome_seq{$genome_name}->{$genome_id},$start_pos-1,$len);
                        $coding_seq{$genome_name}->{$genome_id}->{$locus_id}=reverse($coding_seq{$genome_name}->{$genome_id}->{$locus_id});
                        $coding_seq{$genome_name}->{$genome_id}->{$locus_id}=~tr/ATGCatgc/TACGtacg/;
                }
                $cds_index++;                   
        } #### END OF CDS FEATURE BLOCK  
                               
      } ### END OF FEATURE BLOCK ####
    
    }  #### END OF SEQUENCE BLOCK

  return(\%genome_seq,\%coding_seq,\%protein_seq,\%feature);
}


sub find_stop{

      my($coding_seq)=(shift);

      my $i=0;
      my $stop=0;
      my $j=1;
      my $seq='';

      my $num_codon=int(length($coding_seq)/3);

      while($j<=$num_codon){

           my $code=substr($coding_seq,$i,3);
           
           if($code=~/(TAA|TAG|TGA|taa|tag|tga)/ and $j!=$num_codon){
              $stop=1;
              last;
           }          
           $seq=$seq.$code;
           $i=$i+3;
           $j++;
      }
    return($stop,$seq);
}

######### print genome sequene file #########

sub print_genome_seq_file{

  my(%genome_seq)=%{(shift)};
  my($genome_name)=(shift);
  my($genome_dir)=(shift);

  my $genome_file_name= $genome_dir."/".$genome_name.".fasta";
 
  open(GENOME_FILE,">$genome_file_name");

  foreach my $genome_name(keys %genome_seq){
     
       my %genome_id=%{$genome_seq{$genome_name}};

       foreach my $genome_id(keys %genome_id){          
          print GENOME_FILE ">$genome_id\n$genome_id{$genome_id}\n";
       }
  } 
  close GENOME_FILE;
}

sub print_coding_seq_file{

  my(%coding_seq)=%{(shift)};
  my(%protein_seq)=%{(shift)};
  my($genome_name)=(shift);
  my($cds_dir)=(shift);
  my($protein_dir)=(shift);

    my $coding_seq_file=$cds_dir."/".$genome_name.".fasta";
    my $protein_seq_file=$protein_dir."/".$genome_name.".fasta";

    open(CDS_FILE,">$coding_seq_file");
    open(PROTEIN_FILE,">$protein_seq_file");
  
    foreach my $genome_name(keys %coding_seq){
     
       my %genome_id=%{$coding_seq{$genome_name}};

       foreach my $genome_id(keys %genome_id){

           my %cds=%{$genome_id{$genome_id}};

           foreach my $locus_id(keys %cds){
              print CDS_FILE ">$locus_id\n$cds{$locus_id}\n";
              print PROTEIN_FILE ">$locus_id\n$protein_seq{$genome_name}->{$genome_id}->{$locus_id}\n";            
           }
       }
    }

   close CDS_FILE;
   close PROTEIN_FILE;
}

sub print_feature_file{

  my(%feature)=%{(shift)};
  my($genome_name)=(shift);
  my($feature_dir)=(shift);
  my(%genome_seq)=%{(shift)};
  my($gff_dir)=(shift);

  my $feature_file=$feature_dir."/".$genome_name.".txt";

  my $gff_file=$gff_dir."/".$genome_name.".gff";
 
  open(FEATURE,">$feature_file");

  open(GFF,">$gff_file");
      
   my $gff_line=1;

  my $strand='';

  my %gff_protein_seq=();
      
  print FEATURE "#feature_id\tfeature_type\tprotein_id\tdbxref\tgenome_id\tgenome_type\tgenome_name\tgenome_length\tfeature_start\tfeature_end\tnuc_lenth\taa_length\tstrand\tframe\tindex_on_genome\tdescription\n";

  foreach my $genome_id(keys %feature){

     my %cds_index=%{$feature{$genome_id}};
    
     foreach my $index(sort{$a<=>$b} keys %cds_index){

         my %cds_feature=%{$cds_index{$index}};

         if(!defined($cds_feature{genome_id}) or !defined($cds_feature{cds_start})){
           next;
         }

         print FEATURE "$cds_feature{locus_tag}\t";
         print FEATURE "CDS\t";
         print FEATURE "$cds_feature{protein_id}\t";
         print FEATURE "$cds_feature{db_xref}\t";
         print FEATURE "$cds_feature{genome_id}\t";
         print FEATURE "contig\t";
         print FEATURE "$cds_feature{organism_name}\t";
         print FEATURE "$cds_feature{genome_length}\t";  
         print FEATURE "$cds_feature{cds_start}\t";
         print FEATURE "$cds_feature{cds_end}\t";
         print FEATURE "$cds_feature{nuc_length}\t";
         print FEATURE "$cds_feature{protein_length}\t";
         if($cds_feature{strand} eq 1){
         print FEATURE "+\t";
         $strand="+";
         }elsif($cds_feature{strand} eq -1){
         print FEATURE "-\t";
         $strand="-";
         }
         print FEATURE "$cds_feature{frame}\t";
         print FEATURE "$index\t";
         print FEATURE "$cds_feature{product_description}\n";   

         if($gff_line==1){
            print GFF "##gff-version 3\n";
            print GFF "# organismn $cds_feature{organism_name}\n";
         } 

         print GFF "$cds_feature{genome_id}\tGenBank\tregion\t1\t$cds_feature{genome_length}\t.\t$strand\t$cds_feature{frame}\tID=$cds_feature{genome_id};Name=$cds_feature{genome_id};Note=$cds_feature{organism_name}\n";
         print GFF "$cds_feature{genome_id}\tGenBank\tgene\t$cds_feature{cds_start}\t$cds_feature{cds_end}\t.\t$strand\t$cds_feature{frame}\tID=$cds_feature{locus_tag};Name=$cds_feature{locus_tag}\n";
         print GFF "$cds_feature{genome_id}\tGenBank\tmRNA\t$cds_feature{cds_start}\t$cds_feature{cds_end}\t.\t$strand\t$cds_feature{frame}\tID=$cds_feature{locus_tag}".".t01".";Parent=$cds_feature{locus_tag}\n";
         print GFF "$cds_feature{genome_id}\tGenBank\tCDS\t$cds_feature{cds_start}\t$cds_feature{cds_end}\t.\t$strand\t$cds_feature{frame}\tID=$cds_feature{locus_tag}".".p01".";Parent=$cds_feature{locus_tag}".".t01".";Name=$cds_feature{locus_tag};codon_start=$cds_feature{frame};product=$cds_feature{product_description};protein_id=$cds_feature{protein_id};translation=length.$cds_feature{protein_length}\n";
         print GFF "$cds_feature{genome_id}\tGenBank\texon\t$cds_feature{cds_start}\t$cds_feature{cds_end}\t.\t$strand\t$cds_feature{frame}\tParent=$cds_feature{locus_tag}".".t01\n";
        
         $gff_protein_seq{$cds_feature{locus_tag}}=$cds_feature{protein_sequence};
         $gff_line++;
     }
  }
   
  print GFF "##FASTA\n";

  foreach my $genome_name(keys %genome_seq){
     
       my %genome_id=%{$genome_seq{$genome_name}};

       foreach my $genome_id(keys %genome_id){  
        
          print GFF ">$genome_id\n$genome_id{$genome_id}\n";
       }
  } 

  foreach my $protein_id(keys %gff_protein_seq){

       print GFF ">$protein_id".".p01"."\n".$gff_protein_seq{$protein_id}."\n";
  }

  close FEATURE;
  close GFF;
}
