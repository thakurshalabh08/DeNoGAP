#!/usr/bin/perl

use strict;
use warnings;
use Env;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Range;
use Tie::File;


my $blast_alignment_file=(shift);
my $feature_dir=(shift);
my $cds_dir=(shift);
my $protein_dir=(shift);
my $genome_dir=(shift);
my $project_dir=(shift);
my $verify_cds=(shift);
my $verify_protein=(shift);
my $verify_feature=(shift);
my $verify_genbank=(shift);
my $verify_gff=(shift);
my $evalue_cutoff=(shift);
my $identity_cutoff=(shift);
my $query_coverage_cutoff=(shift);
my $seq_len_cutoff=(shift);

chomp($seq_len_cutoff);

mkdir($project_dir);
mkdir($verify_cds);
mkdir($verify_protein);
mkdir($verify_feature);
mkdir($verify_genbank);
mkdir($verify_gff);

##### PARSE BLAST ALIGNMENT FILE #####
my $report = Bio::SearchIO->new(-format => 'blast', -file =>"$blast_alignment_file"); 
    
    my %gene_annotation=();

    while (my $result = $report->next_result) {

       my $query=$result->query_name;
       $query=~s/(_1)$//g;

       while( my $hit = $result->next_hit ) {   
          
          while( my $hsp = $hit->next_hsp ) {
          
              my $hit_desc=$hit->description; 
              my $evalue=$hit->significance; 
              my $identity=$hsp->percent_identity;
              my $qstart=$hsp->start('query');
              my $qend=$hsp->end('query');
              my $sstart=$hsp->start('hit'); 
              my $send=$hsp->end('hit');
              my $qlen=$result->query_length;
              my $slen=$hit->length;

              my $qcoverage=(($qend-$qstart)+1/$qlen)*100;
              my $scoverage=(($send-$sstart)+1/$slen)*100;
              
              $hit_desc=~s/(OS=)(.+)//g;

              if(($identity>=$identity_cutoff) and ($qcoverage>=$query_coverage_cutoff) and ($evalue <= $evalue_cutoff)){
                  $gene_annotation{$query}=$hit_desc;
              }
          }
       }
    }


############# FEATURE ##################

opendir(FEATURE,"$feature_dir/MULTIPLE");
my @feature_file=readdir(FEATURE);

my %verified_seq_id=();

foreach my $feature_file(@feature_file){

        if($feature_file=~/^\.+$/){next;}

        my $genome_name=$feature_file;
           $genome_name=~s/(\.)(.+)//g;

        my @feature=();
        tie @feature, 'Tie::File', "$feature_dir/MULTIPLE/$feature_file";   
        
        my @feature_s=();
        tie @feature_s, 'Tie::File', "$feature_dir/SINGLE/$genome_name.feature.single.txt";  
        
        push(@feature,@feature_s);


        print "Verifying $genome_name sequences\n";

        open(V_FEATURE,">$verify_feature/$genome_name.txt");    

        open(V_GFF,">$verify_gff/$genome_name.gff");

        my $gff_line=1;

        foreach my $feature_line(@feature){

              if($feature_line=~/\#/){
                next;
              }
  
              my @feature_column=split("\t",$feature_line);
              
              my $seq_id=$feature_column[0];
              my $seq_len=$feature_column[10];

          my $feature_id=$feature_column[0];
          my $feature_type=$feature_column[1];
          my $protein_id=$feature_column[2];
          my $dbxref=$feature_column[3];
          my $genome_id=$feature_column[4];
          my $genome_type=$feature_column[5];
          my $genome_name=$feature_column[6];
          my $genome_length=$feature_column[7];
          my $feature_start=$feature_column[8];
          my $feature_end=$feature_column[9];
          my $nuc_length=$feature_column[10];
          my $aa_length=$feature_column[11];
          my $strand=$feature_column[12];
          my $frame=$feature_column[13];
          my $index_on_genome=$feature_column[14];
          my $description=$feature_column[15];

              if($gene_annotation{$seq_id}){
                 
                  my $gene_description=$gene_annotation{$seq_id};
 
                  pop(@feature_column);
                  push(@feature_column,$gene_description);

                  my $feature_line=join("\t",@feature_column);

                  print V_FEATURE "$feature_line\n"; 

                  $verified_seq_id{$seq_id}=$genome_name;

                  if($gff_line==1){
                   print V_GFF "##gff-version 3\n";
                   print V_GFF "# organismn $genome_name\n";
                  }

                  print V_GFF "$genome_id\tGenBank\tregion\t1\t$genome_length\t.\t$strand\t$frame\tID=$genome_id;Name=$genome_id;Note=$genome_name\n";
                  print V_GFF "$genome_id\tGenBank\tgene\t$feature_start\t$feature_end\t.\t$strand\t$frame\tID=$feature_id;Name=$feature_id\n";
                  print V_GFF "$genome_id\tGenBank\tmRNA\t$feature_start\t$feature_end\t.\t$strand\t$frame\tID=$feature_id".".t01".";Parent=$feature_id\n";
                  print V_GFF "$genome_id\tGenBank\tCDS\t$feature_start\t$feature_end\t.\t$strand\t$frame\tID=$feature_id".".p01".";Parent=$feature_id".".t01".";Name=$feature_id;codon_start=$frame;product=$gene_description;protein_id=$protein_id;translation=length.$aa_length\n";
                  print V_GFF "$genome_id\tGenBank\texon\t$feature_start\t$feature_end\t.\t$strand\t$frame\tParent=$feature_id".".t01\n";
                
                  $gff_line++;

              }elsif($seq_len>=$seq_len_cutoff){
                
                  pop(@feature_column);
                  push(@feature_column,"hypothetical protein");

                  my $feature_line=join("\t",@feature_column);

                  print V_FEATURE "$feature_line\n"; 

                  $verified_seq_id{$seq_id}=$genome_name;

                  if($gff_line==1){
                   print V_GFF "##gff-version 3\n";
                   print V_GFF "# organismn $genome_name\n";
                  }

                  print V_GFF "$genome_id\tGenBank\tregion\t1\t$genome_length\t.\t$strand\t$frame\tID=$genome_id;Name=$genome_id;Note=$genome_name\n";
                  print V_GFF "$genome_id\tGenBank\tgene\t$feature_start\t$feature_end\t.\t$strand\t$frame\tID=$feature_id;Name=$feature_id\n";
                  print V_GFF "$genome_id\tGenBank\tmRNA\t$feature_start\t$feature_end\t.\t$strand\t$frame\tID=$feature_id".".t01".";Parent=$feature_id\n";
                  print V_GFF "$genome_id\tGenBank\tCDS\t$feature_start\t$feature_end\t.\t$strand\t$frame\tID=$feature_id".".p01".";Parent=$feature_id".".t01".";Name=$feature_id;codon_start=$frame;product=hypothetical protein;protein_id=$protein_id;translation=length.$aa_length\n";
                  print V_GFF "$genome_id\tGenBank\texon\t$feature_start\t$feature_end\t.\t$strand\t$frame\tParent=$feature_id".".t01\n";
                 
                  $gff_line++;
              }         
        }


       my $seqio_obj_genome = Bio::SeqIO->new(-file => "$genome_dir/$genome_name.fasta", -format => "fasta");

       while(my $seq_obj= $seqio_obj_genome->next_seq){
    
          my $seq=$seq_obj->seq;
          my $seq_id=$seq_obj->display_id;
          
          if($seq_id=~/\|/){
             $seq_id=~s/\|$//g;
             my @contig=split(/\|/,$seq_id);
             $seq_id=pop(@contig);   
          }

           print V_GFF ">$seq_id\n$seq\n";

       }
        close V_FEATURE;
        close V_GFF;
}


#### PROTEIN ####

opendir(PROTEIN,"$protein_dir/MULTIPLE");
my @protein_file=readdir(PROTEIN);

foreach my $protein_file(@protein_file){

     if($protein_file=~/^\.+$/){next;}

    my $seqio_obj = Bio::SeqIO->new(-file => "$protein_dir/MULTIPLE/$protein_file", -format => "fasta");

    my $genome_name=$protein_file;
       $genome_name=~s/(\.)(.+)//g;

    open(V_PROTEIN,">$verify_protein/$genome_name.fasta");   

    open(V_GFF,">>$verify_gff/$genome_name.gff"); 

    while(my $seq_obj= $seqio_obj->next_seq){
    
       my $seq=$seq_obj->seq;
       my $seq_id=$seq_obj->display_id;
       my $seq_desc=$seq_obj->desc;

       my $gene_id=$seq_id;
       $gene_id=~s/(_1)$//g;

       $seq_desc=~/(chromosome=)(.+)/;

       my $seq_desc_temp=$1.$2;

       if($gene_annotation{$gene_id} or length($seq)>=$seq_len_cutoff){

          if($gene_annotation{$gene_id}){
            $seq_desc="product=$gene_annotation{$gene_id}".$seq_desc_temp;
          }

          print V_PROTEIN ">$seq_id $seq_desc\n$seq\n";

          print V_GFF ">$seq_id".".p01"."\n".$seq."\n"; 
       }   
   }
   
   my $seqio_obj_s = Bio::SeqIO->new(-file => "$protein_dir/SINGLE/$genome_name.aa.single.fasta", -format => "fasta");
   
    while(my $seq_obj= $seqio_obj_s->next_seq){
    
       my $seq=$seq_obj->seq;
       my $seq_id=$seq_obj->display_id;
       my $seq_desc=$seq_obj->desc;

       my $gene_id=$seq_id;
       $gene_id=~s/(_1)$//g;

       $seq_desc=~/(chromosome=)(.+)/;

       my $seq_desc_temp=$1.$2;

       if($gene_annotation{$gene_id} or length($seq)>=$seq_len_cutoff){

          if($gene_annotation{$gene_id}){
            $seq_desc="product=$gene_annotation{$gene_id}".$seq_desc_temp;
          }

          print V_PROTEIN ">$seq_id $seq_desc\n$seq\n";

          print V_GFF ">$seq_id".".p01"."\n".$seq."\n";
       }   
   }

   close V_PROTEIN;
   close V_GFF;
}

#### CDS ###########

opendir(CDS,"$cds_dir/MULTIPLE");
my @cds_file=readdir(CDS);

foreach my $cds_file(@cds_file){

     if($cds_file=~/^\.+$/){next;}

    my $seqio_obj = Bio::SeqIO->new(-file => "$cds_dir/MULTIPLE/$cds_file", -format => "fasta");

    my $genome_name=$cds_file;
       $genome_name=~s/(\.)(.+)//g;

    open(V_CDS,">$verify_cds/$genome_name.fasta");    

    while(my $seq_obj= $seqio_obj->next_seq){
    
       my $seq=$seq_obj->seq;
       my $seq_id=$seq_obj->display_id;
       my $seq_desc=$seq_obj->desc;

       my $aa_length=length($seq)/3;

       $seq_desc=~/(chromosome=)(.+)/;

       my $seq_desc_temp=$1.$2;

       if($gene_annotation{$seq_id} or $aa_length>=$seq_len_cutoff){

          if($gene_annotation{$seq_id}){
            $seq_desc="product=$gene_annotation{$seq_id}".$seq_desc_temp;
          }

          print V_CDS ">$seq_id $seq_desc\n$seq\n";
       }   
    }
    
    my $seqio_obj_s = Bio::SeqIO->new(-file => "$cds_dir/SINGLE/$genome_name.single.fasta", -format => "fasta");
    
    while(my $seq_obj= $seqio_obj_s->next_seq){
    
       my $seq=$seq_obj->seq;
       my $seq_id=$seq_obj->display_id;
       my $seq_desc=$seq_obj->desc;

       my $aa_length=length($seq)/3;

       $seq_desc=~/(chromosome=)(.+)/;

       my $seq_desc_temp=$1.$2;

       if($gene_annotation{$seq_id} or $aa_length>=$seq_len_cutoff){

          if($gene_annotation{$seq_id}){
            $seq_desc="product=$gene_annotation{$seq_id}".$seq_desc_temp;
          }

          print V_CDS ">$seq_id $seq_desc\n$seq\n";
       }   
    }

   close V_CDS;
}

##### GenBank Files #####

opendir(V_FEATURE,"$verify_feature");
my @verify_feature_file=readdir(V_FEATURE);

foreach my $feature_file(@verify_feature_file){

        if($feature_file=~/^\.+$/){next;}

        my @feature=();
        tie @feature, 'Tie::File', "$verify_feature/$feature_file";   

        my $genome_name=$feature_file;
           $genome_name=~s/(\.)(.+)//g;
           
        my $genome_file="$genome_dir/$genome_name.fasta";   

        print "Creating GenBank file for $genome_name\n";

        system("perl CreateGBK.pl $genome_name $genome_file $verify_cds/$genome_name.fasta $verify_protein/$genome_name.fasta $verify_feature/$genome_name.txt $verify_genbank/$genome_name.gbk"); 

}


