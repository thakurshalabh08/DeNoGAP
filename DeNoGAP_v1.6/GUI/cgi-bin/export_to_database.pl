#!/usr/bin/perl -w
###### ABOUT: This Script export gene information for selected gene id and stores in database ############
###### AUTHOR:Shalabh Thakur###################################################################


use strict;
use warnings;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Env;
use FindBin qw($Bin);
use DBI;

print "Content-type:text/html\n\n";

my $db_dir=param('db_dir');
my $cdb_name=param('c_database');
my $hdb_name=param('h_database');
my $odb_name=param('o_database');
my $locus_id=param('locus_tag');
my $protein_id=param('protein_id');
my $dbxref_id=param('dbxref');
my $species_name=param('species_name');
my $genome_abbreviation=param('abbreviation');
my $genome_id=param('genome_id');
my $genome_type=param('genome_type');
my $genome_length=param('genome_length');
my $index_genome=param('index_genome');
my $feature_start=param('feature_start');
my $feature_end=param('feature_end');
my $strand=param('strand');
my $frame=param('frame');
my $gene_name=param('gene_name');
my $feature_type=param('feature_type');
my $aa_len=param('protein_length');
my $nuc_len=param('nucleotide_length');
my $product_description=param('product_description');
my $comment=param('comments');
my $homolog_group=param('homolog_group');
my $ortholog_group=param('ortholog_group');
my $hmm_group=param('hmm_group');
my $go_id_list=param('list_go_id_name');
my $go_category_list=param('list_go_category_name');
my $go_description_list=param('list_go_desc_name');
my $pfam_id_list=param('list_pfam_id_tag');
my $pfam_name_list=param('list_pfam_name_tag');
my $pfam_desc_list=param('list_pfam_desc_tag');
my $pfam_start_list=param('list_pfam_start_tag');
my $pfam_end_list=param('list_pfam_end_tag');
my $pfam_significance_list=param('list_pfam_significance');
my $interpro_id_list=param('list_interpro_id_tag');
my $interpro_name_list=param('list_interpro_name_tag');
my $pathway_id_list=param('list_pathway_id_tag');
my $pathway_name_list=param('list_pathway_name_tag');
my $signalp_start=param('signalp_start_tag');
my $signalp_end=param('signalp_end_tag');
my $tmhmm_start_list=param('list_tmhmm_start_tag');
my $tmhmm_end_list=param('list_tmhmm_end_tag');
my $phobius_domain_list=param('list_phobius_domain_tag');
my $phobius_description_list=param('list_phobius_description_tag');
my $phobius_start_list=param('list_phobius_start_tag');
my $phobius_end_list=param('list_phobius_end_tag');

my @pfam_id_tag_list=split(";",$pfam_id_list);
my @pfam_name_tag_list=split(";",$pfam_name_list);
my @pfam_desc_tag_list=split(";",$pfam_desc_list);
my @pfam_start_tag_list=split(";",$pfam_start_list);
my @pfam_end_tag_list=split(";",$pfam_end_list);
my @pfam_sig_tag_list=split(";",$pfam_significance_list);

my @go_id_tag=split(";",$go_id_list);
my @go_category_tag=split(";",$go_category_list);
my @go_description_tag=split(";",$go_description_list);

my @interpro_id_tag=split(";",$interpro_id_list);
my @interpro_name_tag=split(";",$interpro_name_list);

my @pathway_id_tag=split(";",$pathway_id_list);
my @pathway_name_tag=split(";",$pathway_name_list);

my @tmhmm_start_tag=split(";",$tmhmm_start_list);
my @tmhmm_end_tag=split(";",$tmhmm_end_list);

my @phobius_domain_list=split(";",$phobius_domain_list);
my @phobius_description_list=split(";",$phobius_description_list);
my @phobius_start_list=split(";",$phobius_start_list);
my @phobius_end_list=split(";",$phobius_end_list);


######## set general annotation for gene id #######

#my $feature_stmt="SELECT * from GeneFeature WHERE feature_id='$locus_id' and genome_name='$species_name'";
#my ($feature_data)=get_record($db_dir,$cdb_name,$feature_stmt);

my $del_old_record="Delete from GeneFeature where feature_id='$locus_id' and genome_name='$species_name'";
my ($delete_data)=set_record($db_dir,$cdb_name,$del_old_record);

my $update_new_record="INSERT into GeneFeature (feature_id, feature_type, protein_id, dbxref, genome_id, genome_type, genome_name, genome_length, feature_start, feature_end, nuc_length, aa_length, strand, frame, index_in_genome, description , comment) VALUES ('$locus_id','$feature_type','$protein_id','$dbxref_id','$genome_id','$genome_type','$species_name','$genome_length','$feature_start','$feature_end','$nuc_len','$aa_len','$strand','$frame','$index_genome','$product_description','$comment')";
my ($add_data)=set_record($db_dir,$cdb_name,$update_new_record);


##### set domain annotation for protein id ######

for(my $i=0;$i<scalar(@pfam_id_tag_list);$i++){

    my $id_tag_name=$pfam_id_tag_list[$i];
    my $name_tag=$pfam_name_tag_list[$i];
    my $desc_tag=$pfam_desc_tag_list[$i];
    my $start_tag=$pfam_start_tag_list[$i];
    my $end_tag=$pfam_end_tag_list[$i];
    my $significance_tag=$pfam_sig_tag_list[$i];

    my $pfam_id_value=param($id_tag_name);
    my $pfam_name_value=param($name_tag);
    my $pfam_desc_value=param($desc_tag);
    my $pfam_start_value=param($start_tag);
    my $pfam_end_value=param($end_tag);
    my $significance_value=param($significance_tag);

my $del_old_domain_record="Delete from DomainAnnotation where protein_id='$protein_id' and domain_id='$pfam_id_value'";
my ($delete_domain)=set_record($db_dir,$cdb_name,$del_old_domain_record);

  if($pfam_name_value ne '' and $pfam_desc_value ne '' and $pfam_start_value ne '' and $pfam_end_value ne '' and $significance_value ne ''){
    my $update_domain_record="INSERT into DomainAnnotation (protein_id,genome_name,seq_len,domain_id,domain_name,domain_start,domain_end,significance_value,description) Values ('$protein_id','$genome_abbreviation','$aa_len','$pfam_id_value','$pfam_name_value','$pfam_start_value','$pfam_end_value','$significance_value','$pfam_desc_value')";
    my ($update_domain)=set_record($db_dir,$cdb_name,$update_domain_record);
  }
}

#### set go annotation for protein id ####

for(my $j=0;$j<scalar(@go_id_tag);$j++){

     my $go_id_tag=$go_id_tag[$j];
     my $go_category_tag=$go_category_tag[$j];
     my $go_description_tag=$go_description_tag[$j];

     my $go_id_value=param($go_id_tag);
     my $go_category_value=param($go_category_tag);
     my $go_desc_value=param($go_description_tag);

     my $del_old_go_record="Delete from GOAnnotation where protein_id='$protein_id' and go_id='$go_id_value'";
     my ($delete_go)=set_record($db_dir,$cdb_name,$del_old_go_record);

     if($go_category_value ne '' and $go_desc_value ne ''){
         my $update_go_record="INSERT into GOAnnotation (protein_id,genome_name,go_id,go_category,go_description) Values ('$protein_id','$genome_abbreviation','$go_id_value','$go_category_value','$go_desc_value')";
         my ($update_go)=set_record($db_dir,$cdb_name,$update_go_record);
     }
}

#### set interpro annotation for protein id ####

for(my $j=0;$j<scalar(@interpro_id_tag);$j++){

     my $interpro_id_tag=$interpro_id_tag[$j];
     my $interpro_name_tag=$interpro_name_tag[$j];

     my $interpro_id_value=param($interpro_id_tag);
     my $interpro_name_value=param($interpro_name_tag);

     my $del_old_interpro_record="Delete from InterProAnnotation where protein_id='$protein_id' and interpro_id='$interpro_id_value'";
     my ($delete_interpro)=set_record($db_dir,$cdb_name,$del_old_interpro_record);

     if($interpro_name_value ne ''){
        my $update_interpro_record="INSERT into InterProAnnotation (protein_id,genome_name,interpro_id,interpro_name) Values ('$protein_id','$genome_abbreviation','$interpro_id_value','$interpro_name_value')";
        my ($update_go)=set_record($db_dir,$cdb_name,$update_interpro_record);
     }
}

##### set pathway annotation for protein id #####

for(my $j=0;$j<scalar(@pathway_id_tag);$j++){

     my $pathway_id_tag=$pathway_id_tag[$j];
     my $pathway_name_tag=$pathway_name_tag[$j];

     my $pathway_id_value=param($pathway_id_tag);
     my $pathway_name_value=param($pathway_name_tag);

     my $del_old_pathway_record="Delete from PathwayAnnotation where protein_id='$protein_id' and pathway_id='$pathway_id_value'";
     my ($delete_pathway)=set_record($db_dir,$cdb_name,$del_old_pathway_record);

     if($pathway_name_value ne ''){
         my $update_pathway_record="INSERT into PathwayAnnotation (protein_id,genome_name,pathway_id,pathway_name) Values ('$protein_id','$genome_abbreviation','$pathway_id_value','$pathway_name_value')";
         my ($pathway_go)=set_record($db_dir,$cdb_name,$update_pathway_record);
     }
}

#### set Signalp annotation ###

my $del_old_signalp_record="Delete from SignalPAnnotation where protein_id='$protein_id'";
my ($delete_signalp)=set_record($db_dir,$cdb_name,$del_old_signalp_record);

if($signalp_start!=0 and $signalp_end!=0){

    my $update_signalp_record="INSERT into SignalPAnnotation (protein_id,genome_name,domain_start,domain_end) Values ('$protein_id','$genome_abbreviation','$signalp_start','$signalp_end')";
    my ($signalp_go)=set_record($db_dir,$cdb_name,$update_signalp_record);
}


#### set tmhmm annotation #####

my $del_old_tmhmm_record="Delete from TMHMMAnnotation where protein_id='$protein_id'";
my ($delete_tmhmm)=set_record($db_dir,$cdb_name,$del_old_tmhmm_record);

for(my $j=0;$j<scalar(@tmhmm_start_tag);$j++){

     my $tmhmm_start_tag=$tmhmm_start_tag[$j];
     my $tmhmm_end_tag=$tmhmm_end_tag[$j];

     my $tmhmm_start_value=param($tmhmm_start_tag);
     my $tmhmm_end_value=param($tmhmm_end_tag);

     if($tmhmm_start_value!=0 and $tmhmm_end_value!=0){

         my $update_tmhmm_record="INSERT into TMHMMAnnotation (protein_id,genome_name,domain_start,domain_end) Values ('$protein_id','$genome_abbreviation','$tmhmm_start_value','$tmhmm_end_value')";
         my ($tmhmm_go)=set_record($db_dir,$cdb_name,$update_tmhmm_record);
     }
}

#### set phobius annotation #####

for(my $j=0;$j<scalar(@phobius_domain_list);$j++){

     my $phobius_domain_tag=$phobius_domain_list[$j];
     my $phobius_description_tag=$phobius_description_list[$j];
     my $phobius_start_tag=$phobius_start_list[$j];
     my $phobius_end_tag=$phobius_end_list[$j];

     my $phobius_domain_value=param($phobius_domain_tag);
     my $phobius_description_value=param($phobius_description_tag);
     my $phobius_start_value=param($phobius_start_tag);
     my $phobius_end_value=param($phobius_end_tag);

     if($phobius_start_value!=0 and $phobius_end_value!=0){

         my $del_old_phobius_record="Delete from PhobiusAnnotation where protein_id='$protein_id' and domain_name='$phobius_domain_value'";
         my ($delete_phobius)=set_record($db_dir,$cdb_name,$del_old_phobius_record);

         my $update_phobius_record="INSERT into PhobiusAnnotation (protein_id,genome_name,domain_name,domain_description,domain_start,domain_end) Values ('$protein_id','$genome_abbreviation','$phobius_domain_value','$phobius_description_value','$phobius_start_value','$phobius_end_value')";
         my ($phobius_go)=set_record($db_dir,$cdb_name,$update_phobius_record);
     }else{

         my $del_old_phobius_record="Delete from PhobiusAnnotation where protein_id='$protein_id' and domain_name='$phobius_domain_value'";
         my ($delete_phobius)=set_record($db_dir,$cdb_name,$del_old_phobius_record);
     }
}

print qq*
<span
style="font-weight: bold;">
Annotation is updated for $locus_id. Go back and refresh Gene Information page to see updated information.
</span>
*;

#### get record from database ####
sub set_record {

   my $db_dir=(shift);
   my $db_name=(shift);
   my $sql_stmt=(shift);

   chdir($db_dir);

   my $dbh=DBI->connect("dbi:SQLite:dbname=$db_name","","",{RaiseError =>1}) or die $DBI::errstr; 
     
   my $stmt=$dbh->prepare("$sql_stmt");
   $stmt->execute();

   $stmt->finish;
   $dbh->disconnect();
   chdir($Bin);
}

#### get record from database ####
sub get_record {

   my $db_dir=(shift);
   my $db_name=(shift);
   my $sql_stmt=(shift);

   my @row_data=();

   chdir($db_dir);
   my $dbh=DBI->connect("dbi:SQLite:dbname=$db_name","","",{RaiseError =>1}) or die $DBI::errstr;   
   my $stmt=$dbh->prepare("$sql_stmt");
   $stmt->execute();

   while(my @row=$stmt->fetchrow_array()){       
      push(@row_data,\@row);
   }

   $stmt->finish;
   $dbh->disconnect();
   chdir($Bin);
   return(\@row_data);
}


