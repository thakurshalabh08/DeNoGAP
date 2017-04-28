#!/usr/bin/perl -w
###### ABOUT: This Script find gene information for selected gene id ############
###### AUTHOR:Shalabh Thakur###################################################################


use strict;
use warnings;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Env;
use FindBin qw($Bin);
use DBI;

print "Content-type:text/html\n\n";

my $project_name=param('project_name');
my $db_dir=param('db_dir');
my $db_name=param('db_name');
my $hdb_name=param('hdb_name');
my $odb_name=param('odb_name');
my $report_dir=param('report_dir');
my $gene_id=param('gene_detail_radio');
my $genome_code=param('genome_name');

my %gene_annotation=();

###### get feature details #####
$gene_annotation{locus_tag}=$gene_id;

### get full name of the species /strain ###
my $species_stmt="SELECT * from OrganismInfo WHERE abbreviation='$genome_code'";
my ($species_data)=get_record($db_dir,$db_name,$species_stmt);

 if(scalar(@{$species_data})>0){
      foreach my $row(@{$species_data}){
          shift(@{$row});          
          $gene_annotation{genome_full_name}=shift(@{$row});
          $gene_annotation{species}=shift(@{$row});
          $gene_annotation{species_type}=shift(@{$row});
          $gene_annotation{species_abbreviation}=shift(@{$row});
          shift(@{$row});
      }
 }

### get Taxonomy data ####
my $taxonomy_stmt="SELECT * from Taxonomy WHERE species='$gene_annotation{species}'";
my ($taxonomy_data)=get_record($db_dir,$db_name,$taxonomy_stmt);

### get Genomic Feature ###
my $feature_stmt="SELECT * from GeneFeature WHERE feature_id='$gene_id' and genome_name='$gene_annotation{genome_full_name}'";
my ($feature_data)=get_record($db_dir,$db_name,$feature_stmt);

 if(scalar(@{$feature_data})>0){

      foreach my $row(@{$feature_data}){
          shift(@{$row}); 
          shift(@{$row});            
          $gene_annotation{feature_type}=shift(@{$row});
          $gene_annotation{protein_id}=shift(@{$row});
          $gene_annotation{dbxref_id}=shift(@{$row});
          $gene_annotation{genome_id}=shift(@{$row});
          $gene_annotation{genome_type}=shift(@{$row});
          shift(@{$row}); 
          $gene_annotation{genome_length}=shift(@{$row});
          $gene_annotation{feature_start}=shift(@{$row});
          $gene_annotation{feature_end}=shift(@{$row});
          $gene_annotation{nuc_len}=shift(@{$row});
          $gene_annotation{aa_len}=shift(@{$row}); 
          $gene_annotation{strand}=shift(@{$row});         
          $gene_annotation{index_on_genome}=shift(@{$row}); 
          $gene_annotation{frame}=shift(@{$row});          
          $gene_annotation{product_description}=shift(@{$row});  
          $gene_annotation{gene_name}="Not available";
          $gene_annotation{comment}=shift(@{$row});      
      }
 }

##### get comparative genomics details #####
my $hmm_group_stmt="SELECT * from GenetoSuperFamily WHERE gene_id='$gene_id' and genome_name='$gene_annotation{species_abbreviation}'";
my ($hmm_group_data)=get_record($db_dir,$hdb_name,$hmm_group_stmt);

 if(scalar(@{$hmm_group_data})>0){
      foreach my $row(@{$hmm_group_data}){
         shift(@{$row});
         shift(@{$row});
         shift(@{$row});
         $gene_annotation{hmm_group}=shift(@{$row});
         $gene_annotation{homolog_group}=shift(@{$row});
      }
 }

my $ortholog_group_stmt="SELECT * from MapGeneIdtoGeneFamily WHERE gene_id='$gene_id' and species_abbreviation='$gene_annotation{species_abbreviation}'";
my ($ortho_group_data)=get_record($db_dir,$odb_name,$ortholog_group_stmt);

if(scalar(@{$ortho_group_data})>0){
      foreach my $row(@{$ortho_group_data}){
         shift(@{$row});
         $gene_annotation{ortholog_group}=shift(@{$row});
         $gene_annotation{ortholog_group}=~s/\://g;
      }
}

#### get annotation details #######
my $pfam_stmt="SELECT * from DomainAnnotation where protein_id='$gene_annotation{protein_id}' and genome_name='$gene_annotation{species_abbreviation}'";
my ($pfam_data)=get_record($db_dir,$db_name,$pfam_stmt);

if(scalar(@{$pfam_data})>0){     
  $gene_annotation{pfam_domain}=$pfam_data;
}else{
  $gene_annotation{pfam_domain}="No annotation";
}

my $go_stmt="SELECT * from GOAnnotation where protein_id='$gene_annotation{protein_id}' and genome_name='$gene_annotation{species_abbreviation}'";
my ($go_data)=get_record($db_dir,$db_name,$go_stmt);

if(scalar(@{$go_data})>0){     
  $gene_annotation{go_annotation}=$go_data;
}else{
  $gene_annotation{go_annotation}="No annotation";
}

my $ipr_stmt="SELECT * from InterProAnnotation where protein_id='$gene_annotation{protein_id}' and genome_name='$gene_annotation{species_abbreviation}'";
my ($ipr_data)=get_record($db_dir,$db_name,$ipr_stmt);

if(scalar(@{$ipr_data})>0){     
  $gene_annotation{interpro_annotation}=$ipr_data;
}else{
  $gene_annotation{interpro_annotation}="No annotation";
}

my $pathway_stmt="SELECT * from PathwayAnnotation where protein_id='$gene_annotation{protein_id}' and genome_name='$gene_annotation{species_abbreviation}'";
my ($pathway_data)=get_record($db_dir,$db_name,$pathway_stmt);

if(scalar(@{$pathway_data})>0){     
  $gene_annotation{pathway_annotation}=$pathway_data;
}else{
  $gene_annotation{pathway_annotation}="No annotation";
}

my $signalp_stmt="SELECT * from SignalPAnnotation where protein_id='$gene_annotation{protein_id}' and genome_name='$gene_annotation{species_abbreviation}'";
my ($signalp_data)=get_record($db_dir,$db_name,$signalp_stmt);

if(scalar(@{$signalp_data})>0){     
  $gene_annotation{signalp_annotation}=$signalp_data;
}else{
  $gene_annotation{signalp_annotation}="No annotation";
}

my $tmhmm_stmt="SELECT * from TMHMMAnnotation where protein_id='$gene_annotation{protein_id}' and genome_name='$gene_annotation{species_abbreviation}'";
my ($tmhmm_data)=get_record($db_dir,$db_name,$tmhmm_stmt);

if(scalar(@{$tmhmm_data})>0){     
  $gene_annotation{tmhmm_annotation}=$tmhmm_data;
}else{
  $gene_annotation{tmhmm_annotation}="No annotation";
}

my $phobius_stmt="SELECT * from PhobiusAnnotation where protein_id='$gene_annotation{protein_id}' and genome_name='$gene_annotation{species_abbreviation}'";
my ($phobius_data)=get_record($db_dir,$db_name,$phobius_stmt);

if(scalar(@{$phobius_data})>0){     
  $gene_annotation{phobius_annotation}=$phobius_data;
}else{
  $gene_annotation{phobius_annotation}="No annotation";
}

### get sequence details ######


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



print qq*
<html>
<head>
<meta content="text/html; charset=ISO-8859-1"
http-equiv="content-type">
<title></title>
</head>
<body>
<form method="post" name="Gene_Information" action="edit_gene_annotation.pl" target="_blank"> <br>
<table style="text-align: left; width: 70%; height: 70%; margin-left: auto; margin-right: auto;" border="1" cellpadding="2" cellspacing="2">
<tbody>
<tr>
 <td colspan="2" rowspan="1" style="vertical-align: top; width: 305px;">
   <input name="Edit" value="Edit" type="submit">
 </td>
</tr>
<tr>
<td colspan="2" rowspan="1"
style="vertical-align: top; width: 305px;">
<div style="text-align: center;"><big style="font-weight: bold;"><big><big>GENE
: $gene_id</big></big></big><br>
<br>
</div>
</td>
</tr>
<tr>
<td colspan="2" rowspan="1"
style="width: 305px; background-color: rgb(204, 204, 204);"><big
style="font-weight: bold;"><big>Genomic Feature:</big></big><br>
</td>
</tr>
<tr>
<td colspan="2" rowspan="1" style="width: 305px;"></td>
</tr>
<tr>
<td style="width: 305px; font-weight: bold;">Locus Tag:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_id<br>
<input name="locus_tag" id="locus_tag" type="hidden" value="$gene_id">
<input name="db_dir" id="db_dir" type="hidden" value="$db_dir">
<input name="c_database" id="c_database" type="hidden" value="$db_name">
<input name="h_database" id="h_database" type="hidden" value="$hdb_name">
<input name="o_database" id="o_database" type="hidden" value="$odb_name">
</td>
</tr>
<tr>
<td style="width: 305px; font-weight: bold;">Protein Id:<br>
</td>
<td style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{protein_id}<br>
<input name="protein_id" id="protein_id" type="hidden" value="$gene_annotation{protein_id}">
</td>
</tr>
<tr>
<td style="width: 305px; font-weight: bold;">External database Id:<br>
</td>
<td style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{dbxref_id}<br>
<input name="dbxref" id="dbxref" type="hidden" value="$gene_annotation{dbxref_id}">
</td>
</tr>
<tr>
<td style="width: 305px; font-weight: bold;">Species / Strain
Name:<br>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{genome_full_name}<br>
<input name="species_name" id="species_name" type="hidden" value="$gene_annotation{genome_full_name}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Species Type:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{species_type}<br>
<input name="species_type" id="species_type" type="hidden" value="$gene_annotation{species_type}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Species Abbreviation:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{species_abbreviation}<br>
<input name="abbreviation" id="abbreviation" type="hidden" value="$gene_annotation{species_abbreviation}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Genome ID: <br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{genome_id}<br>
<input name="genome_id" id="genome_id" type="hidden" value="$gene_annotation{genome_id}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Genome Type: <br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{genome_type}<br>
<input name="genome_type" id="genome_type" type="hidden" value="$gene_annotation{genome_type}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Genome Length: <br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{genome_length}<br>
<input name="genome_length" id="genome_length" type="hidden" value="$gene_annotation{genome_length}">
</td>
</tr>
<tr>
<tr>
<td style="font-weight: bold;">Index on Genome: <br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{index_on_genome}<br>
<input name="index_genome" id="index_genome" type="hidden" value="$gene_annotation{index_on_genome}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Genomic Location:<br>
</td>
<td style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{feature_start} : $gene_annotation{feature_end} ($gene_annotation{strand}$gene_annotation{frame})<br>
<input name="feature_start" id="feature_start" type="hidden" value="$gene_annotation{feature_start}">
<input name="feature_end" id="feature_end" type="hidden" value="$gene_annotation{feature_end}">
<input name="strand" id="strand" type="hidden" value="$gene_annotation{strand}">
<input name="frame" id="frame" type="hidden" value="$gene_annotation{frame}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Gene Name:<br>
</td>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{gene_name}<br>
<input name="gene_name" id="gene_name" type="hidden" value="$gene_annotation{gene_name}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Feature Type:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{feature_type}<br>
<input name="feature_type" id="feature_type" type="hidden" value="$gene_annotation{feature_type}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Protein Length:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{aa_len}<br>
<input name="aa_len" id="aa_len" type="hidden" value="$gene_annotation{aa_len}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Nucleotide Length:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{nuc_len}<br>
<input name="nuc_len" id="nuc_len" type="hidden" value="$gene_annotation{nuc_len}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Product Description:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{product_description}<br>
<input name="product_description" id="product_description" type="hidden" value="$gene_annotation{product_description}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Comments:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{comment}<br>
<input name="comment" id="comment" type="hidden" value="$gene_annotation{comment}">
</td>
</tr>
<tr>
<td colspan="2" rowspan="1"
style="background-color: rgb(204, 204, 204);"><big
style="font-weight: bold;"><big>Comparative Genomics Information</big>:</big><br>
</td>
</tr>
<tr>
</tr>
<tr>
<td style="font-weight: bold;">Homolog Group:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{homolog_group}<br>
<input name="homolog_group" id="homolog_group" type="hidden" value="$gene_annotation{homolog_group}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">Ortholog Group:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{ortholog_group}<br>
<input name="ortholog_group" id="ortholog_group" type="hidden" value="$gene_annotation{ortholog_group}">
</td>
</tr>
<tr>
<td style="font-weight: bold;">HMM Model Group:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{hmm_group}<br>
<input name="hmm_group" id="hmm_group" type="hidden" value="$gene_annotation{hmm_group}">
</td>
</tr>
<tr>
<td colspan="2" rowspan="1"
style="background-color: rgb(204, 204, 204);"><big
style="font-weight: bold;"><big>Annotation</big>:</big><br>
</td>
</tr>
<tr valign="center">
<td style="font-weight: bold;">GO Annotation:<br>
</td>
<td>
<br>
<br>
<table style="vertical-align: center; width: 100%; height: 15%; margin-left: auto; margin-right: auto;" border="0" cellpadding="2" cellspacing="2">
<tbody>
<tr><th> GO ID </th><th>GO Category</th><th>GO Description</th></tr>
*;

if($gene_annotation{go_annotation} ne "No annotation"){

      my $go_tag_count=1;

      my $go_id_list='';
      my $go_category_list='';
      my $go_description_list='';

      foreach my $row(@{$gene_annotation{go_annotation}}){
         shift(@{$row});
         shift(@{$row});
         my $go_id=shift(@{$row});
         my $go_category=shift(@{$row});
         my $go_description=shift(@{$row});

          print qq*
            <tr align="center">
              <td>$go_id              
              </td> 
              <td>$go_category
              </td> 
              <td> $go_description
              </td> 
            </tr>   
              
          *;

          $go_id_list=$go_id_list.";".$go_id;
          $go_category_list=$go_category_list.";".$go_category;
          $go_description_list=$go_description_list.";".$go_description;
      }

      print qq*
      <tr>
        <td>
           <input name="GO_ID" id="GO_ID" type="hidden" value="$go_id_list">
           <input name="GO_CATEGORY" id="GO_CATEGORY" type="hidden" value="$go_category_list">
           <input name="GO_DESC" id="GO_DESC" type="hidden" value="$go_description_list">
        </td>
      <tr>
      *;
}else{
   print qq*
      <tr align="center">
        <td>No annotation</td>
       </tr>     
       <tr>
        <td>
           <input name="GO_ID" id="GO_ID" type="hidden" value="No annotation">
           <input name="GO_CATEGORY" id="GO_CATEGORY" type="hidden" value="No annotation">
           <input name="GO_DESC" id="GO_DESC" type="hidden" value="No annotation">
        </td>
      <tr>    
      *;
}
print qq*
</tbody>
</table>
<br>
<br>
</td>
</tr>
<tr valign="center">
<td style="font-weight: bold;">PFam:<vr>
</td>
<td style="vertical-align: top; width: 70%; text-align: left;"><br><br>
<table style="vertical-align: center; width: 100%; height:15%; margin-left: auto; margin-right: auto;" border="0" cellpadding="5" cellspacing="5">
<tbody>
<tr><th> PFam ID</th><th>Pfam Name</th><th>Description</th><th>Start</th><th>End</th><th>Significance</th></tr>
*;

if($gene_annotation{pfam_domain} ne "No annotation"){

      my $pfam_tag_count=1;

      my $pfam_id_list='';
      my $pfam_name_list='';
      my $pfam_description_list='';
      my $pfam_start_list='';
      my $pfam_end_list='';
      my $pfam_significance_list='';

      foreach my $row(@{$gene_annotation{pfam_domain}}){
         shift(@{$row});
         shift(@{$row});
         shift(@{$row});
         my $pfam_id=shift(@{$row});
         my $pfam_name=shift(@{$row});
         my $pfam_start=shift(@{$row});
         my $pfam_end=shift(@{$row});
         my $significance=shift(@{$row});
         my $pfam_description=shift(@{$row});
         
         print qq*
            <tr align="center">
              <td>$pfam_id              
              </td> 
              <td>$pfam_name
              </td> 
              <td>$pfam_description
              </td>
              <td>$pfam_start
              </td>
              <td>$pfam_end              
              </td>
              <td>$significance              
              </td>
            </tr>                    
          *;

         $pfam_id_list=$pfam_id_list.";".$pfam_id;
         $pfam_name_list=$pfam_name_list.";".$pfam_name;
         $pfam_description_list=$pfam_description_list.";".$pfam_description;
         $pfam_start_list=$pfam_start_list.";".$pfam_start;
         $pfam_end_list=$pfam_end_list.";".$pfam_end;
         $pfam_significance_list=$pfam_significance_list.";".$significance;
      }
     
      print qq*
        <tr>
           <td>
              <input name="PFAM_ID" id="PFAM_ID" type="hidden" value="$pfam_id_list">
              <input name="PFAM_NAME" id="PFAM_NAME" type="hidden" value="$pfam_name_list">
              <input name="PFAM_DESC" id="PFAM_DESC" type="hidden" value="$pfam_description_list">
              <input name="PFAM_START" id="PFAM_START" type="hidden" value="$pfam_start_list">
              <input name="PFAM_END" id="PFAM_END" type="hidden" value="$pfam_end_list">
              <input name="PFAM_SIGNIFICANCE" id="PFAM_SIGNIFICANCE" type="hidden" value="$pfam_significance_list">
           </td>
        <tr>
      *;
}else{
   print qq*
      <tr align="center">
        <td>No annotation</td>
       </tr> 
        <tr>
           <td>
              <input name="PFAM_ID" id="PFAM_ID" type="hidden" value="No annotation">
              <input name="PFAM_NAME" id="PFAM_NAME" type="hidden" value="No annotationt">
              <input name="PFAM_DESC" id="PFAM_DESC" type="hidden" value="No annotation">
              <input name="PFAM_START" id="PFAM_START" type="hidden" value="No annotation">
              <input name="PFAM_END" id="PFAM_END" type="hidden" value="No annotation">
              <input name="PFAM_SIGNIFICANCE" id="PFAM_SIGNIFICANCE" type="hidden" value="No annotation">
           </td>
        <tr>        
   *;
}
print qq*
</tbody>
</table>
<br>
<br>
</td>
</tr>
<tr valign="center">
<td style="font-weight: bold;">InterPro:<br>
</td>
<td style="vertical-align: top; width: 70%; text-align: left;"><br><br>
<table style="vertical-align: center; width: 100%; height:15%; margin-left: auto; margin-right: auto;" border="0" cellpadding="5" cellspacing="5">
<tbody>
<tr><th>InterPro ID</th><th>Description</th></tr>
*;

if($gene_annotation{interpro_annotation} ne "No annotation"){
      my $prev_ipr_id='';

      my $interpro_tag_count=1;
      my $interpro_id_list='';
      my $interpro_name_list='';

      foreach my $row(@{$gene_annotation{interpro_annotation}}){
         shift(@{$row});
         shift(@{$row});
         my $interpro_id=shift(@{$row});
         my $interpro_name=shift(@{$row});
   
         if($interpro_id ne $prev_ipr_id){
            print qq*
              <tr align="center">
                <td>$interpro_id
                </td> 
                <td>$interpro_name
                </td>
              </tr>                    
            *;
           $prev_ipr_id=$interpro_id;

           $interpro_id_list=$interpro_id_list.";".$interpro_id;
           $interpro_name_list=$interpro_name_list.";".$interpro_name;
         }
      }
     
      print qq*
      <tr>
        <td>
           <input name="INTERPRO_ID" id="INTERPRO_ID" type="hidden" value="$interpro_id_list">
           <input name="INTERPRO_NAME" id="INTERPRO_NAME" type="hidden" value="$interpro_name_list">
        </td>
      <tr>
      *;
       
}else{
   print qq*
      <tr align="center">
        <td>No annotation</td>
        <td>
           <input name="INTERPRO_ID" id="INTERPRO_ID" type="hidden" value="No annotation">
           <input name="INTERPRO_NAME" id="INTERPRO_NAME" type="hidden" value="No annotation">
        </td>
      </tr>         
   *;
}
print qq*
</tbody>
</table>
<br>
<br>
</td>
</tr>
<tr valign="center">
<td style="font-weight: bold;">Pathway:<br>
</td>
<td style="vertical-align: top; width: 70%; text-align: left;"><br><br>
<table style="vertical-align: center; width: 100%; height:15%; margin-left: auto; margin-right: auto;" border="0" cellpadding="5" cellspacing="5">
<tbody>
<tr><th>Pathway ID</th><th>Description</th></tr>
*;

if($gene_annotation{pathway_annotation} ne "No annotation"){
      my $prev_path_id='';

      my $path_tag_count=1;

      my $pathway_id_list='';
      my $pathway_name_list='';

      foreach my $row(@{$gene_annotation{pathway_annotation}}){
         shift(@{$row});
         shift(@{$row});
         my $pathway_id=shift(@{$row});
         my $pathway_name=shift(@{$row});
   
         if($pathway_id ne $prev_path_id){
            print qq*
              <tr align="center">
                <td>$pathway_id                
                </td> 
                <td>$pathway_name
                </td>
              </tr>                    
            *;

           $prev_path_id=$pathway_id;

           $pathway_id_list=$pathway_id_list.";".$pathway_id;
           $pathway_name_list=$pathway_name_list.";".$pathway_name;
         }
      }
      
      print qq*
      <tr>
        <td>
           <input name="PATHWAY_ID" id="PATHWAY_ID" type="hidden" value="$pathway_id_list">
           <input name="PATHWAY_NAME" id="PATHWAY_NAME" type="hidden" value="$pathway_name_list">
        </td>
      <tr>
      *;
}else{
   print qq*
      <tr align="center">
        <td>No annotation</td>
        <td>
           <input name="PATHWAY_ID" id="PATHWAY_ID" type="hidden" value="No annotation">
           <input name="PATHWAY_NAME" id="PATHWAY_NAME" type="hidden" value="No annotation">
        </td>
       </tr>         
      *;
}
print qq*
</tbody>
</table>
<br>
<br>
</td>
</tr>

<tr valign="center">
<td style="font-weight: bold;">SignalP:<br>
</td>
<td style="vertical-align: top; width: 70%; text-align: left;"><br><br>
<table style="vertical-align: center; width: 100%; height:15%; margin-left: auto; margin-right: auto;" border="0" cellpadding="5" cellspacing="5">
<tbody>
<tr><th>Signal Peptide</th><th>Start</th><th>End</th></tr>
*;

if($gene_annotation{signalp_annotation} ne "No annotation"){

      my $signalp_tag_count=1;

      my $signalp_start_list='';
      my $signalp_end_list='';

      foreach my $row(@{$gene_annotation{signalp_annotation}}){
         shift(@{$row});
         shift(@{$row});
         my $start=shift(@{$row});
         my $end=shift(@{$row});

            print qq*
              <tr align="center">
                <td>Yes                
                </td>
                <td>$start               
                </td> 
                <td>$end
                </td>
              </tr>                    
            *;

        $signalp_start_list=$start;
        $signalp_end_list=$end;
      }
      
      print qq*
      <tr>
        <td>
           <input name="SIGNALP_PEPTIDE" id="SIGNALP_START" type="hidden" value="Yes">
           <input name="SIGNALP_START" id="SIGNALP_START" type="hidden" value="$signalp_start_list">
           <input name="SIGNALP_END" id="SIGNALP_END" type="hidden" value="$signalp_end_list">
        </td>
      <tr>
      *;
}else{
   print qq*
      <tr align="center">
                <td>No                
                </td>
                <td>0               
                </td> 
                <td>0
                </td>
        <td>
           <input name="SIGNALP_PEPTIDE" id="SIGNALP_START" type="hidden" value="No">
           <input name="SIGNALP_START" id="SIGNALP_START" type="hidden" value="0">
           <input name="SIGNALP_END" id="SIGNALP_END" type="hidden" value="0">
        </td>
       </tr>         
      *;
}
print qq*
</tbody>
</table>
<br>
<br>
</td>
</tr>

<tr valign="center">
<td style="font-weight: bold;">TMHMM:<br>
</td>
<td style="vertical-align: top; width: 70%; text-align: left;"><br><br>
<table style="vertical-align: center; width: 100%; height:15%; margin-left: auto; margin-right: auto;" border="0" cellpadding="5" cellspacing="5">
<tbody>
<tr><th>Domain Name</th><th>Start</th><th>End</th></tr>
*;

if($gene_annotation{tmhmm_annotation} ne "No annotation"){

      my $tmhmm_tag_count=1;

      my $tm_start_list='';
      my $tm_end_list='';

      foreach my $row(@{$gene_annotation{tmhmm_annotation}}){
         shift(@{$row});
         shift(@{$row});
         my $start=shift(@{$row});
         my $end=shift(@{$row});

            print qq*
              <tr align="center">
                <td>transmembrane                
                </td>
                <td>$start               
                </td> 
                <td>$end
                </td>
              </tr>                    
            *;

        $tm_start_list=$tm_start_list.";".$start;
        $tm_end_list=$tm_end_list.";".$end;
      }
      
      print qq*
      <tr>
        <td>
           <input name="TMHMM_START" id="TMHMM_START" type="hidden" value="$tm_start_list">
           <input name="TMHMM_END" id="TMHMM_END" type="hidden" value="$tm_end_list">
        </td>
      <tr>
      *;
}else{
   print qq*
      <tr align="center">
           <td>No transmembrane                
           </td>
           <td>0              
           </td> 
           <td>0
           </td>
        <td>
           <input name="TMHMM_START" id="TMHMM_START" type="hidden" value="0">
           <input name="TMHMM_END" id="TMHMM_END" type="hidden" value="0">
        </td>
       </tr>         
      *;
}
print qq*
</tbody>
</table>
<br>
<br>
</td>
</tr>

<tr valign="center">
<td style="font-weight: bold;">Phobius:<br>
</td>
<td style="vertical-align: top; width: 70%; text-align: left;"><br><br>
<table style="vertical-align: center; width: 100%; height:15%; margin-left: auto; margin-right: auto;" border="0" cellpadding="5" cellspacing="5">
<tbody>
<tr><th>Domain Name</th><th>Start</th><th>End</th><th>Description</th></tr>
*;

if($gene_annotation{phobius_annotation} ne "No annotation"){

      my $phobius_tag_count=1;

      my $phobius_domain_list='';
      my $phobius_start_list='';
      my $phobius_end_list='';
      my $phobius_description_list='';

      foreach my $row(@{$gene_annotation{phobius_annotation}}){
         shift(@{$row});
         shift(@{$row});
         my $domain_name=shift(@{$row});
         my $description=shift(@{$row});
         my $start=shift(@{$row});
         my $end=shift(@{$row});

            print qq*
              <tr align="center">
                <td>$domain_name                
                </td>
                <td>$start               
                </td> 
                <td>$end
                </td>
                <td>$description
                </td>
              </tr>                    
            *;

        $phobius_domain_list=$phobius_domain_list.";".$domain_name;
        $phobius_description_list=$phobius_description_list.";".$description;
        $phobius_start_list=$phobius_start_list.";".$start;
        $phobius_end_list=$phobius_end_list.";".$end;
      }
      
      print qq*
      <tr>
        <td>
           <input name="PHOBIUS_DOMAIN" id="PHOBIUS_DOMAIN" type="hidden" value="$phobius_domain_list">
           <input name="PHOBIUS_START" id="PHOBIUS_START" type="hidden" value="$phobius_start_list">
           <input name="PHOBIUS_END" id="PHOBIUS_END" type="hidden" value="$phobius_end_list">
           <input name="PHOBIUS_DESCRIPTION" id="PHOBIUS_DESCRIPTION" type="hidden" value="$phobius_description_list">
        </td>
      <tr>
      *;
}else{
   print qq*
      <tr align="center">
                <td>No domain                
                </td>
                <td>0              
                </td> 
                <td>0
                </td>
                <td>No annotation
                </td>
        <td>
           <input name="PHOBIUS_DOMAIN" id="PHOBIUS_DOMAIN" type="hidden" value="No domain">
           <input name="PHOBIUS_START" id="PHOBIUS_START" type="hidden" value="0">
           <input name="PHOBIUS_END" id="PHOBIUS_END" type="hidden" value="0">
           <input name="PHOBIUS_DESCRIPTION" id="PHOBIUS_DESCRIPTION" type="hidden" value="No annotation">
        </td>
       </tr>         
      *;
}
print qq*
</tbody>
</table>
<br>
<br>
</td>
</tr>

</tbody>
</table>
<br>
</form>
</body>
</html>
*;

