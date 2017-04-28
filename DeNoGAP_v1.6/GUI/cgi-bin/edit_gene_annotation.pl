#!/usr/bin/perl -w
###### ABOUT: This Script edit gene information for selected gene id and stores in database ############
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
my $aa_len=param('aa_len');
my $nuc_len=param('nuc_len');
my $product_description=param('product_description');
my $comment=param('comment');
my $homolog_group=param('homolog_group');
my $ortholog_group=param('ortholog_group');
my $hmm_group=param('hmm_group');
my $go_id_list=param('GO_ID');
my $go_category_list=param('GO_CATEGORY');
my $go_description_list=param('GO_DESC');
my $pfam_id_list=param('PFAM_ID');
my $pfam_name_list=param('PFAM_NAME');
my $pfam_desc_list=param('PFAM_DESC');
my $pfam_start_list=param('PFAM_START');
my $pfam_end_list=param('PFAM_END');
my $pfam_significance_list=param('PFAM_SIGNIFICANCE');
my $interpro_id_list=param('INTERPRO_ID');
my $interpro_name_list=param('INTERPRO_NAME');
my $pathway_id_list=param('PATHWAY_ID');
my $pathway_name_list=param('PATHWAY_NAME');
my $signalp_peptide=param('SIGNALP_PEPTIDE');
my $signalp_start=param('SIGNALP_START');
my $signalp_end=param('SIGNALP_END');
my $tmhmm_start_list=param('TMHMM_START');
my $tmhmm_end_list=param('TMHMM_END');
my $phobius_domain_list=param('PHOBIUS_DOMAIN');
my $phobius_description_list=param('PHOBIUS_DESCRIPTION');
my $phobius_start_list=param('PHOBIUS_START');
my $phobius_end_list=param('PHOBIUS_END');




my @go_id=split(";",$go_id_list);
my @go_category=split(";",$go_category_list);
my @go_desc=split(";",$go_description_list);

my @pfam_id=split(";",$pfam_id_list);
my @pfam_name=split(";",$pfam_name_list);
my @pfam_desc=split(";",$pfam_desc_list);
my @pfam_start=split(";",$pfam_start_list);
my @pfam_end=split(";",$pfam_end_list);
my @pfam_significance=split(";",$pfam_significance_list);

my @interpro_id=split(";",$interpro_id_list);
my @interpro_name=split(";",$interpro_name_list);

my @pathway_id=split(";",$pathway_id_list);
my @pathway_name=split(";",$pathway_name_list);

my @tmhmm_start=split(";",$tmhmm_start_list);
my @tmhmm_end=split(";",$tmhmm_end_list);

my @phobius_domain=split(";",$phobius_domain_list);
my @phobius_description=split(";",$phobius_description_list);
my @phobius_start=split(";",$phobius_start_list);
my @phobius_end=split(";",$phobius_end_list);

print qq*
<html>
<head>
<meta content="text/html; charset=ISO-8859-1"
http-equiv="content-type">
<title>Edit_Gene</title>
</head>
<body>
<form method="post" action="export_to_database.pl" name="Edit_Gene_Information">
<div style="text-align: center;"></div>
<table style="text-align: left; width: 1336px; height: 1486px;"
border="0" cellpadding="1" cellspacing="1">
<tbody>
<tr align="center">
<td colspan="2" rowspan="1"
style="vertical-align: top; width: 1145px;"><span
style="font-weight: bold;">EDIT GENE ANNOTATION</span><br>
</td>
</tr>
<tr>
<td style="width: 165px;"><span style="font-weight: bold;">Locus
Tag :</span><br>
</td>
<td style="width: 1145px;"><input maxlength="100" size="30"
name="locus_tag" value="$locus_id"></td>
</tr>
<tr>
<td style="width: 165px;"><span style="font-weight: bold;">Protein Id :</span><br>
</td>
<td style="width: 1145px;"><input maxlength="100" size="30" name="protein_id" value="$protein_id"></td>
</tr>
<tr>
<td style="width: 165px;"><span style="font-weight: bold;">External database Id :</span><br>
</td>
<td style="width: 1145px;"><input maxlength="100" size="100" name="dbxref" value="$dbxref_id"></td>
</tr>
<tr>
<td style="width: 165px; font-weight: bold;">Species / Strain
Name :<br>
</td>
<td style="width: 1145px; font-weight: bold;"><input
maxlength="100" size="100" name="species_name" value="$species_name"><br>
</td>
</tr>
<tr>
<td style="width: 165px; font-weight: bold;">Species
Abbreviation : <br>
</td>
<td style="width: 1145px; font-weight: bold;"><input
maxlength="30" size="30" name="abbreviation" value="$genome_abbreviation"><br>
</td>
</tr>
<tr>
<td style="width: 165px; font-weight: bold;">Genome ID :<br>
</td>
<td style="width: 1145px; font-weight: bold;"><input
maxlength="50" size="30" name="genome_id" value="$genome_id"><br>
</td>
</tr>
<tr>
<td style="width: 165px; font-weight: bold;">Genome Type :<br>
</td>
<td style="width: 1145px; font-weight: bold;"><input
maxlength="50" size="30" name="genome_type" value="$genome_type"><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Genome Length :<br>
</td>
<td style="font-weight: bold; width: 1145px;"><input
maxlength="30" size="30" name="genome_length" value="$genome_length"><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Index on Genome :<br>
</td>
<td style="font-weight: bold; width: 1145px;"><input
maxlength="30" size="30" name="index_genome" value="$index_genome"><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Feature Start : <br>
</td>
<td style="font-weight: bold; width: 1145px;"><input
maxlength="30" size="30" name="feature_start" value="$feature_start"><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Feature End : <br>
</td>
<td style="font-weight: bold; width: 1145px;"><input
maxlength="30" size="30" name="feature_end" value="$feature_end"><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Strand :<br>
</td>
<td style="font-weight: bold; width: 1145px;">
*;

if($strand eq "+"){
print qq*
<select size="2" name="strand">
<option selected="selected">+</option>
<option>-</option>
</select>
*;
}elsif($strand eq "-"){
print qq*
<select size="2" name="strand">
<option>+</option>
<option selected="selected">-</option>
</select>
*;
}else{
print qq*
<select size="2" name="strand">
<option selected="selected">+</option>
<option>-</option>
</select>
*;
}

print qq*
<br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Frame : <br>
</td>
<td style="font-weight: bold; width: 1145px;"><input
maxlength="30" size="30" name="frame" value="$frame"><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Gene Name : <br>
</td>
<td style="font-weight: bold; width: 1145px;"><input
maxlength="50" size="30" name="gene_name" value="$gene_name"><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Feature Type :<br>
</td>
<td style="font-weight: bold; width: 1145px;"><input
maxlength="50" size="30" name="feature_type" value="$feature_type"><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Protein Length : <br>
</td>
<td style="font-weight: bold; width: 1145px;"><input
maxlength="30" size="30" name="protein_length" value="$aa_len"><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Nucleotide Length :<br>
</td>
<td style="font-weight: bold; width: 1145px;"><input
maxlength="30" size="30" name="nucleotide_length" value="$nuc_len"><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Product
Description : <br>
</td>
<td style="font-weight: bold; width: 1194px;"><textarea cols="50" rows="5" name="product_description">$product_description</textarea><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Homolog Group :<br>
</td>
<td style="font-weight: bold; width: 1145px;"><label><input
readonly="readonly" maxlength="50" size="30" name="homolog_group" value="$homolog_group"></label><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">Ortholog Group :<br>
</td>
<td style="font-weight: bold; width: 1145px;"><input readonly="readonly"
maxlength="50" size="30" name="ortholog_group" value="$ortholog_group"><br>
</td>
</tr>
<tr>
<td style="font-weight: bold; width: 165px;">HMM Group :<br>
</td>
<td style="font-weight: bold; width: 1145px;"><input readonly="readonly"
maxlength="50" size="30" name="hmm_group" value="$hmm_group"><br>
</td>
</tr>
*;

my $list_go_id_tag_name='';
my $list_go_category_tag_name='';
my $list_go_description_tag_name='';

for(my $i=0;$i<scalar(@go_id);$i++){

   my $go_id=$go_id[$i];
   my $go_category=$go_category[$i];
   my $go_description=$go_desc[$i];

   if($go_id eq ''){
      next;
   }

print qq*
<tr>
<td
style="font-weight: bold; width: 129px; vertical-align: bottom; height: 46px;">GO
Id :<br>
</td>
<td
style="font-weight: bold; width: 1194px; vertical-align: bottom; height: 46px;"><input
maxlength="50" size="30" name="go_id_$i" value="$go_id">&nbsp;&nbsp;&nbsp;&nbsp; GO
Category : <input name="go_category_$i" size="30" maxlength="100" value="$go_category">&nbsp;&nbsp;&nbsp;&nbsp;
GO Description :&nbsp;&nbsp; <input name="go_description_$i" size="80" maxlength="200" value="$go_description"><br>
</td>
*;
  $list_go_id_tag_name=$list_go_id_tag_name.";"."go_id_$i";
  $list_go_category_tag_name=$list_go_category_tag_name.";"."go_category_$i";
  $list_go_description_tag_name=$list_go_description_tag_name.";"."go_description_$i";  
}

print qq*
<td>
<input name="list_go_id_name" id="list_go_id_name" type="hidden" value="$list_go_id_tag_name">
<input name="list_go_category_name" id="list_go_category_name" type="hidden" value="$list_go_category_tag_name">
<input name="list_go_desc_name" id="list_go_desc_name" type="hidden" value="$list_go_description_tag_name">
</td>
</tr>
*;

my $list_pfam_id_tag_name='';
my $list_pfam_name_tag_name='';
my $list_pfam_description_tag_name='';
my $list_pfam_start_tag_name='';
my $list_pfam_end_tag_name='';
my $list_pfam_significance_name='';

for(my $j=0;$j<scalar(@pfam_id);$j++){

     my $pfam_id=$pfam_id[$j];
     my $pfam_name=$pfam_name[$j];
     my $pfam_description=$pfam_desc[$j];
     my $pfam_start=$pfam_start[$j];
     my $pfam_end=$pfam_end[$j];
     my $pfam_significance=$pfam_significance[$j];

     if($pfam_id eq ''){
       next;
     }

print qq*

<tr>
<td
style="font-weight: bold; width: 129px; height: 46px; vertical-align: bottom;">Pfam
Id :<br>
</td>
<td
style="font-weight: bold; width: 1194px; height: 46px; vertical-align: bottom;"><input
maxlength="50" size="30" name="pfam_id_$j" value="$pfam_id">&nbsp;&nbsp;&nbsp;&nbsp; Pfam
Name :&nbsp;&nbsp; <input maxlength="100" size="30" name="pfam_name_$j" value="$pfam_name">&nbsp;&nbsp;&nbsp;&nbsp;
Pfam Description :&nbsp;<input maxlength="200" size="80" name="pfam_description_$j" value="$pfam_description"><br>
</td>
</tr>
<tr>
<td style="vertical-align: bottom; height: 46px;"><span
style="font-weight: bold;">Pfam Start :</span><br>
</td>
<td style="vertical-align: bottom; height: 46px;"><span
style="font-weight: bold;"><input maxlength="50" size="30"
name="pfam_start_$j" value="$pfam_start"></span>&nbsp;&nbsp;&nbsp;&nbsp; <span
style="font-weight: bold;">Pfam End :&nbsp;&nbsp; &nbsp;&nbsp; <input
maxlength="50" size="30" name="pfam_end_$j" value="$pfam_end"></span>&nbsp;&nbsp;&nbsp;&nbsp; <span
style="font-weight: bold;">Significance :&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;<input
maxlength="50" size="26" name="significance_$j" value="$pfam_significance"></span><br>
</td>

*;
 $list_pfam_id_tag_name=$list_pfam_id_tag_name.";"."pfam_id_$j";
 $list_pfam_name_tag_name=$list_pfam_name_tag_name.";"."pfam_name_$j";
 $list_pfam_description_tag_name=$list_pfam_description_tag_name.";"."pfam_description_$j";
 $list_pfam_start_tag_name=$list_pfam_start_tag_name.";"."pfam_start_$j";
 $list_pfam_end_tag_name=$list_pfam_end_tag_name.";"."pfam_end_$j";
 $list_pfam_significance_name=$list_pfam_significance_name.";"."significance_$j";
}

print qq*
<td>
<input name="list_pfam_id_tag" id="list_pfam_id_tag" type="hidden" value="$list_pfam_id_tag_name">
<input name="list_pfam_name_tag" id="list_pfam_name_tag" type="hidden" value="$list_pfam_name_tag_name">
<input name="list_pfam_desc_tag" id="list_pfam_desc_tag" type="hidden" value="$list_pfam_description_tag_name">
<input name="list_pfam_start_tag" id="list_pfam_start_tag" type="hidden" value="$list_pfam_start_tag_name">
<input name="list_pfam_end_tag" id="list_pfam_end_tag" type="hidden" value="$list_pfam_end_tag_name">
<input name="list_pfam_significance" id="list_pfam_significance" type="hidden" value="$list_pfam_significance_name">
</td>
</tr>
*;

my $list_interpro_id='';
my $list_interpro_name='';

for(my $k=0;$k<scalar(@interpro_id);$k++){

    my $interpro_id=$interpro_id[$k];
    my $interpro_description=$interpro_name[$k];

    if($interpro_id eq ''){
      next;
    }
    
print qq*

<tr>
<td
style="font-weight: bold; width: 129px; height: 46px; vertical-align: bottom;">InterPro
Id :<br>
</td>
<td
style="font-weight: bold; width: 1194px; height: 46px; vertical-align: bottom;"><input
maxlength="100" size="30" name="interpro_id_$k" value="$interpro_id">&nbsp;&nbsp;&nbsp;&nbsp;
InterPro Description : <input maxlength="100" size="80" name="interpro_description_$k" value="$interpro_description"><br>
</td>

*;
 $list_interpro_id=$list_interpro_id.";"."interpro_id_$k";
 $list_interpro_name=$list_interpro_name.";"."interpro_description_$k";
}

print qq*
<td>
<input name="list_interpro_id_tag" id="list_interpro_id_tag" type="hidden" value="$list_interpro_id">
<input name="list_interpro_name_tag" id="list_interpro_name_tag" type="hidden" value="$list_interpro_name">
</td>
</tr>
*;

my $list_pathway_id='';
my $list_pathway_description='';

for(my $p=0;$p<scalar(@pathway_id);$p++){

    my $pathway_id=$pathway_id[$p];
    my $pathway_description=$pathway_name[$p];

    if($pathway_id eq ''){
      next;
    }

print qq*
<tr>
<td
style="font-weight: bold; width: 129px; height: 46px; vertical-align: bottom;">Pathway
Id :<br>
</td>
<td
style="font-weight: bold; width: 1194px; height: 46px; vertical-align: bottom;"><input
maxlength="100" size="30" name="pathway_id_$p" value="$pathway_id">&nbsp;&nbsp;&nbsp;&nbsp;
Pathway Description : <input maxlength="100" size="80" name="pathway_description_$p" value="$pathway_description"><br>
</td>

*;
$list_pathway_id=$list_pathway_id.";"."pathway_id_$p";
$list_pathway_description=$list_pathway_description.";"."pathway_description_$p"; 

}

print qq*
<td>
<input name="list_pathway_id_tag" id="list_pathway_id_tag" type="hidden" value="$list_pathway_id">
<input name="list_pathway_name_tag" id="list_pathway_name_tag" type="hidden" value="$list_pathway_description">
</td>
</tr>
*;

print qq*
<tr>
<td
style="font-weight: bold; width: 129px; height: 46px; vertical-align: bottom;">Signal Peptide Start:<br>
</td>
<td
style="font-weight: bold; width: 1194px; height: 46px; vertical-align: bottom;"><input
maxlength="50" size="30" name="signalp_start" value="$signalp_start">&nbsp;&nbsp;&nbsp;&nbsp;
Signal Peptide End: &nbsp;&nbsp;&nbsp;&nbsp;<input maxlength="50" size="30" name="signalp_end" value="$signalp_end"><br>
</td>
<td>
<input name="signalp_start_tag" id="signalp_start_tag" type="hidden" value="$signalp_start">
<input name="signalp_end_tag" id="signalp_end_tag" type="hidden" value="$signalp_end">
</td>
</tr>
*;

my $list_tmhmm_start='';
my $list_tmhmm_end='';

for(my $p=0;$p<scalar(@tmhmm_start);$p++){

    my $tmhmm_start=$tmhmm_start[$p];
    my $tmhmm_end=$tmhmm_end[$p];

print qq*
<tr>
<td
style="font-weight: bold; width: 129px; height: 46px; vertical-align: bottom;">Transmembrane Start :<br>
</td>
<td
style="font-weight: bold; width: 1194px; height: 46px; vertical-align: bottom;"><input
maxlength="50" size="30" name="tmhmm_start_$p" value="$tmhmm_start">&nbsp;&nbsp;&nbsp;&nbsp;
Transmembrane End : <input maxlength="50" size="30" name="tmhmm_end_$p" value="$tmhmm_end"><br>
</td>

*;
$list_tmhmm_start=$list_tmhmm_start.";"."tmhmm_start_$p";
$list_tmhmm_end=$list_tmhmm_end.";"."tmhmm_end_$p"; 

}

print qq*
<td>
<input name="list_tmhmm_start_tag" id="list_tmhmm_start_tag" type="hidden" value="$list_tmhmm_start">
<input name="list_tmhmm_end_tag" id="list_tmhmm_end_tag" type="hidden" value="$list_tmhmm_end">
</td>
</tr>
<br>
*;

my $list_phobius_domain='';
my $list_phobius_description='';
my $list_phobius_start='';
my $list_phobius_end='';


for(my $p=0;$p<scalar(@phobius_domain);$p++){

    my $phobius_domain_name=$phobius_domain[$p];
    my $phobius_description=$phobius_description[$p];
    my $phobius_start=$phobius_start[$p];
    my $phobius_end=$phobius_end[$p];

print qq*
<tr>
<td
style="font-weight: bold; width: 129px; height: 46px; vertical-align: bottom;">Phobius Domain :<br>
</td>
<td
style="font-weight: bold; width: 1194px; height: 46px; vertical-align: bottom;"><input
maxlength="100" size="30" name="phobius_domain_$p" value="$phobius_domain_name">&nbsp;&nbsp;&nbsp;&nbsp;
Phobius Start : <input maxlength="10" size="10" name="phobius_start_$p" value="$phobius_start">&nbsp;&nbsp;&nbsp;&nbsp;
Phobius End : <input maxlength="10" size="10" name="phobius_end_$p" value="$phobius_end">&nbsp;
Phobius Description : <input maxlength="200" size="60" name="phobius_description_$p" value="$phobius_description"><br>
</td>

*;
$list_phobius_domain=$list_phobius_domain.";"."phobius_domain_$p";
$list_phobius_start=$list_phobius_start.";"."phobius_start_$p";
$list_phobius_end=$list_phobius_end.";"."phobius_end_$p";
$list_phobius_description=$list_phobius_description.";"."phobius_description_$p";

}

print qq*
<td>
<input name="list_phobius_domain_tag" id="list_phobius_domain_tag" type="hidden" value="$list_phobius_domain">
<input name="list_phobius_description_tag" id="list_phobius_description_tag" type="hidden" value="$list_phobius_description">
<input name="list_phobius_start_tag" id="list_phobius_start_tag" type="hidden" value="$list_phobius_start">
<input name="list_phobius_end_tag" id="list_phobius_end_tag" type="hidden" value="$list_phobius_end">
</td>
</tr>
*;

print qq*
<tr>
<td
style="vertical-align: middle; font-weight: bold; text-align: left;">Comments
:<br>
</td>
<td
style="vertical-align: middle; font-weight: bold; text-align: left;"><textarea
cols="98" rows="10" name="comments">$comment</textarea></td>
</tr>
<tr>
<td colspan="2" rowspan="1" style="vertical-align: top;"> 
   <input name="db_dir" id="db_dir" type="hidden" value="$db_dir">
   <input name="c_database" id="c_database" type="hidden" value="$cdb_name">
   <input name="h_database" id="h_database" type="hidden" value="$hdb_name">
   <input name="o_database" id="o_database" type="hidden" value="$odb_name">
   <input name="submit_database" value="Submit to database" type="submit"><br>
</td>
</tr>
</tbody>
</table>
<br>
</form>
</body>
</html>

*;




