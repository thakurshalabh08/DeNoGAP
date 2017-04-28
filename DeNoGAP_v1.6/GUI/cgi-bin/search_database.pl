#!/usr/bin/perl -w
###### ABOUT: This Script Search Database ############
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
my $analysis=param('analysis');
my $core_threshold=param('core_gene_define');
my $show_result_list=param('show_result_list');
my $genome_with_homolog=param('genome_with_homolog');
my $genome_without_homolog=param('genome_without_homolog');
my $keyword=param('keyword');

#### list of genome name ####
$genome_with_homolog=~s/^\://g;
$genome_without_homolog=~s/^\://g;

my @genome_with_homolog=split(":",$genome_with_homolog);
my @genome_without_homolog=split(":",$genome_without_homolog);

my $result_id={};

##### Get Result for the Gene #######
if($analysis eq "core_gene" or $analysis eq "variable_gene"){

   $result_id=get_gene($db_dir,$odb_name,$analysis);

}
elsif($analysis eq "unique_gene"){
   ######## get column name #######
   my $sql_get_column_name="PRAGMA table_info(Profile)";

   my($column_data)=get_record($db_dir,$odb_name,$sql_get_column_name);

   shift(@{$column_data});
   
   my @result_id=();
   my %unique_gene_count=();

   ### For every genome selected ####
   foreach my $genome_name(@genome_with_homolog){

       my $unique_gene_sql=$genome_name."=1";

       ##### gene_count=1 and for other genomes gene_count=0
       foreach my $row(@{$column_data}){               
              my @column=@{$row};
              my $column_name=$column[1];   

              if($column_name ne $genome_name){
                $unique_gene_sql=$unique_gene_sql." AND ".$column_name."=0";
              } 
       }

       ##### Get unique gene list for the genome and add to main result_id array #####
       ##### also count number of unique gene in each genome and add to unique_gene_count hash ###
       #print $unique_gene_sql;
       my $sql_stmt_unique="SELECT id ,$genome_name FROM Profile WHERE ($unique_gene_sql)";
       my($unique_row_data)=get_record($db_dir,$odb_name,$sql_stmt_unique);      

        if(scalar(@{$unique_row_data})>0){
             my $gene_count=1;
             foreach my $row(@{$unique_row_data}){
               my @column=@{$row};
               my $count_genome=0;
               my $id=shift(@column);                            
               push(@result_id,$id);               
               $gene_count++;
             }           
           $unique_gene_count{$genome_name}=$gene_count;
       }
  }
$result_id=\@result_id;  
}

my $number_group=scalar(@{$result_id});

my($result_gene)=get_gene_info($db_dir,$db_name,$odb_name,$result_id,\@genome_with_homolog,$keyword);


######### FIND GENES ################
sub get_gene{

   my($db_dir)=(shift);
   my($odb_name)=(shift);
   my($analysis)=(shift);

   my $homolog_in_genome='';
   my $homolog_notin_genome='';
   my $show_column='';
   my $sql_stmt='';
   my $define_core_genome=0;
   my @result_id=();

   if(scalar(@genome_with_homolog)>=1) {
     $homolog_in_genome=join("=1 OR ",@genome_with_homolog);
     $homolog_in_genome=$homolog_in_genome."=1";
     $show_column=join(", ",@genome_with_homolog);
     $define_core_genome=int((($core_threshold/100)*scalar(@genome_with_homolog)));
   }

   if(scalar(@genome_without_homolog)>0){
      $homolog_notin_genome=join("=0 AND ",@genome_without_homolog);
      $homolog_notin_genome=$homolog_notin_genome."=0";
   } 

   if($homolog_notin_genome eq ''){     
     $sql_stmt="SELECT id ,$show_column FROM Profile WHERE ($homolog_in_genome)";
   }else{
     $sql_stmt="SELECT id ,$show_column FROM Profile WHERE ($homolog_in_genome) AND ($homolog_notin_genome)";
   }

   my($row_data)=get_record($db_dir,$odb_name,$sql_stmt); 
 
   if(scalar(@{$row_data})>0){
      foreach my $row(@{$row_data}){
          my @column=@{$row};
          my $count_genome=0;
          my $id=shift(@column);
             $id=~s/\://g;
        
          foreach(@column){
            $count_genome=$count_genome + $_;
          }

        if($analysis eq "core_gene"){ 
          if($count_genome>=$define_core_genome){
              push(@result_id,$id);
          }
        }elsif($analysis eq "variable_gene"){          
          if($count_genome>1 and $count_genome<$define_core_genome){              
              push(@result_id,$id);
          }elsif($define_core_genome==2){
              push(@result_id,$id);
          }
        }
      }
   }
   return(\@result_id);
}

#### get gene ids in each family #####
sub get_gene_info {

   my($db_dir)=(shift);
   my($db_name)=(shift);
   my($odb_name)=(shift);
   my(@result_id)=@{(shift)};
   my(@genome_with_homolog)=@{(shift)};
   my($keyword)=(shift);
   my $sql_stmt='';
   my $genome='';

   if(scalar(@genome_with_homolog)>1){
     $genome=join("' , '",@genome_with_homolog);
   }else{
     $genome=$genome_with_homolog[0];
   }
   
   my $group_ids=join("' , '",@result_id);
        
   $sql_stmt="SELECT * from MapGeneIdtoGeneFamily WHERE genefamily_id IN ('$group_ids') and species_abbreviation IN ('$genome')";

   my($result_gene)=get_record($db_dir,$odb_name,$sql_stmt);

   #### If gene list by keyword #####
   if($keyword ne ''){
      my @row_data=@{$result_gene};
      my @row_data2=@{$result_gene};
      my $list_gene_id='';      
           
      foreach my $row(sort{$a->[1] cmp $b->[1]}@row_data){
         my @row=@{$row}; 
         my $index_id=$row[0];
         my $group_id=$row[1];
         my $gene_id=$row[2];
         my $genome_abbrv=$row[3];     
         $list_gene_id=$list_gene_id.",'".$gene_id."'";         
      }
      $list_gene_id=~s/^\,//g;
 
      my $keyword_sql_stmt="Select Distinct(feature_id),description from GeneFeature where (feature_id IN ($list_gene_id) or protein_id IN ($list_gene_id)) and description LIKE '%$keyword%'";
  
      my($gene_with_keyword)=get_record($db_dir,$db_name,$keyword_sql_stmt); 
      my %gene_id_keyword=();
      my %group_id_keyword=();
      my @keyword_result_gene=();     

      if(scalar(@{$gene_with_keyword})>=1){           
           foreach my $row(@{$gene_with_keyword}){           
             my $feature_id=shift(@{$row}); 
             my $description=shift(@{$row});                                      
             $gene_id_keyword{$feature_id}=$description;             
           }
      } 

      my $prev_group_id='';
      my @group_With_keyword=();
      my $has_keyword=0;
      
      #### Search for group with atleast one gene having assigned keyword #### 
      foreach my $row(sort{$a->[1] cmp $b->[1]}@row_data){
         my @row=@{$row}; 
         my $index_id=$row[0];
         my $group_id=$row[1];
         my $gene_id=$row[2];
         my $genome_abbrv=$row[3];   
       
         if($group_id ne $prev_group_id){

            if($has_keyword eq 1){
              @keyword_result_gene=(@keyword_result_gene , @group_With_keyword);
            }
            $prev_group_id=$group_id;
            $has_keyword=0;
            @group_With_keyword=();

             if($gene_id_keyword{$gene_id}){                          
                $has_keyword=1;
             }
             push(@group_With_keyword,$row);

         }else{
            if($gene_id_keyword{$gene_id}){                          
               $has_keyword=1;
            }
            push(@group_With_keyword,$row);
         }
      } 
   
      $result_gene=\@keyword_result_gene;
   } 
  
   return($result_gene);
}

sub get_genome_name {
 
   my($db_dir)=(shift);
   my($db_name)=(shift);
   my($genome_abbrv)=(shift);
   my $genome_name='';

   my $sql_stmt="SELECT Distinct(genome_name),abbreviation from OrganismInfo WHERE abbreviation IN ('$genome_abbrv')";
   my ($row_data)=get_record($db_dir,$db_name,$sql_stmt);

   return($row_data);
}

sub get_cluster_id {
 
   my($db_dir)=(shift);
   my($hdb_name)=(shift);
   my($gene_id)=(shift);
   my($genome_abbrv)=(shift);

   my $hmm_group_stmt="SELECT * from GenetoSuperFamily WHERE gene_id='$gene_id' and genome_name='$genome_abbrv'";
   my ($hmm_group_data)=get_record($db_dir,$hdb_name,$hmm_group_stmt);
  
   return($hmm_group_data);
}

sub get_gene_description{

   my($db_dir)=(shift);
   my($db_name)=(shift);
   my($gene_id)=(shift);
   my($genome_abbrv)=(shift);
   my($keyword)=(shift);
   my $genome_description='';
   my $sql_stmt='';

   if($keyword eq ''){
     $sql_stmt="SELECT description from GeneFeature WHERE feature_id='$gene_id' or protein_id='$gene_id'";
   }else{     
     $sql_stmt="SELECT description from GeneFeature WHERE (feature_id='$gene_id' or protein_id='$gene_id') and description LIKE '%$keyword%'";
   }

   my ($row_data)=get_record($db_dir,$db_name,$sql_stmt);

    if(scalar(@{$row_data})>0){
      foreach my $row(@{$row_data}){
          $genome_description=shift(@{$row});
      }
    }else{
        $genome_description="Unannotated protein sequence";
    }
  return($genome_description);
}

sub get_record {

   my $db_dir=(shift);
   my $db_name=(shift);
   my $sql_stmt=(shift);

   my @row_data=();

   #print "$sql_stmt\n";

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


##################################### HTML CODE STARTS HERE FOR GENE FAMILY WISE SORTING TABLE ###############################
if($show_result_list eq "Gene Family ID" and $analysis ne "variable_gene" and $analysis ne "unique_gene"){

my %group_data=();

print qq*
<html>
<head>
<meta content="text/html; charset=ISO-8859-1"
http-equiv="content-type">
<title>result_analysis</title>

<link rel="stylesheet" type="text/css" href="jquery/jquery_datatable/jquery.dataTables.css">
<link rel="stylesheet" type="text/css" href="jquery/jquery_datatable/media/css/jquery.dataTables.css">
<link rel="stylesheet" type="text/css" href="jquery/jquery_datatable/extensions/TableTools/css/dataTables.tableTools.css">

<script type="text/javascript" src="jquery/jquery_tablesorter/jquery-latest.js"></script>
<script type="text/javascript" src="jquery/jquery_tablesorter/jquery.tablesorter.js"></script>
<script type="text/javascript" src="jquery/jquery_datatable/jquery.js"></script>
<script type="text/javascript" src="jquery/jquery_datatable/jquery-1.11.1.min.js"></script>
<script type="text/javascript" src="jquery/jquery_datatable/jquery.dataTables.min.js"></script>

<script type="text/javascript" src="jquery/jquery_datatable/media/js/jquery.js"></script>
<script type="text/javascript" src="jquery/jquery_datatable/media/js/jquery.dataTables.js"></script>
<script type="text/javascript" src="jquery/jquery_datatable/dataTables.tableTools.js"></script>


<script type="text/javascript">

    function change_genome(genome_name_id){

       var genome_name=document.getElementById(genome_name_id).options[document.getElementById(genome_name_id).selectedIndex].value;
        
       var gene_info_list=document.getElementsByClassName(genome_name);

       var table = \$('#result_table').DataTable(); 
      
       table.clear().draw();     

       for(var i=0;i<gene_info_list.length;i=i+4){

           var group_id_obj=gene_info_list[i].id;
           var homolog_group_id_obj=gene_info_list[i+1].id;
           var hmm_group_id_obj=gene_info_list[i+2].id;
           var description_obj=gene_info_list[i+3].id;
           var gene_id_obj=gene_info_list[i+4].id;
           var genome_name_obj=gene_info_list[i+5].id;

           var group_id_value=document.getElementById(group_id_obj).value;
           var homolog_group_id_value=document.getElementById(homolog_group_id_obj).value;
           var hmm_group_id_value=document.getElementById(hmm_group_id_obj).value;
           var description_value=document.getElementById(description_obj).value;
           var gene_id_value=document.getElementById(gene_id_obj).value;
           var genome_name_value=document.getElementById(genome_name_obj).value;

           var radioHtml = '<input name="gene_detail_radio" id="detail" type="radio" value='+gene_id_value+'>'+
                           '<input name="project_name" id="project_name" type="hidden" value="$project_name">'+
                           '<input name="db_name" id="db_name" type="hidden" value="$db_name">'+
                           '<input name="hdb_name" id="hdb_name" type="hidden" value="$hdb_name">'+
                           '<input name="odb_name" id="odb_name" type="hidden" value="$odb_name">'+
                           '<input name="db_dir" id="db_dir" type="hidden" value="$db_dir">'+
                           '<input name="report_dir" id="report_dir" type="hidden" value="$report_dir">';                
   
           table.row.add( [group_id_value,
                           homolog_group_id_value,
                           hmm_group_id_value,
                           gene_id_value,
                           genome_name_value,
                           description_value,
                           radioHtml
                          ] ).draw();
       }
    }

</script>

<script>
\$.fn.dataTable.TableTools.defaults.aButtons = [ "copy", "csv", "xls" ];

\$(document).ready(function() {
    \$('#result_table').DataTable( {"lengthMenu": [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, "All"]],
        dom: 'T<"clear">lfrtip',
        tableTools: {
            "sSwfPath": "jquery/jquery_datatable/extensions/TableTools/swf/copy_csv_xls.swf"
        }
    } );
} );
</script>

</head>

<body>
<div style="text-align: center;"><big><big><big>Search Result</big></big></big><br>
<br>
<form method="post" name="result_set" action="get_gene_detail.pl" target="_blank">
<div style="text-align: center;"></div>
<div style="text-align: left;">
</div>

<div>
<table id="genome" class="tablesorter" style="text-align: left; width: 100%; margin-left: auto; margin-right: auto;" border="0" cellpadding="1" cellspacing="1">
<tbody>
<tr>
<td colspan="4" rowspan="1" style="vertical-align: top;">
<span style="font-weight: bold;"> Genome Name:
<select name="genome_name" id="genome_name" onchange="change_genome(this.id)">
*;
##### selection list of genome names for faimlywise representation #####
my $list_genome=join("' , '",@genome_with_homolog);   #### Comman-sperated vales of genome names in which homolog is present #####
my $row_genome_name=get_genome_name($db_dir,$db_name,$list_genome); 

my $initial_selected_genome='';

if(scalar(@{$row_genome_name})>0){
     my $count_genome=1;
     foreach my $row(sort{$a->[1] cmp $b->[1]} @{$row_genome_name}){
            my $genome_name=shift(@{$row});
            my $genome_abbreviation=shift(@{$row});

            if($count_genome==1){
               print qq*
                  <option selected="selected" value="$genome_abbreviation">$genome_name</option>				 
               *;
               $initial_selected_genome=$genome_abbreviation;
               $count_genome++;
            }else{
               print qq*
                  <option value="$genome_abbreviation">$genome_name</option>
               *; 
            }
     }
} 
print qq*
</select>
</span> 
</td>
<td rowspan="1" style="vertical-align: top; text-align: right;"><input name="Submit_Gene_ID" value="Show Gene Information" type="submit"><br></td>
</td>
</tr>
</tbody>
</table>
<br>
<br>
</div>

<table id="result_table" class="tablesorter" style="text-align: left; width: 100%; margin-left: auto; margin-right: auto;" border="1" cellpadding="1" cellspacing="1">
<thead>
<tr>
<th
style="vertical-align: top; width: 162px; font-weight: bold; text-align: center;">Ortholog Group
ID<br>
</th>
<th
style="vertical-align: top; width: 162px; font-weight: bold; text-align: center;">Homolog Group
ID<br>
</th>
<th
style="vertical-align: top; width: 162px; font-weight: bold; text-align: center;">HMM Group
ID<br>
</th>
<th
style="vertical-align: top; font-weight: bold; width: 335px; text-align: center;">Gene
ID<br>
</th>
<th
style="vertical-align: top; font-weight: bold; width: 335px; text-align: center;">Genome Name<br>
</th>
<th
style="vertical-align: middle; font-weight: bold; width: 632px; text-align: center; white-space: nowrap;">&nbsp;Gene
Description<br>
</th>
<th
style="vertical-align: top; width: 48px; font-weight: bold; text-align: center;">Details<br>
</th>
</tr>
</thead>
<tbody>
*;

###### load result set in the html table #######
my @row_data=@{$result_gene};

if(scalar(@row_data)>0){
 
  my $prev_group_id='';
  my $gene_count=1;
  my $genome_name='';
  my $gene_description='';
  my $selected_genome_abbrv='';
  my $detail_gene='';  
  my @gene_feature=();
  my %genome_abbrv=();

  open(EXPORT_RESULT,">$Bin/../tmp/exported_result.txt");

  system("chmod 777 $Bin/../tmp/exported_result.txt");

  my @export_result=();

  #### sort by group id and genome name #####
  foreach my $row(sort{$a->[1] cmp $b->[1] || $a->[3] cmp $b->[3]}@row_data){

         shift(@{$row});
         my $group_id=shift(@{$row});
         my $gene_id=shift(@{$row});
         my $genome_abbrv=shift(@{$row});
         my $genome_name='';
         my $hmm_group='';
         my $homolog_group='';

         my $row_genome=get_genome_name($db_dir,$db_name,$genome_abbrv); 
         my $description=get_gene_description($db_dir,$db_name,$gene_id,$genome_abbrv,$keyword);
         my $hmm_group_data=get_cluster_id($db_dir,$hdb_name,$gene_id,$genome_abbrv);

         if(scalar(@{$row_genome})>0){
            foreach my $row(@{$row_genome}){
               $genome_name=shift(@{$row});
            }
         }

          if(scalar(@{$hmm_group_data})>0){
             foreach my $row(@{$hmm_group_data}){
                 shift(@{$row});
                 shift(@{$row});
                 shift(@{$row});
                 $hmm_group=shift(@{$row});
                 $homolog_group=shift(@{$row});
             }
          }

         if($hmm_group eq ''){
           $hmm_group = "NA"
         }
         if($homolog_group eq ''){
           $homolog_group = "NA"
         }

         $genome_abbrv{$genome_abbrv}=$genome_abbrv;
 
         $group_data{$genome_abbrv}->{$group_id}->{$genome_name}->{$gene_id}=$hmm_group.":".$homolog_group.":".$description;

         if($genome_abbrv eq $initial_selected_genome){

            print qq*
                  <tr>
                  <td style="vertical-align: top; width: 162px; text-align: left;" id="group_$group_id" class="group_id">$group_id<br>
                  </td>
                  <td style="vertical-align: top; width: 162px; text-align: left;" id="homolog_group_$group_id" class="group_id">$homolog_group<br>
                  </td>
                  <td style="vertical-align: top; width: 162px; text-align: left;" id="hmm_group_$group_id" class="group_id">$hmm_group<br>
                  </td>          
                  <td style="vertical-align: top; width: 162px; text-align: left;" id="listgene_$group_id" class="group_id">$gene_id<br>
                  </td>
                 <td style="vertical-align: top; width: 162px; text-align: left;" id="genome_$group_id" class="group_id">$genome_name<br>
                  </td>
                  <td style="vertical-align: top; width: 162px; text-align: left;" id="description_$group_id" class="group_id">$description<br>
                  </td>
                  <td style="vertical-align: top; width: 48px; text-align: center;" id="detail_$group_id" class="group_id">
                     <input checked="checked" name="gene_detail_radio" id="detail_$group_id" type="radio" value="$gene_id">
                     <input name="project_name" id="project_name" type="hidden" value="$project_name">
                     <input name="db_name" id="db_name" type="hidden" value="$db_name">
                     <input name="hdb_name" id="hdb_name" type="hidden" value="$hdb_name">
                     <input name="odb_name" id="odb_name" type="hidden" value="$odb_name">
                     <input name="db_dir" id="db_dir" type="hidden" value="$db_dir">
                     <input name="report_dir" id="report_dir" type="hidden" value="$report_dir">
                  </td>
                  </tr>
            *;
         }
  }
}          
print qq*
</tbody>
</table>
<br>
<table id="temp_data_table" style="text-align: left; width: 100%; margin-left: auto; margin-right: auto;" border="0" cellpadding="0" cellspacing="0">
<tr>
<td>
*;
foreach my $abbrv(keys %group_data){

    my %list_group_id=%{$group_data{$abbrv}};
  
    foreach my $group_id(keys %list_group_id){
 
       my %list_genome_name=%{$list_group_id{$group_id}};

       foreach my $genome_name(keys %list_genome_name){
     
          my %list_gene_id=%{$list_genome_name{$genome_name}};

           foreach my $gene_id(keys %list_gene_id){  

               my $value=$list_gene_id{$gene_id};
               my @column_value=split(":",$value);
               my $hmm_group=$column_value[0];
               my $homolog_group=$column_value[1];
               my $description=$column_value[2];
        
                 print qq*
                    <input class="$abbrv" name="hidden_genome_name" id="group_$gene_id" value="$group_id" type="hidden">
                    <input class="$abbrv" name="hidden_genome_name" id="homolog_group_$gene_id" value="$homolog_group" type="hidden">
                    <input class="$abbrv" name="hidden_genome_name" id="hmm_group_$gene_id" value="$hmm_group" type="hidden">
                    <input class="$abbrv" name="hidden_genome_name" id="desc_$gene_id" value="$description" type="hidden">
                    <input class="$abbrv" name="hidden_genome_name" id="$gene_id" value="$gene_id" type="hidden">
                    <input class="$abbrv" name="hidden_genome_name" id="genome_$gene_id" value="$genome_name" type="hidden">
                 *; 
           }
       }
    }
}
print qq*
</td>
</tr>
</table>
</form>
</div>
</body>
</html>
*;
}

##########################################################################
########## Show result sorted by gene and strains (for variable genes) #####

if($show_result_list eq "Gene ID" or $analysis eq "variable_gene" or $analysis eq "unique_gene"){
my $number_row=scalar(@{$result_gene});
print qq*
<html>
<head>
<meta content="text/html; charset=ISO-8859-1"
http-equiv="content-type">
<title>result_analysis</title>
<link rel="stylesheet" type="text/css" href="jquery/jquery_datatable/jquery.dataTables.css">
<link rel="stylesheet" type="text/css" href="jquery/jquery_datatable/media/css/jquery.dataTables.css">
<link rel="stylesheet" type="text/css" href="jquery/jquery_datatable/extensions/TableTools/css/dataTables.tableTools.css">

<script type="text/javascript" src="jquery/jquery_tablesorter/jquery-latest.js"></script>
<script type="text/javascript" src="jquery/jquery_tablesorter/jquery.tablesorter.js"></script>
<script type="text/javascript" src="jquery/jquery_datatable/jquery.js"></script>
<script type="text/javascript" src="jquery/jquery_datatable/jquery-1.11.1.min.js"></script>
<script type="text/javascript" src="jquery/jquery_datatable/jquery.dataTables.min.js"></script>

<script type="text/javascript" src="jquery/jquery_datatable/media/js/jquery.js"></script>
<script type="text/javascript" src="jquery/jquery_datatable/media/js/jquery.dataTables.js"></script>
<script type="text/javascript" src="jquery/jquery_datatable/dataTables.tableTools.js"></script>

<script>
\$.fn.dataTable.TableTools.defaults.aButtons = [ "copy", "csv", "xls" ];

\$(document).ready(function() {
    \$('#result_table').DataTable( {"lengthMenu": [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, "All"]],
        dom: 'T<"clear">lfrtip',
        tableTools: {
            "sSwfPath": "jquery/jquery_datatable/extensions/TableTools/swf/copy_csv_xls.swf"
        }
    } );
} );
</script>

</head>
<body>
<div style="text-align: center;"><big><big><big>Search Result</big></big></big><br>
<br>
<form method="post" name="result_set" action="get_gene_detail.pl" target="_blank">
<div style="text-align: right;">
 <input name="Submit_Gene_ID" value="Show Gene Information" type="submit"><br><br>
</div>
<table id="result_table" class="tablesorter" style="text-align: left; width: 100%; margin-left: auto; margin-right: auto;"
border="1" cellpadding="1" cellspacing="1">
<thead>
<tr>
<th
style="vertical-align: top; width: 162px; font-weight: bold; text-align: center;">Group
ID<br>
</th>
<th
style="vertical-align: top; width: 171px; font-weight: bold; text-align: center;">Gene
ID<br>
</th>
<th
style="vertical-align: top; font-weight: bold; width: 335px; text-align: center;">Genome
Name<br>
</th>
<th
style="vertical-align: middle; font-weight: bold; width: 632px; text-align: center; white-space: nowrap;">&nbsp;Gene
Description<br>
</th>
<th
style="vertical-align: top; width: 48px; font-weight: bold; text-align: center;">Details<br>
</th>
</tr>
</thead>
<tbody>
*;

###### load result set in the html table #######
my @row_data=@{$result_gene};

if(scalar(@row_data)>0){

  if($analysis eq "variable_gene" or $analysis eq "unique_gene"){
    @row_data=sort{$a->[3] cmp $b->[3] || $a->[1] cmp $b->[1]}@row_data
  }else{
    @row_data=sort{$a->[1] cmp $b->[1] || $a->[3] cmp $b->[3]}@row_data
  }
 
  foreach my $row(@row_data){
         shift(@{$row});
         my $group_id=shift(@{$row});
         my $gene_id=shift(@{$row});
         my $genome_abbrv=shift(@{$row});
         my $genome='';
         my $hmm_group='';
         my $homolog_group='';
       
        my $row_genome=get_genome_name($db_dir,$db_name,$genome_abbrv); 
        my $description=get_gene_description($db_dir,$db_name,$gene_id,$genome_abbrv,$keyword);
        my $hmm_group_data=get_cluster_id($db_dir,$hdb_name,$gene_id,$genome_abbrv);

        if(scalar(@{$row_genome})>0){
           foreach my $row(@{$row_genome}){
              $genome=shift(@{$row});              
           }
        }

         if(scalar(@{$hmm_group_data})>0){
             foreach my $row(@{$hmm_group_data}){
                 shift(@{$row});
                 shift(@{$row});
                 shift(@{$row});
                 $hmm_group=shift(@{$row});
                 $homolog_group=shift(@{$row});
             }
          }

               print qq*
                  <tr>
                    <td style="vertical-align: top; width: 162px; text-align: left;" id="group_$group_id" class="group_id">$group_id<br>
                    </td>
                    <td style="vertical-align: top; width: 162px; text-align: left;" id="homolog_group_$group_id" class="group_id">$homolog_group<br>
                    </td>
                    <td style="vertical-align: top; width: 162px; text-align: left;" id="hmm_group_$group_id" class="group_id">$hmm_group<br>
                    </td>          
                    <td style="vertical-align: top; width: 171px; text-align: left;" class="gene_id">$gene_id<br>
                    </td>
                    <td style="vertical-align: top; width: 335px;" id="genome_$group_id">$genome<br>
                        <input name="genome_name" id="genome_name" type="hidden" value="$genome_abbrv">
                    </td>
                    <td style="vertical-align: top; width: 632px;" id="desc_$group_id">$description<br>
                    </td>
                    <td style="vertical-align: top; width: 48px; text-align: center;">
                    <input checked="checked" name="gene_detail_radio" id="detail_$group_id" value="$gene_id" type="radio">
                     <input name="project_name" id="project_name" type="hidden" value="$project_name">
                     <input name="db_name" id="db_name" type="hidden" value="$db_name">
                     <input name="hdb_name" id="hdb_name" type="hidden" value="$hdb_name">
                     <input name="odb_name" id="odb_name" type="hidden" value="$odb_name">
                     <input name="db_dir" id="db_dir" type="hidden" value="$db_dir">
                     <input name="report_dir" id="report_dir" type="hidden" value="$report_dir">
                    </td>
                  </tr>
                 *;     
   }
}
print qq*
</tbody>
</table>
<br>
</form>
</div>
</body>
</html>
*;
}


