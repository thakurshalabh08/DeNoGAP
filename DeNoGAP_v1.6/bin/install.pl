#!/usr/bin/perl -w
###### ABOUT: This Script install the require perl modules and software to run the genomics pipeline ####
###### AUTHOR:Shalabh Thakur###################################################################

#### INSTRUCTIONS TO INSTALL THE PACKAGE #####

### Required Input: <path to home directory>, directory in which all required softwares should be installed #####
### Permissions: Installation of Perl module and some external program would need root permission ########
### USAGE: perl install.pl <complete path of home directory> ####

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use Getopt::Long;

my $HOME=undef;

chomp($HOME);


########## URL FOR EXTERNAL PROGRAM ########
########## USER Can change this to get newer versions of the programs ######

my $muscle="http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz";
my $kalign="http://msa.sbc.su.se/downloads/kalign/current.tar.gz";
my $mcl="http://micans.org/mcl/src/mcl-12-135.tar.gz";
my $hmmer="http://selab.janelia.org/software/hmmer3/3.1b1/hmmer-3.1b1-linux-intel-x86_64.tar.gz";
my $phylip="http://evolution.gs.washington.edu/phylip/download/phylip-3.696.tar.gz";
my $glimmer="http://ccb.jhu.edu/software/glimmer/glimmer302b.tar.gz";
my $prodigal="http://prodigal.googlecode.com/files/prodigal.v2_60.linux";
my $fraggenescan="http://omics.informatics.indiana.edu/mg/get.php?software=FragGeneScan1.16.tar.gz";
my $genemark="http://opal.biology.gatech.edu";
my $sqlite="https://sqlite.org/2014/sqlite-autoconf-3080401.tar.gz";
my $emboss="sudo apt-get install emboss";
my $blast="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.31+-x64-linux.tar.gz";
my $iprscan="ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.9-50.0/interproscan-5.9-50.0-64-bit.tar.gz";
my $iprscan_md="ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.9-50.0/interproscan-5.9-50.0-64-bit.tar.gz.md5";

####################################################################################################################

######## Do not change anything beyond this line ###################################################################

my @perl_std_module=("FindBin",
                     "Env",
                     "Exporter",
                     "Getopt::Long",
                     "File::Basename",
                     "File::Copy",
                     "Parallel::ForkManager",
                     "List::MoreUtils",
                     "List::Util",
                     "File::Path",
                     "Hash::Merge",
                     "DBI",
                     "Statistics::Basic",
                     "Statistics::R",
                     "Tie::File",
                     "Sort::Fields");

my @bioperl_module= ("Bio::Perl", "Bio::SeqIO", "Bio::Seq", "Bio::SearchIO", "Bio::Tools::Phylo::Phylip::ProtDist","Bio::AlignIO","Bio::SeqFeature::Generic","Bio::Annotation::Collection","Bio::Annotation::Comment");

my @interproscan_module=("CGI",
                         "English",
                         "File::Spec::Functions",
                         "FileHandle",
                         "IO::Scalar",
                         "IO::String",
                         "Mail::Send",
                         "Sys::Hostname",
                         "URI::Escape",
                         "XML::Parser",
                         "XML::Quote");

my %external_prog=('muscle'=>"$muscle",
                   'kalign'=>"$kalign", 
                   'mcl'=>"$mcl", 
                   'hmmer'=>"$hmmer", 
                   'phylip'=>"$phylip", 
                   'glimmer'=>"$glimmer", 
                   'prodigal'=>"$prodigal", 
                   'fragscan'=>"$fraggenescan", 
                   'genemark'=>"$genemark",
                   'sqlite'=>"$sqlite", 
                   'emboss'=>"$emboss",
                   'blast'=>"$blast",
                   'iprscan'=>"$iprscan",
                   'iprscan_md'=>"$iprscan_md"              
                   );


print "Checking if you have standard perl modules installed\n";

foreach my $mod(@perl_std_module){
       eval("use $mod");
       #if error or not installed
       if($@){
         print "$mod is not found on your system or path\n\n";
         print "Do you want to Install $mod now? (Y | N):";
         my $install_now=<STDIN>;   
           if($install_now=~/Y/i){
              eval {system("sudo cpan install $mod")};             
              warn("Error:".$!) if $@;             
           }
       }else{
         print "$mod already installed on your system\n\n";
       } 
}

print "Checking if you have required Bio-Perl modules installed\n";

foreach my $mod(@bioperl_module){
       eval("use $mod");
       #if error or not installed
       if($@){
         print "$mod is not found on your system or path\n\n";
         print "Do you want to Install $mod now? (Y | N):";
         my $install_now=<STDIN>;   
           if($install_now=~/Y/i){
              eval {system("sudo cpan install $mod")};             
              warn("Error:".$!) if $@;             
           }
       }else{
         print "$mod already installed on your system\n\n";
       } 
}

print "Checking if you have required modules for Interproscan installed\n";

foreach my $mod(@interproscan_module){
       eval("use $mod");
       #if error or not installed
       if($@){
         print "$mod is not found on your system or path\n\n";
         print "Do you want to Install $mod now? (Y | N):";
         my $install_now=<STDIN>;   
           if($install_now=~/Y/i){
              eval {system("sudo cpan install $mod")};             
              warn("Error:".$!) if $@;             
           }
       }else{
         print "$mod already installed on your system\n\n";
       }
}


print "Checking if you have required external program installed\n\n";

foreach my $ext_prog(keys %external_prog){

    ### Install MUSCLE ###
    if($ext_prog eq 'muscle'){

       print "Do you want to Install MUSCLE alignment program on your system? ( Y or N):";
       my $option_install=<STDIN>;
       chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
       }

       $HOME=$Bin;

       print "Installation path for MUSCLE [$HOME/../exe]? (Y or N):";
       my $accept_path=<STDIN>;

       chomp($accept_path);

       SET_MUSCLE_PATH:

       if($accept_path eq "N"){
            print "Enter New installation path for MUSCLE:";
            $HOME=<STDIN>;
            chomp($HOME);
       }elsif($accept_path eq "Y"){          
           $HOME="$HOME/../exe";
       }else{
           
           print "Invalid Option, Try Again\n";

           print "Installation path for MUSCLE [$HOME/../exe]? (Y or N):";
           $accept_path=<STDIN>;
           chomp($accept_path);
   
           goto SET_MUSCLE_PATH;
       }

       eval{
        mkdir("$HOME/muscle");        
        system("wget $external_prog{$ext_prog} -P $HOME/muscle");
        system("chmod -R 777 $HOME/muscle");
        system("tar -xvzf $HOME/muscle/*.tar.gz -C $HOME/muscle");
        my @program_path=split(/\//,$external_prog{$ext_prog});
        my $program_name=pop(@program_path);
           $program_name=~s/(\.)(tar)(\.)(gz)//g;
           chomp($program_name);
        system("mv $HOME/muscle/$program_name $HOME/muscle/muscle");
        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/muscle)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/muscle\n";
           close PROFILE;
        }else{
           print "\n\nDeNoGAP did not recognize the shell you are using.  Sorry, but you will have to add MUSCLE program to your path manually.\n\n";

           print "Add following line to your .bashrc or .cshrc file\n\n";

           print "export PATH=\$PATH:$HOME/muscle/muscle\n\n";
 
           print "press [ENTER] to continue.\n\n";
           <>;
        }
      };
      if($@){
         die($!);
      }

    }
    #### INSTALL KALIGN #####
    elsif($ext_prog eq 'kalign'){

       print "Do you want to Install Kalign alignment program on your system? ( Y or N): ";
       my $option_install=<STDIN>;
       chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
       }

       $HOME=$Bin;

       print "Installation path for Kalign [$HOME/../exe]? (Y or N):";
       my $accept_path=<STDIN>;
       chomp($accept_path);

       SET_KALIGN_PATH:

       if($accept_path eq "N"){
            print "Enter New installation path for Kalign:";
            $HOME=<STDIN>;
            chomp($HOME);
       }elsif($accept_path eq "Y"){          
           $HOME="$HOME/../exe";

       }else{
           
           print "Invalid Option, Try Again\n";

           print "Installation path for Kalign [$HOME/../exe]? (Y or N):";
           $accept_path=<STDIN>;
           chomp($accept_path);
           goto SET_KALIGN_PATH;
       }

      eval{
        mkdir("$HOME/kalign");        
        system("wget $external_prog{$ext_prog} -P $HOME/kalign");
        system("chmod -R 777 $HOME/kalign");
        system("tar -xvzf $HOME/kalign/*.tar.gz -C $HOME/kalign");
        chdir("$HOME/kalign");
        system("./configure");
        system("make");
        system("sudo make install");
        chdir($Bin);
      };
      if($@){
         die($!);
      }
    }
    ##### INSTALL HMMER ######
    elsif($ext_prog eq 'hmmer'){

       print "Do you want to Install HMMER program on your system? ( Y or N): ";
       my $option_install=<STDIN>;
       chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
       }

       $HOME=$Bin;

       print "Installation path for Hmmer [$HOME/../exe]? (Y or N):";
       my $accept_path=<STDIN>;
       chomp($accept_path);

       SET_HMMER_PATH:

       if($accept_path eq "N"){
            print "Enter New installation path for Hmmer:";
            $HOME=<STDIN>;
            chomp($HOME);
       }elsif($accept_path eq "Y"){          
           $HOME="$HOME/../exe";

       }else{
           
           print "Invalid Option, Try Again\n";

           print "Installation path for Hmmer [$HOME/../exe]? (Y or N):";
           $accept_path=<STDIN>;
           chomp($accept_path);

           goto SET_HMMER_PATH;
       }

      eval{
        mkdir("$HOME/hmmer");        
        system("wget $external_prog{$ext_prog} -P $HOME/hmmer");
        system("chmod -R 777 $HOME/hmmer");
        system("tar -xvzf $HOME/hmmer/*.tar.gz -C $HOME/hmmer");

        chdir("$HOME/hmmer");    
        my $program_name=`ls -d */`;
        $program_name=~s/\///g;        
        chomp($program_name);
        print "$program_name\n";
        chdir("$program_name");
        system("./configure --prefix=$HOME/hmmer");
        system("make");
        system("make install");
        chdir($Bin);

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/hmmer/bin)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/hmmer/bin\n";
           close PROFILE;
        }else{
           print "\n\nDeNoGAP did not recognize the shell you are using.  Sorry, but you will have to add HMMER program to your path manually.\n\n";

           print "Add following line to your .bashrc or .cshrc file\n\n";

           print "export PATH=\$PATH:$HOME/hmmer/bin\n\n";

           print "press [ENTER] to continue.\n\n";
           <>; 
        }
      };
      if($@){
         die($!);
      }
    }
    #### INSTALL MCL #####
    elsif($ext_prog eq 'mcl'){

       print "Do you want to Install MCL program on your system? ( Y or N): ";
       my $option_install=<STDIN>;
       chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
       }

       $HOME=$Bin;

       print "Installation path for MCL [$HOME/../exe]? (Y or N):";
       my $accept_path=<STDIN>;
       chomp($accept_path);

       SET_MCL_PATH:

       if($accept_path eq "N"){
            print "Enter New installation path for MCL:";
            $HOME=<STDIN>;
            chomp($HOME);
       }elsif($accept_path eq "Y"){          
           $HOME="$HOME/../exe";

       }else{
           
           print "Invalid Option, Try Again\n";

           print "Installation path for MCL [$HOME/../exe]? (Y or N):";
           $accept_path=<STDIN>;
           chomp($accept_path);

           goto SET_MCL_PATH;
       }

       eval{
        mkdir("$HOME/mcl");        
        system("wget $external_prog{$ext_prog} -P $HOME/mcl");
        system("chmod -R 777 $HOME/mcl");
        system("tar -xvzf $HOME/mcl/*.tar.gz -C $HOME/mcl");
    
        chdir("$HOME/mcl");    
        my $program_name=`ls -d */`;
        $program_name=~s/\///g;        
        chomp($program_name);

        chdir("$program_name");
        system("./configure --prefix=$HOME/mcl");
        system("make");
        system("make install");
        chdir($Bin);

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/mcl/bin)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           print "$ENV{HOME}/.bashrc\n";
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/mcl/bin\n";
           close PROFILE;
        }else{
           print "\n\nDeNoGAP did not recognize the shell you are using.  Sorry, but you will have to add MCL program to your path manually.\n\n";

           print "Add following line to your .bashrc or .cshrc file\n\n";

           print "export PATH=\$PATH:$HOME/mcl/bin\n\n";

           print "press [ENTER] to continue.\n\n";
           <>; 
        }
      };
      if($@){
         die($!);
      }
    }
    #### INSTALL PHYLIP ####
    elsif($ext_prog eq 'phylip'){

       print "Do you want to Install Phylip package on your system? ( Y or N): ";
       my $option_install=<STDIN>;
       chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
       }

       $HOME=$Bin;

       print "Installation path for Phylip [$HOME/../exe]? (Y or N):";
       my $accept_path=<STDIN>;
       chomp($accept_path);

       SET_PHYLIP_PATH:

       if($accept_path eq "N"){
            print "Enter New installation path for Phylip:";
            $HOME=<STDIN>;
            chomp($HOME);
       }elsif($accept_path eq "Y"){          
           $HOME="$HOME/../exe";

       }else{
           
           print "Invalid Option, Try Again\n";

           print "Installation path for Phylip [$HOME/../exe]? (Y or N):";
           $accept_path=<STDIN>;
           chomp($accept_path);

           goto SET_PHYLIP_PATH;
       }

       eval{
        mkdir("$HOME/phylip");        
        system("wget $external_prog{$ext_prog} -P $HOME/phylip");
        system("chmod -R 777 $HOME/phylip");
        system("tar -xvzf $HOME/phylip/*.tar.gz -C $HOME/phylip");

        chdir("$HOME/phylip");
        my $program_name=`ls -d */`;
           $program_name=~s/\///g;        
        chomp($program_name);

        chdir("$program_name/src");
        system("mv Makefile.unx Makefile");
        system("make install");
        chdir($Bin);

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/phylip/$program_name/exe)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/phylip/$program_name/exe\n";
           close PROFILE;
        }else{
           print "\n\nDeNoGAP did not recognize the shell you are using.  Sorry, but you will have to add PHYLIP program to your path manually.\n\n";

           print "Add following line to your .bashrc or .cshrc file\n\n";

           print "export PATH=\$PATH:$HOME/phylip/$program_name/exe\n\n";

           print "press [ENTER] to continue.\n\n";
           <>; 
        }
      };
      if($@){
         die($!);
      }

    }
    #### INSTALL GLIMMER ####
    elsif($ext_prog eq 'glimmer'){

        print "Do you want to Install Glimmer on your system? ( Y or N): ";
        my $option_install=<STDIN>;
        chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
       }

       $HOME=$Bin;

       print "Installation path for Glimmer [$HOME/../exe]? (Y or N):";
       my $accept_path=<STDIN>;
       chomp($accept_path);

       SET_GLIMMER_PATH:

       if($accept_path eq "N"){
            print "Enter New installation path for Glimmer:";
            $HOME=<STDIN>;
            chomp($HOME);
       }elsif($accept_path eq "Y"){          
           $HOME="$HOME/../exe";

       }else{
           
           print "Invalid Option, Try Again\n";

           print "Installation path for Glimmer [$HOME/../exe]? (Y or N):";
           $accept_path=<STDIN>;
           chomp($accept_path);

           goto SET_GLIMMER_PATH;
       }

        eval{
        mkdir("$HOME/glimmer");        
        system("wget $external_prog{$ext_prog} -P $HOME/glimmer");
        system("chmod -R 777 $HOME/glimmer");
        system("tar -xvzf $HOME/glimmer/*.tar.gz -C $HOME/glimmer"); 
        chdir("$HOME/glimmer");        
        my $program_name=`ls -d */`;
           $program_name=~s/\///g;        
        chomp($program_name);
        chdir("$program_name/src"); 
        system("make");
        chdir($Bin);

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/glimmer/$program_name/bin)\n";
           print PROFILE "set path=(\$path $HOME/glimmer/$program_name/scripts)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/glimmer/$program_name/bin\n";
           print PROFILE "export PATH=\$PATH:$HOME/glimmer/$program_name/scripts\n";
           close PROFILE;
        }else{
           print "\n\nDeNoGAP did not recognize the shell you are using.  Sorry, but you will have to add GLIMMER program to your path manually.\n\n";

           print "Add following line to your .bashrc or .cshrc file\n\n";

           print "export PATH=\$PATH:$HOME/glimmer/$program_name/bin\n\n";

           print "export PATH=\$PATH:$HOME/glimmer/$program_name/scripts\n\n";

           print "press [ENTER] to continue.\n\n";
           <>; 
        }
      };
      if($@){
         die($!);
      }

    }
    ### INSTALL PRODIGAL ####
    elsif($ext_prog eq 'prodigal'){

       print "Do you want to Install Prodigal on your system? ( Y or N): ";
       my $option_install=<STDIN>;
       chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
       }

       $HOME=$Bin;

       print "Installation path for Prodigal [$HOME/../exe]? (Y or N):";
       my $accept_path=<STDIN>;
       chomp($accept_path);

       SET_PRODIGAL_PATH:

       if($accept_path eq "N"){
            print "Enter New installation path for Prodigal:";
            $HOME=<STDIN>;
            chomp($HOME);
       }elsif($accept_path eq "Y"){          
           $HOME="$HOME/../exe";

       }else{
           
           print "Invalid Option, Try Again\n";

           print "Installation path for Prodigal [$HOME/../exe]? (Y or N):";
           $accept_path=<STDIN>;
           chomp($accept_path);

           goto SET_PRODIGAL_PATH;
       }

       eval{
        mkdir("$HOME/prodigal");        
        system("wget $external_prog{$ext_prog} -P $HOME/prodigal");
        system("chmod -R 777 $HOME/prodigal");       
        my @program_path=split(/\//,$external_prog{$ext_prog});
        my $program_name=pop(@program_path);
        chomp($program_name);
        system("mv $HOME/prodigal/$program_name $HOME/prodigal/prodigal");

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/prodigal)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/prodigal\n";
           close PROFILE;
        }else{
           print "\n\nDeNoGAP did not recognize the shell you are using.  Sorry, but you will have to add Prodigal program to your path manually.\n\n";

           print "Add following line to your .bashrc or .cshrc file\n\n";

           print "export PATH=\$PATH:$HOME/prodigal/prodigal\n\n";

           print "press [ENTER] to continue.\n\n";
           <>; 
        }
      };
      if($@){
         die($!);
      }

   }
   ##### INSTALL FRAGSCAN #####
   elsif($ext_prog eq 'fragscan'){

       print "Do you want to Install FragGeneScan on your system? ( Y or N): ";
       my $option_install=<STDIN>;
       chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
       }

       $HOME=$Bin;

       print "Installation path for FragGeneScan [$HOME/../exe]? (Y or N):";
       my $accept_path=<STDIN>;
       chomp($accept_path);

       SET_FRAGSCAN_PATH:

       if($accept_path eq "N"){
            print "Enter New installation path for FragGeneScan:";
            $HOME=<STDIN>;
            chomp($HOME);
       }elsif($accept_path eq "Y"){          
           $HOME="$HOME/../exe";

       }else{
           
           print "Invalid Option, Try Again\n";

           print "Installation path for FragGeneScan [$HOME/../exe]? (Y or N):";
           $accept_path=<STDIN>;
           chomp($accept_path);

           goto SET_FRAGSCAN_PATH;
       }

       eval{
        mkdir("$HOME/frag_gene_scan");   
   
        DOWNLOAD_FRAGGENESCAN:

        print "Please Download FragGeneScan manually from $external_prog{$ext_prog}\n\n";
        print "Copy dowloaded FragGeneScan tar file into $HOME/frag_gene_scan directory.\n\n";
        print "Change name of .tar.gz file to fraggenescan.tar.gz and press enter to continue.\n\n";
        <>;  

        unless(-s "$HOME/frag_gene_scan/fraggenescan.tar.gz"){
           goto DOWNLOAD_FRAGGENESCAN;
        }

        system("chmod -R 777 $HOME/frag_gene_scan");
        system("tar -xvzf $HOME/frag_gene_scan/*.tar.gz -C $HOME/frag_gene_scan");

        chdir("$HOME/frag_gene_scan");
        my $program_name=`ls -d */`;
           $program_name=~s/\///g;        
           chomp($program_name);
        chdir($Bin);

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/frag_gene_scan/$program_name)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/frag_gene_scan/$program_name\n";
           close PROFILE;
        }else{
           print "\n\nDeNoGAP did not recognize the shell you are using.  Sorry, but you will have to add FragGeneScan program to your path manually.\n\n";

           print "Add following line to your .bashrc or .cshrc file\n\n";

           print "export PATH=\$PATH:$HOME/frag_gene_scan/$program_name\n\n";

           print "press [ENTER] to continue.\n\n";
           <>; 
        }
      };
      if($@){
         die($!);
      }
   }
   #### INSTALL GENEMARK ####
   elsif($ext_prog eq 'genemark'){

       print "Do you want to Install GeneMark on your system? ( Y or N): ";
       my $option_install=<STDIN>;
       chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
       }

       $HOME=$Bin;

       print "Installation path for GeneMark [$HOME/../exe]? (Y or N):";
       my $accept_path=<STDIN>;
       chomp($accept_path);

       SET_GENEMARK_PATH:

       if($accept_path eq "N"){
            print "Enter New installation path for GeneMark:";
            $HOME=<STDIN>;
            chomp($HOME);
       }elsif($accept_path eq "Y"){          
           $HOME="$HOME/../exe";

       }else{
           
           print "Invalid Option, Try Again\n";

           print "Installation path for GeneMark [$HOME/../exe]? (Y or N):";
           $accept_path=<STDIN>;
           chomp($accept_path);

           goto SET_GENEMARK_PATH;
       }

        eval{
        mkdir("$HOME/genemark"); 

        DOWNLOAD_GENEMARK:

        print "Please download genemark for prokaryotes manually from $external_prog{$ext_prog}\n\n";
        print "Copy dowloaded genemark tar file and genemark key file into $HOME/genemark directory.\n\n";
        print "Change name of .tar.gz file to genemark.tar.gz and press enter to continue.\n\n";
        <>;

        unless(-s "$HOME/genemark/genemark.tar.gz"){
           goto DOWNLOAD_GENEMARK;
        }else{
           
           print "Cannot find compressed tar file to extract required files forGeneMark installation\n\n";
        }

        system("chmod -R 777 $HOME/genemark");

        system("tar -xvzf $HOME/genemark/*.tar.gz -C $HOME/genemark");

        chdir("$HOME/genemark"); 
        my $program_name=`ls -d */`;
        chomp($program_name);   

        COPY_GM_KEY:

        chdir("$HOME/genemark/$program_name/gmsuite");

        system("cp -T gm_key ~/.gm_key");

        print "gm_key copied to ~/ directory\n";

        chdir("$HOME/genemark");

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/genemark/$program_name/gmsuite)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/genemark/$program_name/gmsuite\n";
           close PROFILE;
        }else{
           print "\n\nDeNoGAP did not recognize the shell you are using.  Sorry, but you will have to add GeneMark program to your path manually.\n\n";

           print "Add following line to your .bashrc or .cshrc file\n\n";

           print "export PATH=\$PATH:$HOME/genemark/$program_name/gmsuite\n\n";

           print "press [ENTER] to continue.\n\n";
           <>; 
        }
      };
      if($@){
         die($!);
      }
   }
   ##### INSTALL SQLITE ######
   elsif($ext_prog eq 'sqlite'){

       print "Do you want to Install SQLITE on your system? ( Y or N): ";
       my $option_install=<STDIN>;
          chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
        }

       $HOME=$Bin;

       print "Installation path for SQLite [$HOME/../exe]? (Y or N):";
       my $accept_path=<STDIN>;
       chomp($accept_path);

       SET_SQLITE_PATH:

       if($accept_path eq "N"){
            print "Enter New installation path for SQLite:";
            $HOME=<STDIN>;
            chomp($HOME);
       }elsif($accept_path eq "Y"){          
           $HOME="$HOME/../exe";

       }else{
           
           print "Invalid Option, Try Again\n";

           print "Installation path for SQLite [$HOME/../exe]? (Y or N):";
           $accept_path=<STDIN>;
           chomp($accept_path);

           goto SET_SQLITE_PATH;
       }

     if($option_install eq "Y" or $option_install eq "y"){

        my $path="$HOME/sqlite";
        
        print "The default path of installation for SQLITE is [$path], Give new path if you want to change: ";
        my $option_change_path=<STDIN>;
           chomp($option_change_path);

        if($option_change_path ne ''){            
           $path=$option_change_path; 
        }

        eval{

        mkdir($path);
        system("wget $external_prog{$ext_prog} -P $path");
        system("chmod -R 777 $path"); 
        system("tar -xvzf $path/*.tar.gz -C $path");
        chdir("$path");      
        my $program_name=`ls -d */`;
        chomp($program_name); 
        chdir("$program_name");
        system("./configure --prefix=$path");
        system("make");
        system("make install");
        chdir($Bin);

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $path/bin)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$path/bin\n";
           close PROFILE;
        }else{
           print "\n\nDeNoGAP did not recognize the shell you are using.  Sorry, but you will have to add SQLite to your path manually.\n\n";

           print "Add following line to your .bashrc or .cshrc file\n\n";

           print "export PATH=\$PATH:$path/bin\n\n";

           print "press [ENTER] to continue.\n\n";
           <>; 
        }
      };
      if($@){
         die($!);
      } 

     }
   }
   ##### INSTALL EMBOSS ####
    elsif($ext_prog eq 'emboss'){

       print "Do you want to Install EMBOSS on your system? ( Y or N): ";
       my $option_install=<STDIN>;
       chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
       }


     if($option_install eq "Y" or $option_install eq "y"){
         print "Installing EMBOSS requires root user permission\n\n";
         system($external_prog{$ext_prog});
     }
   }
   ##### INSTALL NCBI BLAST ####
   elsif($ext_prog eq 'blast'){

       print "Do you want to Install LOCAL NCBI BLAST on your system? ( Y or N): ";
       my $option_install=<STDIN>;
          chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
        }

       $HOME=$Bin;

       print "Installation path for NCBI-BLAST [$HOME/../exe]? (Y or N):";
       my $accept_path=<STDIN>;
       chomp($accept_path);

       SET_BLAST_PATH:

       if($accept_path eq "N"){
            print "Enter New installation path for NCBI-BLAST:";
            $HOME=<STDIN>;
            chomp($HOME);
       }elsif($accept_path eq "Y"){          
           $HOME="$HOME/../exe";

       }else{
           
           print "Invalid Option, Try Again\n";

           print "Installation path for NCBI-BLAST [$HOME/../exe]? (Y or N):";
           $accept_path=<STDIN>;
           chomp($accept_path);

           goto SET_BLAST_PATH;
       }

     if($option_install eq "Y" or $option_install eq "y"){

        my $path="$HOME/blast";

        print "The default path of installation for  NCBI BLAST is [$path], Give new path if you want to change: ";
        my $option_change_path=<STDIN>;
           chomp($option_change_path);

        if($option_change_path ne ''){            
           $path=$option_change_path; 
        }

        eval{
        mkdir($path);
        system("wget $external_prog{$ext_prog} -P $path");
        system("chmod -R 777 $path"); 
        system("tar -xvzf $path/*.tar.gz -C $path");   
        my $program_name=`ls -d */`;
        chomp($program_name); 
       
        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $path/$program_name/bin)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$path/$program_name/bin\n";
           close PROFILE;
        }else{
           print "\n\nDeNoGAP did not recognize the shell you are using.  Sorry, but you will have to add BLAST to your path manually.\n\n";

           print "Add following line to your .bashrc or .cshrc file\n\n";

           print "export PATH=\$PATH:$path/$program_name/bin\n\n";

           print "press [ENTER] to continue.\n\n";
           <>; 
        }
      };
      if($@){
         die($!);
      } 

     }
   }
   ####### INSTALL INTERPRO_SCAN ######
   elsif($ext_prog eq 'iprscan'){

       print "Do you want to Install LOCAL InterProScan on your system? ( Y or N): ";
       my $option_install=<STDIN>;
          chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
       }

       $HOME=$Bin;

       print "Installation path for InterProScan [$HOME/../exe]? (Y or N):";
       my $accept_path=<STDIN>;
       chomp($accept_path);

       SET_INTERPROSCAN_PATH:

       if($accept_path eq "N"){
            print "Enter New installation path for InterProScan:";
            $HOME=<STDIN>;
            chomp($HOME);
       }elsif($accept_path eq "Y"){          
           $HOME="$HOME/../exe";

       }else{
           
           print "Invalid Option, Try Again\n";

           print "Installation path for InterProScan [$HOME/../exe]? (Y or N):";
           $accept_path=<STDIN>;
           chomp($accept_path);

           goto SET_INTERPROSCAN_PATH;
       }

     if($option_install eq "Y" or $option_install eq "y"){

        my $path="$HOME/interproscan";

        print "The default path of installation for  InterProScan is [$path], Give new path if you want to change: ";
        my $option_change_path=<STDIN>;
           chomp($option_change_path);

        if($option_change_path ne ''){            
           $path=$option_change_path; 
        }

       eval{
        mkdir($path);
        chdir($path);
        system("wget $external_prog{$ext_prog}");
        system("wget $external_prog{iprscan_md}");
        system("md5sum -c interproscan-5.9-50.0-64-bit.tar.gz.md5");
        system("tar -pxvzf interproscan-5.9-50.0-64-bit.tar.gz");
        chdir($Bin);   

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $path/interproscan-5.9-50.0)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$path/interproscan-5.9-50.0\n";
           close PROFILE;
        }else{
           print "\n\nDeNoGAP did not recognize the shell you are using.  Sorry, but you will have to add InterProScan to your path manually.\n\n";

           print "Add following line to your .bashrc or .cshrc file\n\n";

           print "export PATH=\$PATH:$path/interproscan-5.9-50.0\n\n";

           print "press [ENTER] to continue.\n\n";
           <>; 
        }   
       };
       if($@){
         die($!);
       }
     }
   }     
}






