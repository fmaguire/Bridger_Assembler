[![Build Status](https://travis-ci.org/fmaguire/Bridger_Assembler.svg?branch=master)](https://travis-ci.org/fmaguire/Bridger_Assembler)

```
Rehost by Finlay Maguire to fix some minor bugs - 
All credit goes to original authors below;



================
 Description
================
                                                                                  
 Bridger is an efficient de novo trascriptome assembler for RNA-Seq data. It can assemble 
 all transcripts from short reads (single or paired) without using a reference.  
 The software expects as input RNA-Seq reads in fasta or fastq format, and ouput all
 assembled candidate transcripts in fasta format. Briefly, it works in two step: first,
 Bridger partitions the sequence data into many individual splicing graphs, each capturing 
 the full transcriptional complexity at a given gene or no more than a few genes. Then, 
 Bridger uses a rigorous mathematical model called minimum path cover to search minimal 
 set of paths(transcripts) that can be supported by our data and could explain all observed 
 splicing events of each locus.

 This software is free to use, modify, redistribute without any restrictions, 
 except including the license provided with the distribution. 


================
 Installation
================


 1. Installing Boost

    a) download latest boost and unpack it.

       $ tar zxvf boost_1_47_0.tar.gz

    b) change to the boost directory and run ./bootstrap.sh.

       $ cd  boost_1_47_0
       $ ./bootstrap.sh

    c) Run ./b2 and set the install path. 

       $ ./b2 install --prefix=<YOUR_BOOST_INSTALL_DIRECTORY>

       For example,if you want install boost in /home/czheng/local/boost,the commnd is :
       $ ./b2 install --prefix=/home/czheng/local/boost
     
       If the boost is installed successfully, you would fild two sub-directories in /home/czheng/local/boost/:
       /home/czheng/local/boost/include/ 
       /home/czheng/local/boost/lib/

      Note: The default Boost installation directory is /usr/local. Take note of the boost installation
       directory, beacuase you need to tell the Bridger installer where to find boost later on.

    d) Set the LD_LIBRARY_PATH enviroment variable:
       
       The ~/.bash_profile ($HOME/.bash_profile) or ~/.profile file is executed when you login using console or remotely using ssh. 
       Append the following command to  ~/.bash_profile or ~/.profile file:
       $ export LD_LIBRARY_PATH=/home/czheng/local/boost/lib:$LD_LIBRARY_PATH

       Save and close the file.

       OR
       
       just type the command:
       $ export LD_LIBRARY_PATH=/home/czheng/local/boost/lib:$LD_LIBRARY_PATH  

       Note: please replace "/home/czheng/local/boost/lib" with your own directory "<YOUR_BOOST_INSTALL_DIRECTORY>/lib"
       If you do not set this variable , you would possible see the follwoing error information:
        "error while loading shared libraries: libboost_serialization.so.1.47.0: cannot open shared object file: No such file or directory"

 2. Building Bridger [Make sure Boost has been installed successfully]
  
    a) Unpack the Bridger and change to the Bridger direcotry.

       $ tar zxvf Bridger_r2013-06-02.tar.gz
       $ cd Bridger_r2013-06-02

    b) Configure Bridger. If Boost is installed somewhere other than /usr/local, you will need to tell
       the installer where to find it using --with-boost option.

       $ ./configure --with-boost=/home/czheng/local/boost/
       Note: please replace "/home/czheng/local/boost/" with your own directory "<YOUR_BOOST_INSTALL_DIRECTORY>"
     

    c) Make Bridger.

       $ make
      
      note: If you build boost suffessfully without using --prefix option, the following commands may need before "make":
          export LIBS="-L/home/czheng/boost_1_47_0/stage/lib" (replace "/home/czheng/boost_1_47_0/" with your own directory)
          export CPPFLAGS="-I/home/czheng/boost_1_47_0/"

 3. Test the installation. Test data are provided with sofeware distribution in the sample_test directory.
     $ cd src
     $ ./Assemble -h
     you would see:
     ===============================================================================
     Usage: Assemble [--reads/--kmers] <filename>  [opts] 
     ===============================================================================
     **Required :
     --reads/-i <string>           : the name of the file containing reads

     ** Optional :
     --kmer_length/-k <int>        : length of kmer, default: 25.
     --double_stranded_mode        : set it true if double stranded mode.
     --fr_strand<int>              : strand specific protocol, default: 1 
                                       ( 1 : fr-firststrand, e.g. dUTP, NSR, NNSR 
                                         2 : fr-secondstrand, e.g. Strandard SOLID ) 
     --paired_end                  : set it true if paired reads.
     --min_seed_coverage <int>     : minimum coverage of seed kmer, default: 2.
     --min_seed_entropy <float>    : minimum entropy of seed kmer, default: 1.5
     --min_kmer_coverage <int>     : minimum coverage of kmer used to extend, default: 1.
     --min_kmer_entropy <float>    : minimum entroy of kmer used to extend, default: 0.0
     --min_junction_coverage <int> : minimum of the coverage of a junction, default: 2.
     --min_ratio_non_error <float> : min ratio for low/high alternative extension that is 
                                          not an error, default: 0.05.
     --pair_gap_length             : gap length of paired reads, default: 200.
     --out_dir/-o <string>         : name of directory for output, default : ./RawGraphs 
     --help/-h                     : display the help information.

    ===============================================================================
     

     Note : If you  see the error information like "error while loading shared libraries: libboost_serialization.so.1.47.0: 
               cannot open shared object file: No such file or directory", please set the LD_LIBRARY_PATH variable by command:
              "export LD_LIBRARY_PATH=/home/czheng/local/boost_1_47_0/lib:$LD_LIBRARY_PATH"
 
     Test Bridger with a small data:

     $ cd sample_test/
     $ ./run_Me.sh



===============
 Uasge
===============

 A typical command maybe like this:
 
 ./Bridger.pl --seqType fq --left reads.left.fq --right reads.right.fq --CPU 6
 
 Note:
   "--SS_lib_type" is recommended to be used for strand-specific RNA-Seq data.
 
  For more inforamtion, use --help option. 
  or, visit : http://BridgerRNASeq.sourceforge.net


================
 Changelog
================

 Version r2014-12-01
 - latest version


===============
 Authors
===============

 Zheng Chang designed and wrote Bridger, with substantial technical input
 from Yu Zhang and Cody Ashby.  


================
 Contact
================

 Any questions, problems, bugs are welcome and should be dumped to
 Zheng Chang <changzmaths@gmail.com>

 Created on June 25, 2012.
``` 
