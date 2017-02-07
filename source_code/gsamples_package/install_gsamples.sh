#!/bin/bash


#To install gsample, you need to first download/compile two software packages, 
#Samtools and Boost, setup the corresponding environmental variables, and then
#compile the source code of gsample.

if [ $# -lt 1 ]
then
	echo "Usage: ./install_gsamples.sh <-v verbose> "
	echo "        [-s samtool_version] [-b boost_version] "
	echo "        [-k remove_temporary_files] [-K remove_temporary_folders] "
	echo "[note1]: Set 'verbose' as 1 for verbose compile, 0 for not verbose "
	echo "         if set to 2, then use 'full verbose mode' "
	echo "[note2]: Default value of 'samtool_version' is 1.3.1 "
	echo "         Default value of 'boost_version' is 1.63.0 "
	echo "[note3]: Default value of 'remove_temporary_files' is 1 (remove) "
	echo "         Default value of 'remove_temporary_folders' is 1 (remove) "
        exit
fi
curdir="$(pwd)"


# ----- get arguments ----- #
samtool_version="1.3.1"
boost_version="1.63.0"
remove_temporary_files=1
remove_temporary_folders=1
verbose=1

#-> parse arguments
while getopts ":s:b:k:K:v:" opt;
do
	case $opt in
	#-> optional arguments
	s)
		samtool_version=$OPTARG
		;;
	b)
		boost_version=$OPTARG
		;;
	k)
		remove_temporary_files=$OPTARG
		;;
	K)
		remove_temporary_folders=$OPTARG
		;;
	v)
		verbose=$OPTARG
		;;
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done


#---------
#-> Part I
# Samtools
#---------

#1, Download Samtools from {http://www.htslib.org/} with version 1.2 or higher.
#Compile it to generate the htslib file libhts.a. Set environment variable
#HTSLIB to indicate the directory of libhts.a.  For example, for Unix platforms,
#add the following statement to the file ~/.bash_profile:
#        export HTSLIB="/directory/to/your/htslib/htslib-1.2.1"

cd $curdir
SAMtool_Version=$samtool_version
SAMtool_Root="samtools-"$SAMtool_Version
HTSlib_Root="htslib-"$SAMtool_Version
echo ""
echo "Part I: Download and install Samtools $SAMtool_Version "
echo ""
if [ ! -d $SAMtool_Root ] || [ ! -f $SAMtool_Root/$HTSlib_Root/libhts.a ]
then
	rm -f $SAMtool_Root.tar.bz2*
	echo "1.1 download samtool package"
	if [ $verbose -eq 1 ]
	then
		wget https://github.com/samtools/samtools/releases/download/$SAMtool_Version/$SAMtool_Root.tar.bz2
	else
		wget -q https://github.com/samtools/samtools/releases/download/$SAMtool_Version/$SAMtool_Root.tar.bz2
	fi
	echo "1.2 uncompress samtool package"
	if [ $verbose -eq 2 ]
	then
		tar vxjf $SAMtool_Root.tar.bz2
	else
		tar xjf $SAMtool_Root.tar.bz2
	fi
	rm -f $SAMtool_Root.tar.bz2
	echo "1.3 compile samtool package"
	cd $SAMtool_Root
		if [ $verbose -eq 2 ]
		then
			make
		else
			make 1> samtool_out1 2> samtool_out2
		fi
	cd ../
fi
export HTSLIB=$curdir/$SAMtool_Root/$HTSlib_Root


#----------
#-> Part II
# Boost
#----------

#2, Download Boost from {http://www.boost.org}. Uncompress it
#somewhere~(compiling and installing are not necessary). Set environment
#variable BOOST_HOME to indicate the directory of Boost.  For example, for Unix
#platforms, add the following statement to the file ~/.bash_profile:
#	export BOOST_HOME="/directory/to/your/boost/boost_1_60_0"

cd $curdir
Boost_Version_=$boost_version
Boost_Version=`echo $Boost_Version_ | sed 's/\./_/g'`
Boost_Root="boost_"$Boost_Version
echo ""
echo "Part II: Download Boost $Boost_Version_ "
echo ""
if [ ! -d $Boost_Root ]
then
	rm -f $Boost_Root.tar.gz*
	echo "2.1 download boost library"
	if [ $verbose -eq 1 ]
	then
		wget https://sourceforge.net/projects/boost/files/boost/$Boost_Version_/$Boost_Root.tar.gz
	else
		wget -q https://sourceforge.net/projects/boost/files/boost/$Boost_Version_/$Boost_Root.tar.gz
	fi
	echo "2.2 uncompress boost library"
	if [ $verbose -eq 2 ]
	then
		tar vxzf $Boost_Root.tar.gz
	else
		tar xzf $Boost_Root.tar.gz
	fi
	rm -f $Boost_Root.tar.gz
fi
export BOOST_HOME=$curdir/$Boost_Root


#-----------
#-> Part III
# gsamples
#-----------

#3, Execute the following commands to generate Makefile and compile:\\
#The executable file gsamples will be present at src.

cd $curdir
echo ""
echo "Part III: Install gsamples "
echo ""
if [ ! -f "gsamples" ]
then
	echo "3.1 automatic configuration"
	aclocal
	autoconf
	autoheader
	automake -a
	if [ $verbose -eq 1 ]
	then
		echo "3.2 automatic configuration"
		./configure 
		echo "3.3 compile gsamples"
		make 
	else
		echo "3.2 automatic configuration"
		./configure 1> configure_out1 2> configure_out2
		echo "3.3 compile gsamples"
		make 1> make_out1 2> make_out2
	fi
	mv src/gsamples .
fi

#----------
#-> Part IV
# clean
#----------

#-> 4.1 clean temporary files in 'src/'
cd $curdir
echo ""
echo "Part IV: remove temporary files and folders "
echo ""
echo "4.1 remove compiled files in src/"
cd src/
	make clean 1> ws1 2> ws2
	rm -f ws1 ws2
	rm -f Makefile.in
	rm -f Makefile
	rm -f gmon.out
cd ../

#-> 4.2 clean temporary files and folders in current directory
if [ $remove_temporary_files -eq 1 ]
then
echo "4.2 remove temporary files in current directory"
rm -f aclocal.m4
rm -f config.h.in
rm missing
rm install-sh
rm depcomp
rm -rf autom4te.cache
rm -f configure
rm -f configure_out1 configure_out2
rm -f Makefile.in
rm -f config.status
rm -f Makefile
rm -f make_out1 make_out2
rm -f config.h
rm -f stamp-h1
rm -f config.log
rm -f gmon.out
fi

#-> 4.3 remove SAMtool and Boost
if [ $remove_temporary_folders -eq 1 ]
then
echo "4.3 remove temporary folders in current directory"
rm -rf $SAMtool_Root
rm -rf $Boost_Root 
fi

#----------
#-> Part V
# command
#----------

#4, COMMAND LINE
#The usage of gsamples is: ./gsamples <in-bam-file> <out-sample-file> [-f fasta-dir] [-g gtf-file] [-e min-expression]
# The parameter <in-bam-file> specifies the input reads alignemnt file. It has to be sorted (for example, by samtools).
# The parameter <out-sample-file> specifies the the output sample file.
# The parameter [fasta-dir] is optional but strongly recommended, which specifies directory of all sequences files of all chromosomes. 
# These sequence files will be used as features.
# The parameter [gtf-file] is optional, which specifies the ground-truth expressed transcripts. If this file is given, the true boundaries
# will be contained in the output sample file as true labels; otherwise, the true labels will be -1. When evaluting mode is on (i.e., to
# evalute the performance of this method, this file should be given.
# The parameter [min-expression] is together with [-g gtf-file], which will make gsamples ignore these transcripts in the gtf-file with
# expression abundance under this value.

