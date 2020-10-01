#!/bin/bash
     
#SBATCH -A kmane
#SBATCH -J TrialOnOneXme
#SBATCH -p standard
#SBATCH -w cpn6
#SBATCH --mem=100000
#SBATCH --error="TrialOnOneXme.error"
#SBATCH --output="TrialOnOneXme.out"    
#SBATCH --mail-type=BEGIN,FAIL,END          
#SBATCH --mail-user=kmane@mrc.gm

typeset -F SECONDS
module load R/4.0.2
module load gnu7/7.3.0

function usage 
{
	 echo -e "\nAbout:\t this function is used to infer IBD from biallelic
	 SNPs data."
	 echo -e "\t Note that the group file is used when inferring iR. It 
	 should have 3 columns: fid,iid,group."
	 echo -e "\n"
     printf "Usage:\t %s: <options> in interactive mode   \n" $0  #-s <value>
	 echo -e "\t sbatch Run_isoRelate.sh <options> in batch mode"
	 echo -e "\n"
	 echo -e "Options: \t\n"
     echo -e "\t-p : full path to the ped file                                                  "
     echo -e "\t-m : full path to the map file                                                  "
	 echo -e "\t-o : full path to the output files directory                                    "
     echo -e "\t-c : full path to the folder where the R codes are stored                       "
     echo -e "\t-g : full path to the group file                                                "
	 echo -e "\t-h : to display the help                                                                                       "
	 
     echo -e "\n"
     exit 2
}

function validate
{
	if [[ -f $p ]]
	then
		inform "The ped file is "$p
	else
		problem "required ped file not found! Use -h option to display the help."
        exit
	fi
	
	if [[ -f $m ]]
	then
		inform "The map file is "$m
	else
		problem "required map file not found! Use -h option to display the help."
        exit
	fi
	
	if [[ -d $o ]]
	then
		inform "The output folder is "$o
	else
		problem "required output directory not found! Use -h option to display the help."
        exit
	fi
	
	if [[ -d $c ]]
	then
		inform "The path to the codes folder is "$c
	else
		problem "required codes directory not found! Use -h option to display the help."
        exit
	fi
	
	if [[ -f $g ]]
	then
		inform "The group file is "$g
	else
		problem "required group file not found! Use -h option to display the help."
        exit
	fi
	
}

function inform 
{
  TIMESTAMP=`date +"%a %b %d %X %Y"`
  echo "[$TIMESTAMP] $*" 1>&2
}

function problem 
{
  TIMESTAMP=`date +"%a %b %d %T %Y"`
  echo -e "\n[$TIMESTAMP] *ERROR*: $*" 1>&2
}

while getopts hp:m:o:c:g: name
do
    case $name in
        p) p="$OPTARG";;
		m) m="$OPTARG";;
		o) o="$OPTARG";;
		c) c="$OPTARG";;
		g) g="$OPTARG";;
		h) usage;;
        ?) usage;;
    esac
done
shift $(($OPTIND - 1))
mkdir -p $o
validate

workingDir=`pwd`
${workingDir}/Run_isoRelate.R $p $m $o $c $g



echo -e "\n\nWall time is "$SECONDS" seconds"







