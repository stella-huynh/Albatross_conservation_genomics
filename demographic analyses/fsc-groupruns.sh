#!/bin/bash



###############################################################

display_help() {
 echo -e "\nUsage: $0 [OPTIONS...] \n"
 echo "   -m/--model        Name of demographic model tested."
 echo "   -d/--directory    Directory containing the \".bestlhoods\" files."
 echo "   -o/--outdir       Output directory."
 echo "   -s/--suffix       Suffix to append output file and prefix [Optional]."
 echo "   -h/--help         Display help."
 echo -e "\n" && exit 1
}

###############################################################


POSITIONAL=()

while [[ $# -gt 0 ]]
do
      key="$1"
      case $key in
              -m|--model) #Name of demographic model
              MODEL="$2"
              shift # past argument
              shift # past value
              ;;
              -d|--directory) #Directory containing .bestlhoods file
              wDIR="$2"
              shift # past argument
              shift # past value
              ;;
              -o|--outdir) #Output directory
              oDIR="$2"
              shift
              shift
              ;;
              -s|--suffix) #Suffix to append to filename and prefix
              SUFFIX="$2"
              shift
              shift
              ;;
              -h|--help) #display HELP
              display_help #call funtion to display script commands
              exit 0
              ;;
              --default)
              DEFAULT=YES
              shift # past argument
              ;;
              *) # unknown option
              POSITIONAL+=("$1") # save it in an array for later
              shift # past argument
              ;;
      esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters


###############################################################


if [[ ! -z $SUFFIX ]]; then SUFFIX="_$SUFFIX" ; fi

if [[ -z $oDIR ]]; then
        oDIR="${wDIR}/sum_bestlhoods" && mkdir -p -m775 ${oDIR}
fi

listfiles=$( find $wDIR -type f -not -path "*sum_bestlhoods*" -name "*.bestlhoods" -print | sort -V )
counter=0


for file in $listfiles
do

   PREFIX=$( basename ${file%.bestlhoods}${SUFFIX} )
   i=$( echo $file | grep -oP '(?<=run)[0-9]+' )
   
   #cp $file ${oDIR}/${MODEL}_${PREFIX}_run${i}.bestlhoods
   

   if [[ $counter == 0 ]]; then
        awk -v model="$MODEL" -v group="$PREFIX" \
        '{OFS="\t"} NR==1 {print "MODEL","GROUP",$0} NR>1 {print model,group,$0}' \
                $file \
                > ${oDIR}/${MODEL}_${PREFIX}_all.bestlhoods
   else
        awk -v model="$MODEL" -v group="$PREFIX" \
        '{OFS="\t"} NR>1 {print model,group,$0}' \
                $file \
                >> ${oDIR}/${MODEL}_${PREFIX}_all.bestlhoods
   fi

   counter=$((counter +1))
   
done


