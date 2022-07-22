RED='\033[0;31m'
NC='\033[0m' # No Color
GREEN='\033[0;32m'

for d in */ ; do
    echo "$d looking"
    cd $d
    if [ -d "./EAM_energy" ]
    then
	    cd "./EAM_energy"
	    echo "doing EAM for"
	    pwd
	    lmp -in single_point.run > single_point.out
	    if [ $? -ne 0 ]
	    then
		    echo -e "${RED}FAILED${NC}"
	    else
		    echo -e "${GREEN}DONE${NC}"
	    fi
	    cd ../
    fi
    cd ../
done
