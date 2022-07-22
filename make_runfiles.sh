RED='\033[0;31m'
NC='\033[0m' # No Color
GREEN='\033[0;32m'

for d in */ ; do
    echo "$d looking"
    cd $d
    if [ -d "./EAM_energy" ]
    then
	    cd "./EAM_energy"
	    cp ../../make_run.py .
	    python3 make_run.py
	    cd ../
    fi
    cd ../
done
