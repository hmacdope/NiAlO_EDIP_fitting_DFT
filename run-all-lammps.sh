for d in */ ; do
    echo "$d looking"
    cd $d
    if [ -d "./EAM_energy" ]
    then
	    cd "./EAM_energy"
	    echo "doing EAM for"
	    pwd
	    lmp -in single_point.run > single_point.out
	    cd ../
    fi
    cd ../
done
