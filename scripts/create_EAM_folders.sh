cd ../DATA/

for d in */ ; do
    echo "$d looking"
    cd $d
    if [ ! -d "./EAM_energy" ]
    then
	    mkdir EAM_energy
	    cd "./EAM_energy"
	    echo "made EAM folder"
	    pwd
	    cp ../*.data .
	    cp /store/Terry_project/EAM_energy_base/* .
	    cd ../
    fi
    cd ../
done
