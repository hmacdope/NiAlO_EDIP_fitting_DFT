for d in */ ; do
    echo "$d extracting"
    cd $d
    extract_vasp
    cd ../
done
