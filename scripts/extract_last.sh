RED='\033[0;31m'
NC='\033[0m' # No Color
GREEN='\033[0;32m'
cd ../DATA
for d in */ ; do
    echo "$d looking"
    cd $d
    cp ../extract_last_config.py .
    python3 extract_last_config.py
    cd ../
done
