for FOLDER in $(ls ~/dev/stori/data/new_housekeeping)
do
    echo $FOLDER;
    cd ~/dev/stori/data/new_housekeeping/$FOLDER;
    ls;
    ~/dev/stori/run_analysis.sh;
    cd ~/dev/stori/;
done;