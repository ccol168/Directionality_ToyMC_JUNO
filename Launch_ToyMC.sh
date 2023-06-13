rm Output.*

if [ ! -e CfgExample.txt ]
then
    echo ERROR : No configuration file found
    echo Create a configuration file named CfgExample.txt in this directory
    exit 1
fi

echo Input the name of the run
read RunName

if [ -d "Events/$RunName" ]
then 
    echo WARNING: the directory $RunName already exists. This script will overwrite the old results to the old ones if the directory is not empty.
    while true; do
        read -p "Continue anyway?" yn
        case $yn in
            [Yy]* ) break;;
            [Nn]* ) exit 1;;
            * ) echo "Type y to continue, n to exit";;
        esac
    done
fi

if [ ! -d  "Events/$RunName" ] 
then
    mkdir Events/$RunName
    echo Directory Events/$RunName created
fi

cp CfgExample.txt Events/$RunName/

echo

./Directionality_ToyMC CfgExample.txt Events/$RunName/Output.root Events/$RunName/Output.txt

echo 
echo Launching Ordered_hits.py
echo

python Ordered_hits.py "Events/$RunName/Output.txt" "Events/$RunName/Nhits_prob.png" "Events/$RunName/Ang_dependence.png"