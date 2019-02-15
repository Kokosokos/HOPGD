set -e
if [ $# -le 1 ]
then
        echo "ERROR: Please provide start and finish folder numbers"
        exit
fi
start=$1
fin=$2
for i in $(seq $start $fin)
do
if [ ! -f $i/done ]
then 
	echo "-----------------" $i "---------------------"
	cd $i
#	../../run_abaqus.sh just.inp
	abaqus Job="just" interactive
	cd ..
	echo "y" > $i/done
	echo "=================" $i "====================="
fi
done
