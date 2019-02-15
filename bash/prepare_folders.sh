#
#exec bash
loop()
{

		if [ $# -eq 0 ]
		  then
			export dimMax="${#arr[@]}" #array lenght
			d=1
			local args=()
		else
#			let dim=$1 #Loop over this dimension
			let d=$1+1
			local args=($@)
			args=${args[@]:1:${#args[@]}}
		fi

		#let newdim=dim-1
		local dim=0
		let dim=dimMax-d
		local newdim=0
        let newdim=dimMax-d-1
		if [ "$dim" -ge 0 ]
		then
			local r=(${arr[$dim]})

			for i in ${r[@]}
			do
				local params=($args)
				params+=($i)
				loop $d ${params[@]}

			done
		else
			let k=k+1
			mkdir $k
			params=($args)
			#reverse params array and print
			rparams=()
			for (( idx=${#params[@]}-1 ; idx>=0 ; idx-- ))
			do
				rparams+=(${params[$idx]})
			done
			echo ${rparams[@]} > $k"/params.dat"	
		fi
		d=$1
}

inpName=$1 #template .inp filename
echo "Preparation started"
if [ $# -le 1 ]
then
        echo "ERROR: Please provide .inp template file and parameters file"
        exit
fi
parametersFile=$2
readarray arr < $parametersFile
export k=0
loop  
echo "Total number of folders created: " $k
echo "Putting params.dat in folders"
for d in */
do
        i=1
        cd $d
        cp ../128320_4params_template.inp just.inp
        p=$(cat params.dat)
        p=($p)
        for p1 in ${p[@]};
                do
                sed -i "s/__param_${i}_Value__/$p1/g" just.inp
                let i=i+1
        done 
        cd ..
done
echo "Preparation finished"

