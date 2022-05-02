func="python optimization.py "
betas=1.0
betat=1.0
problem="pnm4"
tmux new-session -d 

#solver="scipy_neldermead "
solver="lbfgs"
i=$((0))
#INPUT=parameter_settings.csv
INPUT=$1
OLDIFS=$IFS
IFS=','
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read rho lc ls lt dor alpha empty
do
	if ((i > 0)); then
		tmux split-window -h
		tmux send-keys -t $i "$func $dor $alpha $betas $betat $problem $solver $rho $lc $ls $lt" C-m
		tmux select-layout tile
        	echo "rho : $rho"
        	echo "lc : $lc"
        	echo "ls : $ls"
        	echo "lt : $lt"
        	echo "dor : $dor"
        	echo "alpha : $alpha"
	fi
	i=$((i+1))
done < $INPUT
IFS=$OLDIFS

