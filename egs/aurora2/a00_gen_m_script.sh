# training
echo "addpath(genpath('masking'));" > a00_check_data_train.m
find ../wav/train/* | grep "\.wav" | sed '/\/clean\//d' | awk '{print "StereoDataCheck(@" $1 "@, @" $1 "@);"}' | sed "s: @../wav/train/[a-z0-9_-]*/: @../wav/train/clean/:g" | sed "s:@:':g" >> a00_check_data_train.m

# testing
for data in testa testb
do
	if [ -e a00_check_data_${data}.m ]
    then
        rm a00_check_data_${data}.m
		echo "addpath(genpath('masking'));" > a00_check_data_${data}.m
    fi

	for id in 1 2 3 4
	do
		for snr in 20 15 10 5 0 -5
		do
			ls ../wav/${data}/clean${id}/*.wav | sed "s:../wav/${data}/clean${id}/::g" | awk -v data=${data} -v id=${id} -v snr=${snr} '{print "StereoDataCheck(@../wav/" data "/clean" id "/" $1 "@, @../wav/" data "/n" id "_snr" snr "/" $1 "@);"}' | sed "s:@:':g" >> a00_check_data_${data}.m
		done
	done
done

# testc
data=testc
if [ -e a00_check_data_${data}.m ]
then
    rm a00_check_data_${data}.m
	echo "addpath(genpath('masking'));" > a00_check_data_${data}.m
fi

for id in 1 2
do	
	for snr in 20 15 10 5 0 -5
	do
		ls ../wav/${data}/clean${id}/*.wav | sed "s:../wav/${data}/clean${id}/::g" | awk -v data=${data} -v id=${id} -v snr=${snr} '{print "StereoDataCheck(@../wav/" data "/clean" id "/" $1 "@, @../wav/" data "/n" id "_snr" snr "/" $1 "@);"}' | sed "s:@:':g" >> a00_check_data_${data}.m
	done
done

