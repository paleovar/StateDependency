#go to the directory of hadcm3 data first that contains the `./surface_temperature folder` using cd ...
#create directory "/stationary/"
#detrend and deseasonalize the runs and write them to the new directory

arr=(xmzka xmzkd xmzki xmzke xmzkg xmzkh xnagd xnage xnagf xnagg xnagh xnagb)

#temperature
mkdir HadCM3/surface_temperature/stationary
for run in ${arr[*]}; do 
	cdo detrend HadCM3/surface_temperature/"$run".nc tmp_"$run".nc
	cdo ymonsub tmp_"$run".nc -ymonmean tmp_"$run".nc HadCM3/surface_temperature/stationary/"$run".nc
	rm tmp*
done

#sea ice
mkdir HadCM3/sea_ice/stationary
for run in ${arr[*]}; do 
	cdo detrend HadCM3/sea_ice/"$run".nc tmp_"$run".nc
	cdo ymonsub tmp_"$run".nc -ymonmean tmp_"$run".nc HadCM3/sea_ice/stationary/"$run".nc
	rm tmp*
done

#precipitation
mkdir HadCM3/precipitation/stationary
for run in ${arr[*]}; do 
	cdo detrend HadCM3/precipitation/"$run".nc tmp_"$run".nc
	cdo ymonsub tmp_"$run".nc -ymonmean tmp_"$run".nc HadCM3/precipitation/stationary/"$run".nc
	rm tmp*
done

#Not needed for our study but maybe also hhelpful:

#sea level pressure
#mkdir HadCM3/sea_level_pressure/stationary
#for run in ${arr[*]}; do 
#	cdo detrend HadCM3/sea_level_pressure/"$run".nc tmp_"$run".nc
#	cdo ymonsub tmp_"$run".nc -ymonmean tmp_"$run".nc HadCM3/precipitation/stationary/"$run".nc
#	rm tmp*
#done

#wind
#mkdir HadCM3/wind/stationary
#for run in ${arr[*]}; do 
#	cdo detrend HadCM3/wind/"$run".nc tmp_"$run".nc
#	cdo ymonsub tmp_"$run".nc -ymonmean tmp_"$run".nc HadCM3/wind/stationary/"$run".nc
#	rm tmp*
#done
