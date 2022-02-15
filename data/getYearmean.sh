#not needed for our study but maybe helpful

arr=(xmzka xmzkd xmzki xmzke xmzkg xmzkh xnagd xnagf xnagg xnage xnagh xnagb)
mkdir HadCM3/sea_ice/yearmean

for run in ${arr[*]}; do 
	cdo -yearmean HadCM3/sea_ice/sea_i"$run".nc HadCM3/sea_ice/yearmean/"$run"_yearmean.nc
done	
