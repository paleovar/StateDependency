#not needed for our study but maybe helpful

arr=(xmzka xmzkd xmzki xmzke xmzkg xmzkh xnagd xnagf xnagg xnage xnagh xnagb)
mkdir HadCM3/sea_ice/fldmean

for run in ${arr[*]}; do 
	cdo -fldmean HadCM3/sea_ice/"$run".nc HadCM3/sea_ice/fldmean/"$run"_fldmean.nc
done	
