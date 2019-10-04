make
#exec > logfile.csv
cd  /gpfs/alpine/proj-shared/csc289/trianglecounting/snap
for d in */; do
	#cd output
	#pwd	
	touch output/${d::-1}".csv"
	
	exec > output/${d::-1}".csv" 
	./prof /gpfs/alpine/proj-shared/csc289/trianglecounting/snap/$d
done
