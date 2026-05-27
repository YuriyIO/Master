#local
./mtx2csr /mnt/c/Programs/C/VS_Code/SLAY/project2/Matrix/494_bus.mtx /mnt/c/Programs/C/VS_Code/SLAY/project2/src/CSR/test2

#polus

wget https://suitesparse-collection-website.herokuapp.com/MM/DIMACS10/hugetric-00000.tar.gz
pscp  -i "C:\Users\79304\.ssh\SKI\YuriiPovzhik.ppk"  C:\Programs\C\VS_Code\SLAY\project2\Matrix\hugetric-00000.tar.gz edu-cmc-sqi23-25@polus.cs.msu.ru:/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/programs/project/Matrix
tar -xvzf hugetric-00000.tar.gz
cd hugetric-00000
cp hugetric-00000.mtx /home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/programs/project/Matrix

awk '
/^%%MatrixMarket/ { gsub("pattern", "real"); print; next } 
/^%/ { print; next } 
NF == 2 { print $1, $2, "1.0"; next } 
{ print }
' /home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/programs/project/Matrix/hugetric_00000.mtx > /home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/programs/project/Matrix/hugetric_00000_weighted.mtx

./mtx2csr /home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/programs/project/Matrix/hugetric_00000_weighted.mtx  /home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/programs/project/src/CSR/hugetric_00000_weighted


