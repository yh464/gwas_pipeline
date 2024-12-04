while true;
do
echo Failed subjects:;
find . -size -296c -size +120c;
echo Timed out subjects:;
find . -size -120c -size +100c;
echo Jobs waiting:;
squeue -u yh464 | wc -l;
echo ;
echo Jobs completed:;
ls | wc -l;
echo;
echo Successful jobs:;
find . -size +296c | wc -l;
echo;
mybalance;
sleep 5m;
done
