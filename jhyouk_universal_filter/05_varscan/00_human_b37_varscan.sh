pwd_mpileup=$1
tumor=$2
normal=$3
log=$2.varscan2.log

set -e

echo "varscan2 start for" $2 > $log
java -Xmx12g -jar /home/users/tools/varscan2.4.2/VarScan.v2.4.2.jar somatic $1/$3.mpileup $1/$2.mpileup $2.varscan --min-var-freq 0.01 --output-vcf 1 &>> $log
echo "varscan2 finish for" $2 >> $log
mv $log $2.varscan2.success
