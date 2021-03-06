name="continent_island_twoway_bottleneck3"
basedir=`pwd`
demographyOptions=" -en 0.01 1 0.005 -en 0.05 1 1 -en 0.01 2 0.005 -en 0.05 2 1 -ej 3 2 1 "
mkdir $name
cd $name

echo "SIMS START"

trajdemognpops -nreps 5000 -npop 2 -s 0.04 0.0 -t 0 1500 -pfinal 0.99 0. -eps .25 1.1 -seed 30341 -npres 10000 -mig 0.002 0.002 | stepftn2 2 > trajectory_$name.txt

echo "SWEED STARTS"

$basedir/SweeD -name $name -grid 10 -input $basedir/mssel_mig2.out -length 100000 -strictPolymorphic -mssel 20 2000 0 20 ./trajectory_$name.txt 5000 -I 2 0 20 0 0 8  -t 8000 -r 8000 10000 $demographyOptions -ej 3 2 1

cp ../run.sh .

cd $basedir
