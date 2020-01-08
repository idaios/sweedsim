name="continent_island_twoway_bottleneck4"
basedir=`pwd`
demographyOptions=" -en 0.01 1 0.05 -en 0.05 1 1 -en 0.01 2 0.05 -en 0.05 2 1 -ej 3 2 1 "
sweedfolder="../sweed/"
trajfolder="../trajdemog/"
reps=1000
mkdir $name

cd $name

echo "SIMS START"

$basedir/$trajfolder/trajdemognpops -nreps $reps -npop 2 -s 0.04 0.0 -t 0 1500 -pfinal 0.99 0. -eps .25 1.1 -seed 30341 -npres 10000 -mig 0.002 0.002 | stepftn2 2 > trajectory_$name.txt

echo "SWEED STARTS"

$basedir/$sweedfolder/SweeD -name $name -grid 10 -input $basedir/mssel.data.out -length 100000 -strictPolymorphic -mssel 20 $reps 0 20 ./trajectory_$name.txt 5000 -I 2 0 20 0 0 8  -t 8000 -r 8000 10000 $demographyOptions -ej 3 2 1

cp ../run.sh .

cd $basedir
