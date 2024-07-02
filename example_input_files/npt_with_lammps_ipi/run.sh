
if ! [ -f "simulation.restart" ]; then
python3 $ipi input_ipi_npt.xml > log.i-pi &
else
python3 $ipi simulation.restart > log.i-pi &
fi
  
sleep 20

mpirun -np $np $lmp_mpi < in1.lmp > /dev/null &
python3 run-ase.py &

wait
