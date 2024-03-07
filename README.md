# PySSTraj

A simple Python script, based on Pandas and MDTraj, for tracking the temporal evolution of secondary structure along a MD trajectory.

```
Arguments
    -p    --top     Topology file. If you want a result per chain, use a pdb file.
    -t    --traj    Trajectory file. Output from MD simulation
    -s    --str     Stride. Value must be an integer.
    -d    --dir     Working directory. Default is current directory.
    -v    --ver     Script version
```

Example:
python PySSTraj.py -p run.gro -t traj_prod.xtc -s 1
