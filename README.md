# For CS Class on Partial Differential Equations

This project includes my **numerical library** 
(numerlicallib) supporting vector and tri-diagonal matrix operations and
and **pde solver library** (pdesolverslib). proj1 is for the heat equation and proj2 is
for the wave equation.

###

To use:
```
git clone https://github.com/aclyde11/CS28515Proj1.git
cd CS28515Proj1
cmake .
make all
cd main
./proj1 -f file.txt
./proj2 -f file.txt
python ../gif_maker.py file.txt //produces gif modeling pde
```

| Command | Description |
| --- | --- |
| -f | file name for point output |
| -n | Number of mesh points, default (default: 15) |
| -x_0 | Left end point over distance (default: 0) |
| -x_n | Right end point over distance (default: 3) |
| -init_value | For this genetics problem, set the starting value of U(x,o) |
| -tmax | Set max simulation time (default: 5) |
| -c0 | Set mass matrix value |
| -k0 | Set stiffness matrix value |
| -dt | Set dt |
| -periods | set number of periods |