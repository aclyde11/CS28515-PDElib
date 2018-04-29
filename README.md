# For CS Class on Partial Differential Equations

This project includes my **numerical library** 
(numerlicallib) supporting vector and tri-diagonal matrix operations and
and **pde solver library** (pdesolverslib), currently supporting equations of the form 
c(x)u_x - (k(x)u_x)_x = F(x,u), such as the heat equation and logisitic function.

###

To use:
```
git clone https://github.com/aclyde11/CS28515Proj1.git
cd CS28515Proj1
cmake .
make all
cd main
./main -f file.txt
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