# For CS Class on Partial Differential Equations


To use:
```
git clone https://github.com/aclyde11/CS28515Proj1.git
cd CS28515Proj1
cmake .
make all
cd main
./main -f file.txt
python ../gif_maker.py file.txt
```

| Command | Description |
| --- | --- |
| -f | file name for point output |
| -n | Number of mesh points, default (default: 15) |
| -x_0 | Left end point over distance (default: 0) |
| -x_n | Right end point over distance (default: 3) |
| -init_value | For this genetics problem, set the starting value of U(x,o) |
| -tmax | Set max simulation time (default: 5) |