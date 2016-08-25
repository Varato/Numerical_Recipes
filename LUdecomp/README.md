####The LinearSolver executable can solve linear systems, i.e. find the solution of **A**x=**b**.

#####Example:

1. Edit test.txt file to construct data
```
#test.txt
2

1 1
1 -1

3
1
```
The “2” in line 1 means the coefficient matrix is 2x2. And the [3, 1] stands for "b" vector.

2. After compiling the relevant source files using ./compile.sh (one can also check inside), use command
```./LinearSolver test.txt``` to get the solution.







