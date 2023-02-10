# Topology Optimization (MATLAB 88 Lines Code)
A topology optimization project that has the good coding style based on the [reference 1](https://link.springer.com/article/10.1007/s00158-010-0594-7) and [reference 2](http://paulino.princeton.edu/courses/cee307/2016/resources/papers/IJNME_MMA1987.pdf).

## 1. Run The Project
execute `src/main.m`

## 2. User Parameter Settings
In `src/main.m`, it has a parameters detail comment.

If you want to use other scenarios, you can add a new function file @ `assets/samples/`, then call that @ `src/main.m`

Also, you can modify the following files if you want to change other parameters:
- The `MMA (Method of Moving Asymptotes)` update method: `src/update_methods/mma/*`
- The `Projection Function`: `src/utils/projection_function.m`
- The `Filter`: `src/utils/init_filter` and `src/utils/filter_.m`
- The `Stiffness Matrix`: `src/utils/init_stiffness_matrix.m`

## 3. References
- [Efficient Topology Optimization in MATLAB using 88 lines of code](https://link.springer.com/article/10.1007/s00158-010-0594-7)
- [The Method of Moving Asymptotes -- A New Method For Structural Optimization](http://paulino.princeton.edu/courses/cee307/2016/resources/papers/IJNME_MMA1987.pdf)