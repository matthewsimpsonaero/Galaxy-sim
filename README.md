# Galaxy Simulation Using MATLAB 

### N body gravitational simulation using a 4th order runge kutta method for numerical integration with minimal introduction of error optimized using MATLAB matrix operations.

RK-4 Processing:

```python
  function z= rk(z,N,timestep,G,M) 
    F1 = State_Space(z,N,G,M); 
    deltax1 = timestep/2*F1; 
    F2 = State_Space(z+deltax1,N,G,M); 
    deltax2 = timestep/2*F2;
    F3 = State_Space(z+ deltax2,N,G,M);
    deltax3 = timestep*F3;
    F4 = State_Space(z+ deltax3,N,G,M);

    z  = z + (1/6)*timestep*(F1+2*F2+2*F3+F4);
  
```

Various input conditions, such as specifying particle locations within a galaxy spiral shape offset to planes using 3D polar rotations can be performed using

``` math
r(\phi) = \frac{A}{log(B tan\frac{\phi}{2N})}
```

Output simulations:

<p align="center">
    <img src="https://github.com/matthewsimpsonaero/Basketball-numerical-method/blob/main/Final.gif" alt="Alt Text" width="400" height="300"/>
</p>
