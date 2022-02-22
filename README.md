Parameters like # of steps, velocity, and viscosity are defined directly in the source file (for now). Some param choices may result in instability.

Barrier position is determined by `if` statement in `c-reference/simulation.c:stream()`.

```
    # edit c-reference/main.c to have the # of sim steps you want
    cd c-reference
    make
    ./fluidsim simulation.fluidresults
    ../tools/render_fluid.py simulation.fluidresults simulation.avi
```
