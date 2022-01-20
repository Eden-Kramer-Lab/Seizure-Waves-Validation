# Seizure-Waves-Validation

## Usage

### Generate simulations
Use the `SCM` class to make a seizing cortical model simulation. There are lots of options in here, but what's shown below will generate the simulations with the same parameters from Figure 1 of [REF](doi). 

The property `base_dir` controls where the simulations are saved. The default behavior is to save everything in a folder called 'SCM' in your current working directory. The full state of the simulation is saved at `t_step` s intervals in files with names like `<basename>_<sim_num>_<t0>.mat`. The simulated MEA data (the small subregion of the simulation that mimics an MEA recording) is saved to `<base_dir>\<label>_Seizure<sim_num>_MEA_<padding[1]>_<padding[2]>.mat`. 

Figure 1A: A fixed source simulation (no ictal wavefront)
```matlab
scm = SCM('FS');
scm.Run();
```

Figure 1E: A simulation where the fixed source location shifts midway through the seizure
```matlab
% Create a fixed source simulation. Results will be saved to
% folders 'SCM' (full state of sim) and 'FS' (MEA)
scm = SCM('SW');
scm.Run();
```

<details><summary><b>Some other methods</b></summary>
  
  * `show_layout()`: Show the simulation layout (locations of the irritative zone, fixed sources, MEA, and IW source)
  * `rotate(theta)`: Rotate the locations of fixed and IW sources by angle `theta` about the center of the simulation
  * `Preview(show_layout=false)`: Generates a set of figures with panels showing excitatory firing rates (`Qe`) and potassium (`K`) states at each save point (defaults to every 1 second; controlled by parameter `t_step`)
  
</details>

### Load MEA data
```matlab
mea = MEA(<path_to_mea_file>);
```

### Identify ictal wavefronts
```matlab
iw = mea.IW();
iw.plot9)
```

### Estimate traveling wave directions
```matlab
```


### Identify intervals with stable traveling wave directions



