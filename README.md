# Seizure-Waves-Validation

## Notes


## Usage

### Generate simulations
Use the `SCM` class to make a seizing cortical model simulation. There are lots of options in here, but what's shown below will generate the simulations with the same parameters from Figure 1 of [REF](doi). 

The property `base_dir` controls where the simulations are saved. The default behavior is to save everything in a folder called 'SCM' in the current working directory (`~/SCM`). The full state of the simulation is saved at 1 s intervals in files with names like 
  
    ~/SCM/LABEL/LABEL_NN_TTT.mat
    
where _NN_ represents the simulation number (`SCM.sim_num`) and _TTT_ is the starting time of the segment. The simulated MEA data (the small subregion of the simulation that mimics an MEA recording) is saved to 

    ~/SCM/LABEL_SeizureNN_MEA_PP_PP.mat
    
where _PP_ indicates the duration of the pre- and post-ictal periods (`SCM.padding`). In both naming schemes, _LABEL_ is any string value (`SCM.label`).

Figure 1A: A fixed source simulation (no ictal wavefront)
```matlab
scm = SCM('FS');

% Use the following to see where the outputs will be saved:
disp(scm.basename);
disp(scm.mea_path);

% Start the simulation
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
Assumes that the MEA data is stored in a .mat file with the following fields:

    BadChannels                 18x1                        144  double              
    Data                   1875001x96                 360000192  int16               
    Duration                     1x1                          8  double              
    Map                         10x10                       800  double              
    Padding                      1x2                         16  double              
    Patient                      1x2                          4  char                
    Position                    96x2                       1536  double              
    SamplingRate                 1x1                          8  double              
    Seizure                      1x1                          2  char                
    Time                         1x1875001              7500004  single              
  
  
```matlab
mea = MEA(<path_to_mea_file>);
```

### Identify ictal wavefronts
```matlab
iw = mea.IW();
iw.plot()
iw.plot2D()
```

### Estimate traveling wave directions
The following will save the wave direction estimates to a folder called
_WaveFits_ in your current directory (`~/WaveFits`). The D-method can take a long time to 
run, so we have provided a helper script ([d10_helper.m](helpers/d10_helper.m))
that runs the analysis in segments. This allows you to pick up where you left off
if your run is interrupted. To do so, simply keep re-running `d10_helper;` (with the same `mea`
data) until the analysis is complete.

```matlab
mea = MEA(<path_to_mea_file>);
s = mkdir(['WaveFits/' mea.Name]);

% M-method direction estimatees
times = mea.get_discharge_times('lfp_cross');
M = mea.max_descent(times, 'halfwin', 0.05);
save(fullfile('WaveFits', mea.Name, 'M'), 'M');

% D-method direction estimates
% This method is slow so we recommend using the helper script (shown below)
step = 0.1;
times = mea.Time(1):step:mea.Time(end);
D = mea.delays(times, ...
        'halfwin', 5, 'fband', [1 13], 'minfreq', 3);
save(fullfile('WaveFits', mea.Name, 'D'), 'D');

% D-method alternative
d10_helper;
  % if this is interrupted, reinstantiate `mea` and run `d10_helper` again.

```


### Identify intervals with stable traveling wave directions

```matlab
% Load pre-saved discharge arrival times
M = load('WaveFits/P1_Seizure1/M.mat').M;
intervals = M.stable_intervals();
scatter(intervals.time, intervals.dir_trend, [], intervals.phase_num);
```

