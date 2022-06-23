function RunSimulation(scm)
% A wrapper for seizing_cortical_field. Runs the simulation in segments and
% saves each segment; loads the state from a previously saved run if found.

    % Get meta parameters
    BASENAME = scm.basename;
    DURATION = scm.duration;
    PADDING = scm.padding;
    SAVE = scm.save_output;
    SIM_NUM = scm.sim_num;
    T_STEP = scm.t_step;
    T0_START = scm.t0_start;

    % Create stimulus map
    if T0_START > -PADDING(1)
        % If not starting a fresh sim, use previously saved <last>
        try
            load(checkname(sprintf('%s_%d_%03d', BASENAME, SIM_NUM, T0_START - 1)), 'last')
        catch ME
            if ~(strcmpi('MATLAB:load:couldNotReadFile', ME.identifier))
                rethrow(ME);
            else
                last = scm.IC;
            end
        end
    else
        % ... otherwise, start a fresh sim
        last = scm.IC;  %Load the initial conditions to start.
    end

    K = PADDING(2) + DURATION;
    fig = [];
    for t0 = scm.t0_start:T_STEP:K-T_STEP  % For each step
        tic
        % Update time offset
        scm.t0 = t0;

        % ... show progress, 
        fprintf('Running %d / %d .. ', t0, floor(K) - 1);  

        % ... run simulation for duration T_STEP,
        [NP, time, last, fig] = ...  
            scm.seizing_cortical_field(min(T_STEP, K - t0), last, fig);

        % ... save the results of this run,
        if SAVE
            fprintf('Saving .. ')
            fname = checkname(sprintf('%s_%d_%03d', BASENAME, SIM_NUM, t0));
            save(fname, 'NP','time','last');
        end
        toc
        % ... update progress.
        fprintf('Done.\n')  
    end

end
