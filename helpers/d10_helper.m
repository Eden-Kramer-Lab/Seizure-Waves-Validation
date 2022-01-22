%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delays runs slowly so it gets interrupted fairly often. Keep
% running this until the full run is complete.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change these to adjust delay parameters
metric = 'D10';
delay_fun =@(times) mea.delays(times, ...
        'halfwin', 5, 'fband', [1 13], 'minfreq', 3);
disp(delay_fun)

step = 0.1;
times = mea.Time(1):step:mea.Time(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some constants you probably don't need to change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = 10;  % number of segments
temp_file = ['temp_' mea.Name metric '.mat'];
disp(temp_file)

% Compute times and how to segment them (where to put save points)
splits = round(linspace(1, numel(times)+1, K + 1));

% If temp_file exists, load it and pick up where you left off.
% Otherwise, run all segments.
try
    load(temp_file)
    which_splits = find(cellfun(@numel, {D.time}) == 0);
catch ME
    if contains(ME.message, 'No such file or directory.')
        which_splits = K:-1:1;
    else
        rethrow(ME)
    end
end
disp(which_splits)

% Compute the delays, saving each segment
for ii = which_splits
    tt = times(splits(ii):splits(ii + 1)-1);
    D(ii) = delay_fun(tt);
    save(temp_file, 'D');
    fprintf('savepoint %d\n', ii)
end

% Compile the delays into one object, save them to the "fits" file,
% and delete the temp_file
DD = WaveProp.resize_obj(D);
s = mkdir(['WaveFits/' mea.Name]);
save(['WaveFits/' mea.Name filesep metric], 'DD');
delete(temp_file);

