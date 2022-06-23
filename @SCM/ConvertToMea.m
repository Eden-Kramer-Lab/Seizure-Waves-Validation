function ConvertToMea(scm)
    if ~scm.save, return, end
	if ~ismember('label', fieldnames(scm)) || isempty(scm.label), scm.label = 'SCM'; end
    
    
	files = dir(sprintf('%s_%d_*mat', scm.basename, scm.sim_num));
	addpath(files(1).folder);
	load(files(1).name, 'last');
	cmap = 1-gray;
	im = round(rescale(last.Ve) * (length(cmap) - 1)) + 1;
	movQ(numel(files) - 1) = im2frame(im, cmap);
	movV(numel(files) - 1) = im2frame(im, cmap);
    movK(numel(files) - 1) = im2frame(im, cmap);
    
	[qe, ve, k, qi, tt, dii, vi] = deal(cell(numel(files) - 1 , 1));

    file_inds = cellfun(@(f) strsplit(f, {'_', '.'}), {files.name}, 'uni', 0);
    file_inds = cellfun(@(f) str2double(strrep(f{end - 1}, 'M', '-')), file_inds);
    [~, file_order] = sort(file_inds, 'ascend');
    
	ii = 1;
    for f = files(file_order)'
		if strfind(f.name, 'info'), continue, end
		load(f.name, 'last', 'NP', 'time');
		disp(f.name)
		disp(ii)
		im = round(rescale(last.Qe) * (length(cmap) - 1)) + 1;
		movQ(ii) = im2frame(im, cmap);
		movV(ii) = im2frame(round(rescale(last.Ve) * (length(cmap) - 1)) + 1, cmap);
        movK(ii) = im2frame(round(rescale(last.K, 0, 1, 'InputMin', 0, 'InputMax', 1.5) * (length(cmap) - 1)) + 1, cmap);
        tt{ii} = time;
		qe{ii} = NP.Qe;
		ve{ii} = NP.Ve;
        ii = ii + 1;
	end
	
	
    qe_mat = cat(1, qe{:});
	time = cat(1, tt{:});
	sample_rate = min(round(1/mean(diff(time))/1e3)*1e3, scm.subsample);
	dt = 1 / sample_rate;
	nt = size(qe_mat, 1);
	inds = interp1(time, 1:nt, time(1):dt:time(end), 'nearest');
	time =@() time(1):dt:time(end);
	qe_mat = qe_mat(inds, :, :);
	
	
	mea = create_mea( ...
		-qe_mat, ...  % 'firing_rate', single(qe_mat), ... (This just takes up space)
		'SamplingRate', sample_rate, ... 
		'Padding', scm.padding, ...
        'Patient', scm.label, ...
        'Seizure', num2str(scm.sim_num), ...
		'Time', time ... 
		);
    mea.GridSize = scm.dimsNP;
	

	fprintf('Saving %s ... ', mea.Path);
	save(scm.mea_path, '-struct', 'mea', '-v7.3');
	m = matfile(sprintf('%s_%d_info', scm.basename, scm.sim_num), 'Writable', true);
	m.Qe_movie = movQ;
	m.Ve_movie = movV;
    m.K_movie = movK;
	
	fprintf('Done.\n')
	
	fprintf('Done.\n')
	
end


