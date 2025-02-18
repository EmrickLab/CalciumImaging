function correct_movements

% setting user parameters:
filename = ['C:\Users\jjemrick\Desktop\2025-02-16 Calcium Imaging Analysis Pipeline\Sample Set\raw.tif']; % path to 2D+t input image (suffix: '.tif' ; type: 'uint8')
stacks = [1,200; 201,400];

% setting more advancd parameters:
image_subsampling = 4; % subsampling input image to speed up computations: I(image_subsampling:image_subsampling:end,image_subsampling:image_subsampling:end,:);
optimizer_max_iterations = 500; % max iterations in affine registration
optimizer_max_step_length_factor = 0.1; % a multiplicative coefficient to the maximal step size in the affine registration; otherwise registration diverges (within alignment)
optimizer_min_step_length_factor = 0.1; % a multiplicative coefficient to the minimal step size in the affine registration; otherwise registration diverges (within alignment)
relaxation_factor = 0.5; % relaxation factor of the registration optimizer (within alignment)
n_pyramid_levels = 3; % n pyramid levels for the affine registration
first_time_point_weight = 1; % how much relative weight to give the first time point image as opposed to the image of the previous time point
histeq_n_bins = 2^8/2; % n bins for histogram equalization - appropriate for uint8 data type
gauss_filt_sd_source = 1; % how much to smooth the source image
gauss_filt_sd_target = 1.5; % how much to smooth the target image
mask_filter_gauss_sd = 40; % sd of gaussian smoothing of the binary mask
perform_masked_alignment = 1; % whether to do masked alignment after the initial translation
optimizer_initial_radius_factor = 0.1; % multiplicative coefficient to the initial radius (between alignment)

% defining required paths:
mask_filename = fullfile(fileparts(filename), 'ROI_mask.mat');
within_stack_tforms_filename = fullfile(fileparts(filename), 'within_stack_tforms.mat');
between_stack_tforms_filename = fullfile(fileparts(filename), 'between_stack_tforms.mat');
final_tforms_filename = fullfile(fileparts(filename), 'final_tforms.mat');
aligned_representative_slices_realigned_filename = fullfile(fileparts(filename), 'stack_representative_slices_realigned.tif');
output_image_filename = regexprep(filename, '\.tif', '_aligned.tif');

% reading raw image and preparing its spatial reference:
fprintf('Loading image... ');
raw_I = bigread4(filename);
fprintf(['Done.', newline, 'Image size: ', num2str(size(raw_I)), newline]);
ref_raw_image = imref2d([size(raw_I,1), size(raw_I,2)]);

% setting up the parameters of the registration:
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = optimizer_max_iterations;
optimizer.MaximumStepLength = optimizer.MaximumStepLength * optimizer_max_step_length_factor;
optimizer.MinimumStepLength = optimizer.MinimumStepLength * optimizer_min_step_length_factor;
optimizer.RelaxationFactor = relaxation_factor;

% setting up the correction of the transformation to the subsampling:
scaling_pre = eye(3);
scaling_pre([1,5]) = 1 / image_subsampling;
scaling_post = eye(3);
scaling_post([1,5]) = image_subsampling;

% asking the user to mark the region of interest:
fid = figure('Name', 'Please mark the region in focus');
slice_for_ROI = imgaussfilt(raw_I(:,:,1), gauss_filt_sd_source);
imshow(slice_for_ROI, []);
fh = imfreehand; %#ok<IMFREEH>
mask = fh.createMask;
mask = imgaussfilt(single(mask), mask_filter_gauss_sd);
close(fid);
save(mask_filename, 'mask');

% allocating memory for the 'within' transformations:
within_tforms = cell(size(stacks,1),1);

I = cell(size(stacks,1),1);

% in each iteration we perform the 'within' alignment of a single stack:
for s = 1 : size(stacks,1)
    
    I{s} = raw_I(:,:,stacks(s,1):stacks(s,2));
    
    within_tforms{s} = cell(diff(stacks(s,:))+1,1);
    within_tforms{s}{1} = affine2d(eye(3));
    
    % in each iteration we perform align one frame to the first frame within the same stack:
    for i = 2 : size(I{s},3)
        
        fprintf(['Within alignment stack ', num2str(s), '/', num2str(size(stacks,1)), ' slice ', num2str(i), '/', num2str(size(I{s},3)), '... ']);
        tic;
        
        % preparing the source and the target images:
        source = uint8(I{s}(:,:,1) * first_time_point_weight + I{s}(:,:,i-1) * (1-first_time_point_weight));
        target = I{s}(:,:,i);
        
        % calling the function that aligns the target to the source:
        within_tforms{s}{i} = align_single_slices(source, target, image_subsampling, gauss_filt_sd_source, gauss_filt_sd_target, histeq_n_bins, mask, ...
            optimizer, metric, n_pyramid_levels, scaling_pre, scaling_post, perform_masked_alignment);
        
        % applying the alignment to the raw slice so we could save and examine it later:
        I{s}(:,:,i) = imwarp(I{s}(:,:,i), ref_raw_image, within_tforms{s}{i}, 'cubic', 'OutputView', ref_raw_image);
        
        toc;
        
    end
    
    fprintf(newline);
    
end

% saving the all the 'within' transformation matrices and the mask file we generated:
save(within_stack_tforms_filename, 'within_tforms');
save(mask_filename, 'mask');

% generating the representative frame from each stack, so now we could align the stacks to one another:
representative_stack_slices = cellfun(@(x) uint8(prctile(x,10,3)), I, 'UniformOutput', false);
representative_stack_slices = cell2mat(permute(representative_stack_slices, [3,2,1]));

% updating the parameters for the alignment 'between' the stacks:
[optimizer, metric] = imregconfig('multimodal');
optimizer.MaximumIterations = optimizer_max_iterations;
optimizer.InitialRadius = optimizer.InitialRadius * optimizer_initial_radius_factor;

% allocating memory for the 'between' transfomration matrices:
between_tforms = cell(size(stacks,1),1);
between_tforms{1} = affine2d(eye(3));

% in each iteration we calculate the alignment 'between' the stacks, based on the representative frames:
realigned_representative_stack_slices = representative_stack_slices;
for s = 2 : size(representative_stack_slices,3)
    
    fprintf(['Between alignment stack ', num2str(s), '/', num2str(size(stacks,1)), '... ']);
    tic;
    
    between_tforms{s} = align_single_slices(representative_stack_slices(:,:,1), representative_stack_slices(:,:,s), image_subsampling, gauss_filt_sd_source, ...
        gauss_filt_sd_target, histeq_n_bins, mask, optimizer, metric, n_pyramid_levels, scaling_pre, scaling_post, perform_masked_alignment);
    realigned_representative_stack_slices(:,:,s) = imwarp(representative_stack_slices(:,:,s), ref_raw_image, between_tforms{s}, 'cubic', 'OutputView', ref_raw_image);
    
    toc;
    
end
fprintf(newline);

% saving the tranformations between the representative slides and the aligned representative slides:
save(between_stack_tforms_filename, 'between_tforms');
FastTiffSave(representative_stack_slices, aligned_representative_slices_realigned_filename);

% calculating the final image by integrating the 'within' and 'between' stack alignments:
final_tforms = cell(length(within_tforms),1);
frame_counter = 0;
aligned_image = raw_I;
for s = 1 : length(within_tforms)
    
    final_tforms{s} = cell(length(within_tforms{s}),1);
    
    for i = 1 : length(within_tforms{s})
        final_tforms{s}{i} = affine2d(between_tforms{s}.T * within_tforms{s}{i}.T);
        frame_counter = frame_counter + 1;
        aligned_image(:,:,frame_counter) = imwarp(aligned_image(:,:,frame_counter), ref_raw_image, final_tforms{s}{i}, 'cubic', 'OutputView', ref_raw_image);
    end
end

% saving the final transformation matrices and the final aligned image file:
fprintf('Saving final image... ');
tic;
save(final_tforms_filename, 'final_tforms');
FastTiffSave(aligned_image, output_image_filename);
disp(['Image is saved to: "', final_tforms_filename, '"']);
disp('done.');
toc;
fprintf(newline);


function tform = align_single_slices(source, target, image_subsampling, gauss_filt_source, gauss_filt_target, histeq_n_bins, mask, optimizer, metric, ...
    n_pyramid_levels, scaling_pre, scaling_post, perform_masked_alignment)

% smoothing and subsampling the source image:
source = imgaussfilt(source, gauss_filt_source);
source = source(image_subsampling : image_subsampling : end, image_subsampling : image_subsampling : end);

% smoothing and subsampling the target image:
target = imgaussfilt(target, gauss_filt_target);
target = target(image_subsampling : image_subsampling : end, image_subsampling : image_subsampling : end);

% equalizing the histogram of the target image to this of the current source:
target = imhistmatch(target, source, histeq_n_bins);

% subsampling the mask of the high contrast ROI:
mask = mask(image_subsampling : image_subsampling : end, image_subsampling : image_subsampling : end);

% setting up the references of the subsampledimages:
ref_subsampled = imref2d(size(source));

% we start alignment by performing initial approximation using translation without masking the source and target:
tform = imregtform(target, ref_subsampled, source, ref_subsampled, 'translation', optimizer, metric, 'PyramidLevels', n_pyramid_levels);

if perform_masked_alignment

    % then, we mask the background (i.e. low contrast) in both source and target:
    meanval = mean(source(mask >= 0.5));
    background = meanval * ones(size(mask));
    
    rounded_tform = affine2d(round(tform.T));
    target = imwarp(target, ref_subsampled, rounded_tform, 'nearest', 'OutputView', ref_subsampled);
    target = uint8(single(target) .* mask + background .* (1-mask));
    target = imwarp(target, ref_subsampled, invert(rounded_tform), 'nearest', 'OutputView', ref_subsampled, 'FillValues', meanval);
    target = int16(target) - int16(meanval);
    
    source = uint8(single(source) .* mask + background .* (1-mask));
    source = int16(source) - int16(meanval);

    % we continue by fine tuning of the alignment based on the now masked source and target:
    tform = imregtform(int16(target), ref_subsampled, int16(source), ref_subsampled, 'translation', optimizer, metric, 'PyramidLevels', n_pyramid_levels, 'InitialTransformation', tform);
    
end

% lastly, since we performed the alignments on subsampled images, we correct the transformation matrix so it would work on the original size image:
tform.T = scaling_pre * tform.T * scaling_post;


