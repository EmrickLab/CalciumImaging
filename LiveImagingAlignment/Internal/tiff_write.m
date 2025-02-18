function tiff_write(output_filename, img, spacing, output_compression)

if ~exist('spacing', 'var')
    spacing = [1,1,1];
end

if ~exist('output_compression', 'var')
    output_compression = 'None';
end

% preparing the 'description' key of the metadata, this mainly adds the Z spacing of the image:
Description = strjoin({ ...
    'ImageJ=1.51p', ...
    ['images=', num2str(size(img,3))], ...
    ['slices=', num2str(size(img,3))], ...
    'unit=\u00B5m', ...
    ['spacing=', num2str(spacing(3))], ...
    'loop=false', ...
    ''}, '\n');

% iteratively saving the slices of the image:
save_ok = false;
while ~save_ok
    imwrite(img(:,:,1), output_filename, 'Resolution', 1/(spacing(1)), 'Description', Description, 'Compression', output_compression);
    for z = 2 : size(img,3)
        imwrite(img(:,:,z), output_filename, 'WriteMode', 'append', 'Resolution', 1/(spacing(1)), 'Description', Description, 'Compression', output_compression);
    end
    save_ok = true;
end

