function FastTiffSave(vid,filename,frame_start,frame_end,imagedesc)
%function FastTiffSave(vid,filename,frame_start,frame_end,image_description)
%Saves Image Stack at speed comparable to ImageJ
%
%Typical implementation: FastTiffSave(vid,filename);
%
%Inputs:
%vid = image stack (n x m x t matrix)
%filename = destination to save stack
%frame_start (optional) = first frame to save (default: 1)
%frame_end (optional) = last frame to save (default: t)
%image_description (optional) = ImageDescription tiff tag as uint8/char
%
%Outputs: None; Video will be saved as .tif
%
%Uses Fast_Tiff_Write by R.Harkes 05-07-2019
%Fast_BigTiff_Write was written as adaption of code by Harkes for saving
%large tiffs
%
%Jan 20, 2020
%J.M.Stujenske
%
%Edit: March 17, 2021: added ability to add an image description, for
%better compatibility with multi-color stacks
%
%Here is code to add an imagej header that will indicate a multi-color
%z-stack given n_z (number of z-stacks). Modify other parameters to your
%preference.
%     imagedesc=[uint8('ImageJ=1.53c'),uint8(10),...
%     uint8(['images=',num2str(n_image_tot)]),uint8(10),...
%     uint8(['channels=',num2str(n_ch)]),uint8(10),... 
%     uint8(['slices=',num2str(n_z)]),uint8(10),... 
%     uint8(['hyperstack=true']),uint8(10),... 
%     uint8(['mode=grayscale']),uint8(10),... 
%     uint8(['unit=micron']),uint8(10),... 
%     uint8(['spacing=',num2str(pixel_size)]),uint8(10),... 
%     uint8(['loop=false']),uint8(10),... 
%     uint8(['min=0.0']),uint8(10),... 
%     uint8(['max=',num2str(ceil(max(vid(:)))),'.0']),uint8(10)];
%
filetype=class(vid);
switch filetype
                    case {'double'}
                        error('64 bit precision not yet supported')
                    case {'single'}
                        bps = 4;sf=3;
                    case {'uint16'}
                        bps = 2;sf=1;
                    case {'uint8'}
                        bps = 1;sf=1;
                    otherwise
                        error('class not supported')
end
im_size=size(vid);
if nargin<5 || isempty(imagedesc)
    imagedesc=[];
else
    imagedesc=uint8(imagedesc);
end
if nargin<4 || isempty(frame_end)
    if length(im_size)>2
    frame_end=im_size(3);
    else
        frame_end=1;
    end
end
if nargin<3 || isempty(frame_start)
    frame_start=1;
end
nFrames=frame_end-frame_start+1;
predicted_filesize=(im_size(1)*im_size(2)*bps+15*12)*nFrames+8;
if predicted_filesize<3.99e9 %make sure less than 3.99GB, to be conservative
    disp('Saving as Regular Tiff.')
    TiffWriter=Fast_Tiff_Write(filename,[],[],imagedesc);
else
    disp('Big Video Detected. Saving as BigTiff.')
    TiffWriter=Fast_BigTiff_Write(filename,[],[],imagedesc);
end
tic;for a=frame_start:frame_end;TiffWriter.WriteIMG(vid(:,:,a)');end;close(TiffWriter);duration=toc;
disp(['File Saved in ',num2str(duration),' seconds'])