function P_rotated = rotate_phantom(P, theta, varargin)
    % 20/05/2025
    % Sofia Pearson
    % Rotate phantom
    
    % inputs:
    % image to rotate - P 
    % angle of rotation - theta
    % varargin: string for interpolation method - 'nearest', 'bilinear', 'bicubic'

    % Define default and allowed options
    validMethods = {'nearest', 'bilinear', 'bicubic'};
    
    if ~isempty(varargin)
        inputStr = lower(string(varargin{1}));  % Convert to lowercase

        if any(strcmp(inputStr, validMethods))
            interp_method = inputStr;
        else
            error('Invalid interpolation method. Choose one of: %s', strjoin(validMethods, ', '));
        end
    else 
            interp_method = 'bilinear';
    end

    
    % to preserve shape, pad image before rotating
    npx = size(P, 1);
    pad_size = ceil(sqrt(2) * npx);  % to accommodate rotation
    padded_img = padarray(P, [pad_size pad_size], 0, 'both');
    
    % Rotate

    % TEST: Remove interpolation:
    % rotated_padded = imrotate(padded_img, theta, interp_method, 'crop'); 
    rotated_padded = imrotate(padded_img, theta, 'crop');
    pad_size = size(rotated_padded,1);
    
    % Crop back to npx size
    start_idx = floor((pad_size - npx) / 2) + 1;
    end_idx = start_idx + npx - 1;
    
    
    % Crop the matrix
    P_rotated = rotated_padded(start_idx:end_idx, start_idx:end_idx);

    %plotting:
    figure,
    tiledlayout(2,2);
    
    nexttile
    imagesc(P), title('phantom'), axis square, colormap gray;
    
    nexttile
    imagesc(padded_img),  title('padded phantom'), axis square, colormap gray;
    
    nexttile
    imagesc(rotated_padded), title('rotated, padded phantom'), axis square, colormap gray;
    
    nexttile
    imagesc(P_rotated),  title('final rotated phantom'), axis square, colormap gray;

end



  



