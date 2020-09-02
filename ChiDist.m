classdef ChiDist < handle
    properties
        % volume properties
        type
        matrix
        image_res
        volume
        sus
    end
    
    methods
        function obj = ChiDist( matrix, image_res, sus, type )
            % Method chi_dist (constructor)
            obj.matrix  = matrix ;
            obj.image_res  = image_res ;
            obj.sus  = sus ;
            
            if nargin < 2
                obj.type = 'generic' ;
            else
                obj.type = type ;
            end
        end
        
        function vol = save(obj, fileName, saveFormat)
            % Get magnitude data
            % fileName: String. Prefix of filename (without extension)
            % saveFormat: 'nifti' (default) or 'mat'
            
            vol = obj.volume;
            
            if ~exist('saveFormat', 'var')
                warning('No save format given - saving to NIfTI')
                saveFormat = 'nifti';
            end
            
            
            if strcmp(fileName(end-3:end), '.nii')
                if ~strcmp(saveFormat, 'nifti')
                    warning('File extension and saveFormat do not match - saving to NIfTI format')
                    saveFormat = 'nifti';
                end
                fileName = fileName(1:end-4);
            elseif strcmp(fileName(end-3:end), '.mat')
                if ~strcmp(saveFormat, 'mat')
                    warning('File extension and saveFormat do not match - saving to MAT format')
                    saveFormat = 'mat';
                end
                fileName = fileName(1:end-4);
            end
            
            switch saveFormat
                case 'nifti'
                    if strcmp(obj.type,'Zubal')
                        nii_vol = make_nii(vol);
                    else
                        nii_vol = make_nii(imrotate(fliplr(vol), -90));
                        
                        % set image resolution in nifti header
                        nii_vol.hdr.dime.pixdim(2) = obj.image_res(1);
                        nii_vol.hdr.dime.pixdim(3) = obj.image_res(2);
                        nii_vol.hdr.dime.pixdim(4) = obj.image_res(3);
                    end
                    
                    save_nii(nii_vol, [fileName '.nii']);
                    
                case 'mat'
                    save([fileName '.mat'], 'vol')
                    
            end
        end
    end
    
end

