classdef FBFest < handle
    properties
        % volume properties
        matrix
        image_res
        volume
        sus % USI (not ppm)
        type
        sus_ext % USI (not ppm)
        b0 % [T]
    end
    
    methods
        function obj = FBFest( sus, image_res, matrix, sus_ext, b0, varargin )
            % Method FBFest (constructor)
            obj.matrix  = matrix ;
            obj.image_res  = image_res ;
            obj.sus  = sus ;
            obj.b0 = b0 ;
            obj.sus_ext = sus_ext;
            
            if nargin > 3
                obj.type = varargin;
            end
            
            obj.volume = zeros(obj.matrix(1), obj.matrix(2), obj.matrix(3));
            obj.calc_dBz();
        end
        
        function obj = calc_dBz(obj)
            %%---------------------------------------------------------------------- %%
            %% Define constants
            %%---------------------------------------------------------------------- %%
            
            % k-space window
            k_max = 1./(2.*obj.image_res);
            interval = 2 * k_max ./ obj.matrix;

            % define k-space grid
            [kx,ky,kz] = ndgrid(-k_max(1):interval(1):k_max(1) - interval(1),-k_max(2):interval(2):k_max(2) - interval(2),-k_max(3):interval(3):k_max(3) - interval(3));
            
            %%---------------------------------------------------------------------- %%
            %% Compute Bdz
            %%---------------------------------------------------------------------- %%
            
            % compute the fourier transform of the susceptibility distribution
            FT_chi = fftshift(fftn(fftshift(obj.sus)));
            
            % calculate the scaling coefficient 'kz^2/k^2'
            k2 = kx.^2 + ky.^2 + kz.^2;
            k_scaling_coeff = kz.^2./k2;
            
            % compute Bdz (the z-component of the magnetic field due to a
            % sphere, relative to B0) FFT. The region of interest is
            % assumed to be surrounded by a region of infinite extent whose
            % susceptibility is equal to the susceptibility on the origin
            % corner of the matrix.
            bdzFFT = obj.b0 * (1/3 - k_scaling_coeff).*FT_chi;
            bdzFFT(k2 == 0) = obj.b0 * obj.sus_ext / 3 * prod(obj.matrix);
            obj.volume = ifftshift(ifftn(ifftshift(bdzFFT)));
            
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
                        nii_vol = make_nii(real(1e6*vol)); % save real part in ppm
                    else
                        nii_vol = make_nii(imrotate(fliplr(real(1e6*vol)), -90)); % save real part in ppm
                    end
                       
                    %set image resolution in nifti header
                    nii_vol.hdr.dime.pixdim(2) = obj.image_res(1);
                    nii_vol.hdr.dime.pixdim(3) = obj.image_res(2);
                    nii_vol.hdr.dime.pixdim(4) = obj.image_res(3);

                    
                    save_nii(nii_vol, [fileName '.nii']);
                    
                case 'mat'
                    save([fileName '.mat'], 'vol')
                    
            end
        end
        
        function vol = plot(obj)
            %%---------------------------------------------------------------------- %%
            %% Plot some results
            %%---------------------------------------------------------------------- %%
            
            % define image grid for plotting
            [x,y,z] = ndgrid(linspace(-(obj.matrix(1)-1)/2,(obj.matrix(1)-1)/2,obj.matrix(1)),linspace(-(obj.matrix(2)-1)/2,(obj.matrix(2)-1)/2,obj.matrix(2)),linspace(-(obj.matrix(3)-1)/2,(obj.matrix(3)-1)/2,obj.matrix(3)));
            
            % plot saggital view of Bdz
            figure;
            imagesc(squeeze(1e6*real(obj.volume(round(obj.matrix(1)/2),:,:))));
            colorbar;
            title('Bdz [ppm]')
            xlabel('z-position')
            ylabel('y-position')
            
            % plot Bdz along the z-axis
            figure;
            plot(squeeze(z(round(obj.matrix(1)/2),round(obj.matrix(2)/2),:))*obj.image_res(3),1e6*squeeze(real(obj.volume(round(obj.matrix(1)/2),round(obj.matrix(2)/2),:))),'-.k','Linewidth',2);
            xlabel('z-position [mm]')
            ylabel('B_{dz} [ppm]')
            
            % plot Bdz along the x-axis
            figure;
            plot(squeeze(x(:,round(obj.matrix(2)/2),round(obj.matrix(3)/2)))*obj.image_res(1),1e6*squeeze(real(obj.volume(:,round(obj.matrix(2)/2),round(obj.matrix(3)/2)))),'-.k','Linewidth',2);
            xlabel('x-position [mm]')
            ylabel('B_{dz} [ppm]')

        end


    end
    
end

