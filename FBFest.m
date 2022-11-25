classdef FBFest < handle
    properties
        % volume properties
        phantom
        matrix
        dim_with_buffer
        image_res
        volume % [T] The field variation is in unit of B0, so needs to be multiplied by the strength field B0 if necessary
        sus % SI unit (not ppm)
        type
        sus_ext % SI unit (not ppm)
    end
    
    methods
        function obj = FBFest(type, sus, image_res, matrix, sus_ext, varargin )
            % Method FBFest (constructor)
            % type      : the type of the phantom : 'spherical',
            % 'cylindrical', 'Zubal', 'Shepp-Logan' or '' for an other phantom
            % sus       : the 3D distribution of susceptibility
            % image_res : a vector of the resolution in the 3 directions
            % matrix    : a vector of the dimensions of the field_view (matching sus dimensions)
            % sus_ext   : a scalar, the susceptibility of the external medium
            %              (the region of interest is assumed to be surrounded
            %              by a region of infinite extent whose susceptibility is sus_ext)
            % varargin  : a dimension with the buffer can be chosen
            obj.type = type;
            obj.matrix  = matrix ;
            obj.image_res  = image_res ;
            obj.sus  = sus ;
            obj.sus_ext = sus_ext;            
            
            if nargin > 5
                obj.dim_with_buffer = varargin{1};
            else
                obj.calc_buffer();
            end
            
            obj.volume = zeros(obj.matrix(1), obj.matrix(2), obj.matrix(3));
            obj.calc_dBz();
            %obj.plot()
        end
        
        function obj = calc_dBz(obj)
            %%---------------------------------------------------------------------- %%
            %% Define constants
            %%---------------------------------------------------------------------- %%
            % k-space window
            k_max = 1./(2.*obj.image_res);
            interval = 2 * k_max ./ obj.dim_with_buffer;

            % define k-space grid
            [kx,ky,kz] = ndgrid(-k_max(1):interval(1):k_max(1) - interval(1),-k_max(2):interval(2):k_max(2) - interval(2),-k_max(3):interval(3):k_max(3) - interval(3));            

            %%---------------------------------------------------------------------- %%
            %% Compute Bdz
            %%---------------------------------------------------------------------- %%            

            % calculate the kernel
            % make sure that the k-space window is the same that the
            % one of the fft of susceptibility
            k2 = kx.^2 + ky.^2 + kz.^2;
            kernel = fftshift(1/3 - kz.^2./k2); % For B0 = 1T
            kernel(1, 1, 1) = 1/3; % for B0 = 1T

            % compute the fourier transform of the susceptibility
            % distribution using the linearity of the FT
            FT_chi = fftn(obj.sus, obj.dim_with_buffer);
            FT_chi(1, 1, 1) = FT_chi(1,1,1) + prod(obj.dim_with_buffer) * obj.sus_ext;

            %FT_chi = fftn(obj.sus - obj.sus_ext, obj.dim_with_buffer); % region of interest
            %FT_chi(1, 1, 1) = FT_chi(1,1,1) + prod(obj.dim_with_buffer) * obj.sus_ext; % external susceptibility at k=0

            % compute Bdz (the z-component of the magnetic field due to a
            % sphere, relative to B0) FFT. 
            bdzFFT = kernel .* FT_chi;
            obj.volume = real(ifftn(bdzFFT));

            %% Truncate the result to the initial dimension of the susceptibility
            obj.volume = obj.volume(1:obj.matrix(1), 1:obj.matrix(2), 1:obj.matrix(3));
            
        end
        
        function obj = calc_buffer(obj)
           switch(obj.type)
               case 'spherical'
                   disp('Calculation dist ROI to have a default buffer size...')
                   [dist_ROI, ] = calc_dist_ROI(obj.sus);
                   disp('ended.')
                   diam_approx = obj.matrix(1) - 2 * dist_ROI;
                   side = max([pow2(nextpow2(5 * diam_approx)), obj.matrix(1)]); 
                   obj.dim_with_buffer = [side, side, side]; % TODO
               case 'cylindrical'
                   obj.dim_with_buffer = 2 * obj.matrix; % TODO
               case 'Zubal'
                   obj.dim_with_buffer = [512, 512, 512]; % [256, 256, 384] voxels added
               case 'Shepp-Logan'
                   obj.dim_with_buffer = 2 * obj.matrix; % todo : optimize
               case ''
                   warning('No type and no buffer dimension given for FBFest estimation - empirical buffer used')
                   obj.dim_with_buffer = 2 * obj.matrix; % TODO
               otherwise
                   warning('No type and no buffer dimension given for FBFest estimation - empirical buffer used')
                   obj.dim_with_buffer = 2 * obj.matrix; % TODO
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
                        nii_vol = make_nii(real(1e6*vol)); % save real part in ppm
                    else
                        nii_vol = make_nii(imrotate(fliplr(real(1e6*vol)), -90)); % save in ppm with a rotation to have a vertical z axis
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
            grid on
            
            % plot Bdz along the x-axis
            figure;
            plot(squeeze(x(:,round(obj.matrix(2)/2),round(obj.matrix(3)/2)))*obj.image_res(1),1e6*squeeze(real(obj.volume(:,round(obj.matrix(2)/2),round(obj.matrix(3)/2)))),'-.k','Linewidth',2);
            xlabel('x-position [mm]')
            ylabel('B_{dz} [ppm]')
            grid on

        end


    end
    
end

