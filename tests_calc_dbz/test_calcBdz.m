%% This script uses the phantoms to test calc_dbz.
% The results are compared to the analytical solution, see Chapter 25 of
% Magnetic Resonance Imaging : Physical Principles and Seauence Design,
% Second Edition, Brown R., Cheng Y.,Haacke E., Thompson M., Venkatsen R.

phantom = "sphere"
b0 = 1; %[T]
sectionz = round(dim(3) / 2) + 1;
sectiony = round(dim(2) / 2) + 1;
sectionx = round(dim(1) / 2) + 1;

switch(phantom)
%%  An anisotropic rectangular susceptibility in a "little" volume
    case "rect" 
        dim = [16, 16, 16];
        res = [1, 1, 1];
        susin = 1;
        susout = 0;
        sus = zeros(dim) + susout;
        sus(7:10, 8:9, 6:11) = susin;
        if param == 0
        else
            % Translation toward down-right to make the matrix symmetric 
            sus = padarray(sus, [1 1 1], susout, 'pre');
            sus = sus(1:dim(1), 1:dim(2), 1:dim(3));
        end
        
%% A sphere
    case "sphere"
        dim = [256, 256, 256];
        res = [1, 1, 1]; % volume unit
        susin = -0.72e-6; 
        susout = -0.36e-6; 
        radius = 12 % volume unit
        spherical_sus = Spherical(dim , res, radius, [susin susout]);
        sus = spherical_sus.volume;
%% A cylinder
    case "cylinder"
        dim = [256, 256, 256];
        res = [1, 1, 1]; % volume unit
        susin = -0.72e-6; 
        susout = -0.36e-6; 
        radius = 12 % volume unit
        theta = 0 % rad, tilt between B0 and y
        phi = pi/2 %angle between x and measure axis in the xy plane (pi/2 for measure along y, 0 for measure along z)
        cylindrical_sus = Cylindrical(dim, res, radius, theta, [susin susout]);
        sus = cylindrical_sus.volume;

%% A sphere with a bigger volume (add a buffer)
    case "sphere_buffer"
        dim_sans_buffer = [256, 256, 256];
        dim = 2 * dim_sans_buffer;
        res = [1, 1, 1]; % volume unit
        susin = -0.72e-6; 
        susout = -0.36e-6; 
        radius = 86 % volume unit
        spherical_sus = Spherical(dim_sans_buffer , res, radius, [susin susout]);
        sus = spherical_sus.volume;
        sus = padarray(sus, dim_sans_buffer / 2, susout, 'post');
        sus = padarray(sus, dim_sans_buffer / 2, susout, 'pre');
    case "zubal"
        dim = [256, 256, 128];
        res = [1, 1, 1]; % volume unit

        spherical_sus = Spherical(dim , res, radius, [susin susout]);
        sus = spherical_sus.volume;
        
end