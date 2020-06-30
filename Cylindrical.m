  classdef Cylindrical < ChiDist
    properties
        R = 5; % radius of sphere in [mm]
        theta = 0; % angle between main axis of cylinder and z-axis (in radians)
    end
    
    methods
       function obj = Cylindrical( matrix, image_res, R, theta, sus )
          % Method Spherical (constructor)
          
          % - Call superclass constructor.
          obj = obj@ChiDist( matrix, image_res, sus, 'Cylindrical' ) ;

          % - Define class properties.
          obj.R = R;
          obj.theta = theta;
          
          obj.cylindrical();
       end
       
       function obj = cylindrical(obj)
          % define image grid
          [x,y,z] = ndgrid(linspace(-(obj.matrix(1)-1)/2,(obj.matrix(1)-1)/2,obj.matrix(1)),linspace(-(obj.matrix(2)-1)/2,(obj.matrix(2)-1)/2,obj.matrix(2)),linspace(-(obj.matrix(3)-1)/2,(obj.matrix(3)-1)/2,obj.matrix(3)));
          
          % in-plane radial position (in [mm])
          r = sqrt((x.*obj.image_res(1)).^2 + (y.*obj.image_res(2)).^2);
          
          obj.volume = zeros(obj.matrix(1), obj.matrix(2), obj.matrix(3));
          
          obj.volume(r <= obj.R ) = obj.sus(1);
          obj.volume(r > obj.R ) = obj.sus(2);
          
          % rotate chi distribution about the y-axis
          t = [cos(obj.theta)  0      -sin(obj.theta)   0
              0             1              0     0
              sin(obj.theta)    0       cos(obj.theta)   0
              0             0              0     1];
          tform = affine3d(t);
          obj.volume = imwarp(obj.volume,tform);
            
       end
       
    end

 end