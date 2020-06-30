  classdef Spherical < ChiDist
    properties
        R = 5; % default radius of sphere (in [mm])
    end
    
    methods
       function obj = Spherical( matrix, image_res, R, sus )
          % Method Spherical (constructor)
          
          % - Call superclass constructor.
          obj = obj@ChiDist( matrix, image_res, sus, 'Spherical' ) ;

          % - Define class properties.
          obj.R = R;
          
          obj.spherical();
       end
       
       function obj = spherical(obj)
          
          % define image grid
          [x,y,z] = ndgrid(linspace(-(obj.matrix(1)-1)/2,(obj.matrix(1)-1)/2,obj.matrix(1)),linspace(-(obj.matrix(2)-1)/2,(obj.matrix(2)-1)/2,obj.matrix(2)),linspace(-(obj.matrix(3)-1)/2,(obj.matrix(3)-1)/2,obj.matrix(3)));
          
          % radial position (in [mm])
          r = sqrt((x.*obj.image_res(1)).^2 + (y.*obj.image_res(2)).^2 + (z.*obj.image_res(3)).^2);
          
          obj.volume = zeros(obj.matrix(1), obj.matrix(2), obj.matrix(3));
          
          obj.volume(r <= obj.R ) = obj.sus(1);
          obj.volume(r > obj.R ) = obj.sus(2);
            
       end
       
    end

 end