  classdef SheppLogan < ChiDist
    properties
    end
    
    methods
       function obj = SheppLogan( matrix, sus )
          % Method SheppLogan (constructor)
          
          % - Call superclass constructor.
          obj = obj@ChiDist( [matrix matrix matrix] , [1 1 1], sus, 'Shepp-Logan' ) ;
          
          obj.shepp_logan(obj.matrix);
       end
       
       function obj = shepp_logan(obj, matrix)       
          % Create a 3D Shepp_Logan volume
          obj.volume = phantom3d(matrix(1)); 

          obj.volume = obj.customize_shepp_logan(obj.volume, obj.sus(1), obj.sus(2), obj.sus(3), obj.sus(4));
       end
       
    end
    
    methods (Access = protected)
        function customVolume = customize_shepp_logan(obj, volume, class1, class2, class3, class4)
            customVolume = volume;
            
            % Set regions to T2
            customVolume(abs(volume-0.2)<0.001) = class1;
            customVolume(abs(volume-0.3)<0.001) = class2;
            customVolume(abs(volume-1)<0.001) = class3;
            customVolume((abs(volume)<0.0001)&volume~=0) = class4;
        end
    end
 end