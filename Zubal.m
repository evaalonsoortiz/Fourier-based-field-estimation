  classdef Zubal < ChiDist
    properties
    end
    
    methods
       function obj = Zubal(zubal_fname)
          % Method Zubal (constructor)
          % must be used with a modified zubal phantom (distribution is not
          % possible, as per the conditions of the original zubal phantom): zubal_fname
          % the original zubal phantom can be found here: http://noodle.med.yale.edu/zubal/data.htm
          
          % absolute susceptibility [ppm] values are from Buch et al. MRM 73:2185?2194 (2015)
          % susceptibility of fat is from https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.27640
          % susceptibility of muscle is from https://www.researchgate.net/figure/Magnetic-susceptibility-of-different-matter_tbl1_288160127
          % susceptibility of venous and arterial blood from https://www.researchgate.net/figure/Magnetic-susceptibility-of-different-matter_tbl1_288160127
          % susceptibility of cartilage is from https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.26596
          sus = struct('sinus',0.2,'teeth',-12.3,'bone',-11.1,'thalamus',-9.02, ...
              'caudate_nucleus',-8.99,'putamen',-8.96,'globus_pallidus',-8.89,'GM',-9.03,'WM',-9.0830,'CSF',-9.05,...
              'air',0.35,'water',-9.05,'fat',-8.39,'muscle',-9.03,'venous_blood',-7.8,'arterial_blood',-9.3,...
              'cartilage',-9.055);
          
          % - Call superclass constructor.
          obj = obj@ChiDist( [256 256 128] ,[1.1 1.1 1.4], sus, 'Zubal') ;
         
          obj.zubal(zubal_fname, obj.sus);
       end
       
       function obj = zubal(obj, zubal_fname, sus)    %#ok<INUSD>
          % load zubal volume
          zubal_phantom = double(niftiread(zubal_fname));
          
          zubal_id = struct('outside_phantom',0,'skin',1,'cerebral_fluid',2,'spinal_cord',3,'skull',4,'spine',5,'dens_of_axis',70, ...
              'jaw_bone',71,'parotid_gland',72,'skeletal_muscle',9,'lacrimal_glands',74,'spinal_canal',75,'hard_palate',76, ...
              'cerebellum',77,'tongue',78,'pharynx',15,'esophagus',16,'horn_of_mandible',81,'nasal_septum',82,'white_matter',83, ...
              'superior_sagittal_sinus',84,'medulla_oblongata',85,'blood_pool',23,'frontal_lobes',89,'bone_marrow',26, ...
              'pons',91,'third_ventricle',92,'trachea',29,'cartilage',30,'occipital_lobes',95,'hippocampus',96,'pituitary_gland',97, ...
              'fat1',98,'fat2',22,'ear_bones',99,'turbinates',100,'caudate_nucleus',101,'zygoma',102,'insula_cortex',103,'sinuses_mouth_cavity',104, ...
              'putamen',105,'optic_nerve',106,'internal_capsule',107,'septum_pellucidium',108,'thalamus',109,'eyeball',110,'corpus_collosum',111, ...
              'special_region_frontal_lobes',112,'cerebral_falx',113,'temporal_lobes',114,'fourth_ventricle',115,'frontal_portion_eyes',116, ...
              'parietal_lobes',117,'amygdala',118,'eye',119,'globus_pallidus',120,'lens',121,'cerebral_aquaduct',122,'lateral_ventricles',123, ...
              'prefrontal_lobes',124,'teeth',125,'sigmoid_sinus',126);
          
          zubal_sus = struct('outside_phantom',obj.sus.air,'skin',obj.sus.water,'cerebral_fluid',obj.sus.CSF,'spinal_cord',(obj.sus.GM+obj.sus.WM)/2,...
              'skull',obj.sus.bone,'spine',obj.sus.bone,'dens_of_axis',obj.sus.bone,'jaw_bone',obj.sus.bone,'parotid_gland',obj.sus.water,...
              'skeletal_muscle',obj.sus.muscle,'lacrimal_glands',obj.sus.water,'spinal_canal',obj.sus.CSF,'hard_palate',obj.sus.bone,...
              'cerebellum',obj.sus.GM,'tongue',obj.sus.muscle,'pharynx',obj.sus.sinus,'esophagus',obj.sus.water,'horn_of_mandible',obj.sus.bone,...
              'nasal_septum',obj.sus.bone,'white_matter',obj.sus.WM,'superior_sagittal_sinus',obj.sus.venous_blood,...
              'medulla_oblongata',(obj.sus.GM+obj.sus.WM)/2,'fat1',obj.sus.fat,'fat2',obj.sus.fat,'blood_pool',(obj.sus.venous_blood+obj.sus.arterial_blood)/2,...
              'frontal_lobes',obj.sus.GM,'bone_marrow',obj.sus.CSF,'pons',obj.sus.WM,'third_ventricle',obj.sus.CSF,'trachea',obj.sus.sinus,...
              'cartilage',obj.sus.cartilage,'occipital_lobes',obj.sus.GM,'hippocampus',obj.sus.GM,'pituitary_gland',obj.sus.water,'ear_bones',obj.sus.bone,...
              'turbinates',(obj.sus.bone+obj.sus.water)/2,'caudate_nucleus',obj.sus.caudate_nucleus,'zygoma',obj.sus.bone,...
              'insula_cortex',obj.sus.GM,'sinuses_mouth_cavity',obj.sus.sinus,'putamen',obj.sus.putamen,'optic_nerve',obj.sus.WM,...
              'internal_capsule',obj.sus.WM,'septum_pellucidium',(obj.sus.WM+obj.sus.GM)/2,'thalamus',obj.sus.thalamus,'eyeball',obj.sus.water,...
              'corpus_collosum',obj.sus.WM,'special_region_frontal_lobes',obj.sus.GM,'cerebral_falx',obj.sus.GM,'temporal_lobes',obj.sus.GM,...
              'fourth_ventricle',obj.sus.CSF,'frontal_portion_eyes',obj.sus.water,'parietal_lobes',obj.sus.GM,'amygdala',obj.sus.GM,...
              'eye',obj.sus.water,'globus_pallidus',obj.sus.globus_pallidus,'lens',obj.sus.water,'cerebral_aquaduct',obj.sus.CSF,...
              'lateral_ventricles',obj.sus.CSF,'prefrontal_lobes',obj.sus.GM,'teeth',obj.sus.teeth,'sigmoid_sinus',obj.sus.venous_blood);

          % create suscceptibility phantom
          obj.volume = zubal_phantom * 1e6;
          
          obj.volume(zubal_phantom==zubal_id.outside_phantom) = zubal_sus.outside_phantom;
          obj.volume(zubal_phantom==zubal_id.skin) = zubal_sus.skin;
          obj.volume(zubal_phantom==zubal_id.cerebral_fluid) = zubal_sus.cerebral_fluid;
          obj.volume(zubal_phantom==zubal_id.spinal_cord) = zubal_sus.spinal_cord;
          obj.volume(zubal_phantom==zubal_id.skull) = zubal_sus.skull;
          obj.volume(zubal_phantom==zubal_id.spine) = zubal_sus.spine;
          obj.volume(zubal_phantom==zubal_id.dens_of_axis) = zubal_sus.dens_of_axis;
          obj.volume(zubal_phantom==zubal_id.jaw_bone) = zubal_sus.jaw_bone;
          obj.volume(zubal_phantom==zubal_id.parotid_gland) = zubal_sus.parotid_gland;
          obj.volume(zubal_phantom==zubal_id.skeletal_muscle) = zubal_sus.skeletal_muscle;
          obj.volume(zubal_phantom==zubal_id.lacrimal_glands) = zubal_sus.lacrimal_glands;
          obj.volume(zubal_phantom==zubal_id.spinal_canal) = zubal_sus.spinal_canal;
          obj.volume(zubal_phantom==zubal_id.hard_palate) = zubal_sus.hard_palate;
          obj.volume(zubal_phantom==zubal_id.cerebellum) = zubal_sus.cerebellum;
          obj.volume(zubal_phantom==zubal_id.tongue) = zubal_sus.tongue;
          obj.volume(zubal_phantom==zubal_id.pharynx) = zubal_sus.pharynx;
          obj.volume(zubal_phantom==zubal_id.esophagus) = zubal_sus.esophagus;
          obj.volume(zubal_phantom==zubal_id.horn_of_mandible) = zubal_sus.horn_of_mandible;
          obj.volume(zubal_phantom==zubal_id.nasal_septum) = zubal_sus.nasal_septum;
          obj.volume(zubal_phantom==zubal_id.white_matter) = zubal_sus.white_matter;
          obj.volume(zubal_phantom==zubal_id.superior_sagittal_sinus) = zubal_sus.superior_sagittal_sinus;
          obj.volume(zubal_phantom==zubal_id.medulla_oblongata) = zubal_sus.medulla_oblongata;
          obj.volume(zubal_phantom==zubal_id.blood_pool) = zubal_sus.blood_pool;
          obj.volume(zubal_phantom==zubal_id.frontal_lobes) = zubal_sus.frontal_lobes;
          obj.volume(zubal_phantom==zubal_id.bone_marrow) = zubal_sus.bone_marrow;
          obj.volume(zubal_phantom==zubal_id.pons) = zubal_sus.pons;
          obj.volume(zubal_phantom==zubal_id.third_ventricle) = zubal_sus.third_ventricle;
          obj.volume(zubal_phantom==zubal_id.trachea) = zubal_sus.trachea;
          obj.volume(zubal_phantom==zubal_id.cartilage) = zubal_sus.cartilage;
          obj.volume(zubal_phantom==zubal_id.occipital_lobes) = zubal_sus.occipital_lobes;
          obj.volume(zubal_phantom==zubal_id.hippocampus) = zubal_sus.hippocampus;
          obj.volume(zubal_phantom==zubal_id.pituitary_gland) = zubal_sus.pituitary_gland;
          obj.volume(zubal_phantom==zubal_id.fat1) = zubal_sus.fat1;
          obj.volume(zubal_phantom==zubal_id.fat2) = zubal_sus.fat2;
          obj.volume(zubal_phantom==zubal_id.ear_bones) = zubal_sus.ear_bones;
          obj.volume(zubal_phantom==zubal_id.turbinates) = zubal_sus.turbinates;
          obj.volume(zubal_phantom==zubal_id.caudate_nucleus) = zubal_sus.caudate_nucleus;
          obj.volume(zubal_phantom==zubal_id.zygoma) = zubal_sus.zygoma;
          obj.volume(zubal_phantom==zubal_id.insula_cortex) = zubal_sus.insula_cortex;
          obj.volume(zubal_phantom==zubal_id.sinuses_mouth_cavity) = zubal_sus.sinuses_mouth_cavity;
          obj.volume(zubal_phantom==zubal_id.putamen) = zubal_sus.putamen;
          obj.volume(zubal_phantom==zubal_id.optic_nerve) = zubal_sus.optic_nerve;
          obj.volume(zubal_phantom==zubal_id.internal_capsule) = zubal_sus.internal_capsule;
          obj.volume(zubal_phantom==zubal_id.septum_pellucidium) = zubal_sus.septum_pellucidium;
          obj.volume(zubal_phantom==zubal_id.thalamus) = zubal_sus.thalamus;
          obj.volume(zubal_phantom==zubal_id.eyeball) = zubal_sus.eyeball;
          obj.volume(zubal_phantom==zubal_id.corpus_collosum) = zubal_sus.corpus_collosum;
          obj.volume(zubal_phantom==zubal_id.special_region_frontal_lobes) = zubal_sus.special_region_frontal_lobes;
          obj.volume(zubal_phantom==zubal_id.cerebral_falx) = zubal_sus.cerebral_falx;
          obj.volume(zubal_phantom==zubal_id.temporal_lobes) = zubal_sus.temporal_lobes;
          obj.volume(zubal_phantom==zubal_id.fourth_ventricle) = zubal_sus.fourth_ventricle;
          obj.volume(zubal_phantom==zubal_id.frontal_portion_eyes) = zubal_sus.frontal_portion_eyes;
          obj.volume(zubal_phantom==zubal_id.parietal_lobes) = zubal_sus.parietal_lobes;
          obj.volume(zubal_phantom==zubal_id.amygdala) = zubal_sus.amygdala;
          obj.volume(zubal_phantom==zubal_id.eye) = zubal_sus.eye;
          obj.volume(zubal_phantom==zubal_id.globus_pallidus) = zubal_sus.globus_pallidus;
          obj.volume(zubal_phantom==zubal_id.lens) = zubal_sus.lens;
          obj.volume(zubal_phantom==zubal_id.cerebral_aquaduct) = zubal_sus.cerebral_aquaduct;
          obj.volume(zubal_phantom==zubal_id.lateral_ventricles) = zubal_sus.lateral_ventricles;
          obj.volume(zubal_phantom==zubal_id.prefrontal_lobes) = zubal_sus.prefrontal_lobes;
          obj.volume(zubal_phantom==zubal_id.teeth) = zubal_sus.teeth;
          obj.volume(zubal_phantom==zubal_id.sigmoid_sinus) = zubal_sus.sigmoid_sinus;
          
       end
       
    end
    
 end