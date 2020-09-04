% An updated, dedicated MRI head phantom made available by techniques
% described in Zubal IG, Harrell CR, Smith EO, Rattner Z, Gindi GR, Hoffer
% PB "Computerized Three-dimensional Segmented Human Anatomy" Medical
% Physics 21(2), Feb 1994, p.299-302"
% - Zubal IG, Harrell CR, Smith EO, Rattner Z, Gindi GR, Hoffer PB
% "Computerized Three-dimensional Segmented Human Anatomy" Medical Physics
% 21(2), Feb 1994, p.299-302.

%head_phantom = double(niftiread('../det_head_u2.nii'));
%head_phantom = imrotate(head_phantom,90);

head_phantom = double(niftiread('zubal_EAO.nii'));

% note: CSF in the midbrain appears to be incorrectly labeled as bone
% marrow, I am keeping the labels here unchanged, but the T2* and
% susceptibility of "bone marrow" are set to that of CSF.
phantom_id = struct('outside_phantom',0,'skin',1,'cerebral_fluid',2,'spinal_cord',3,'skull',4,'spine',5,'dens_of_axis',70, ...
    'jaw_bone',71,'parotid_gland',72,'skeletal_muscle',9,'lacrimal_glands',74,'spinal_canal',75,'hard_palate',76, ...
    'cerebellum',77,'tongue',78,'pharynx',15,'esophagus',16,'horn_of_mandible',81,'nasal_septum',82,'white_matter',83, ...
    'superior_sagittal_sinus',84,'medulla_oblongata',85,'blood_pool',23,'frontal_lobes',89,'bone_marrow',26, ...
    'pons',91,'third_ventricle',92,'trachea',29,'cartilage',30,'occipital_lobes',95,'hippocampus',96,'pituitary_gland',97, ...
    'fat1',98,'fat2',22,'ear_bones',99,'turbinates',100,'caudate_nucleus',101,'zygoma',102,'insula_cortex',103,'sinuses_mouth_cavity',104, ...
    'putamen',105,'optic_nerve',106,'internal_capsule',107,'septum_pellucidium',108,'thalamus',109,'eyeball',110,'corpus_collosum',111, ...
    'special_region_frontal_lobes',112,'cerebral_falx',113,'temporal_lobes',114,'fourth_ventricle',115,'frontal_portion_eyes',116, ...
    'parietal_lobes',117,'amygdala',118,'eye',119,'globus_pallidus',120,'lens',121,'cerebral_aquaduct',122,'lateral_ventricles',123, ...
    'prefrontal_lobes',124,'teeth',125,'sigmoid_sinus',126);

% blood pd from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3023928/
% pd of fat set to 80, which is a guess
% pd of muscle set to 65, which is a guess
% pd of cartilage set to 40, which is a guess
tissue_pd = struct('teeth',20,'bone',20,'GM',82,'WM',70,'CSF',100,...
    'air',0,'water',100,'fat',80,'muscle',65,'blood',83,'cartilage',40);

phantom_pd = struct('outside_phantom',tissue_pd.air,'skin',tissue_pd.water,'cerebral_fluid',tissue_pd.CSF,'spinal_cord',(tissue_pd.GM+tissue_pd.WM)/2,...
    'skull',tissue_pd.bone,'spine',tissue_pd.bone,'dens_of_axis',tissue_pd.bone,'jaw_bone',tissue_pd.bone,'parotid_gland',tissue_pd.water,...
    'skeletal_muscle',tissue_pd.muscle,'lacrimal_glands',tissue_pd.water,'spinal_canal',tissue_pd.CSF,'hard_palate',tissue_pd.bone,...
    'cerebellum',tissue_pd.GM,'tongue',tissue_pd.muscle,'pharynx',tissue_pd.air,'esophagus',tissue_pd.water,'horn_of_mandible',tissue_pd.bone,...
    'nasal_septum',tissue_pd.bone,'white_matter',tissue_pd.WM,'superior_sagittal_sinus',tissue_pd.blood,...
    'medulla_oblongata',(tissue_pd.GM+tissue_pd.WM)/2,'fat1',tissue_pd.fat,'fat2',tissue_pd.fat,'blood_pool',tissue_pd.blood,...
    'frontal_lobes',tissue_pd.GM,'bone_marrow',tissue_pd.CSF,'pons',tissue_pd.WM,'third_ventricle',tissue_pd.CSF,'trachea',tissue_pd.air,...
    'cartilage',tissue_pd.cartilage,'occipital_lobes',tissue_pd.GM,'hippocampus',tissue_pd.GM,'pituitary_gland',tissue_pd.water,'ear_bones',tissue_pd.bone,...
    'turbinates',(tissue_pd.bone+tissue_pd.water)/2,'caudate_nucleus',tissue_pd.GM,'zygoma',tissue_pd.bone,...
    'insula_cortex',tissue_pd.GM,'sinuses_mouth_cavity',tissue_pd.air,'putamen',tissue_pd.GM,'optic_nerve',tissue_pd.WM,...
    'internal_capsule',tissue_pd.WM,'septum_pellucidium',(tissue_pd.WM+tissue_pd.GM)/2,'thalamus',tissue_pd.GM,'eyeball',tissue_pd.water,...
    'corpus_collosum',tissue_pd.WM,'special_region_frontal_lobes',tissue_pd.GM,'cerebral_falx',tissue_pd.GM,'temporal_lobes',tissue_pd.GM,...
    'fourth_ventricle',tissue_pd.CSF,'frontal_portion_eyes',tissue_pd.water,'parietal_lobes',tissue_pd.GM,'amygdala',tissue_pd.GM,...
    'eye',tissue_pd.water,'globus_pallidus',tissue_pd.GM,'lens',tissue_pd.water,'cerebral_aquaduct',tissue_pd.CSF,...
    'lateral_ventricles',tissue_pd.CSF,'prefrontal_lobes',tissue_pd.GM,'teeth',tissue_pd.teeth,'sigmoid_sinus',tissue_pd.blood);

% absolute susceptibility [ppm] values are from Buch et al. MRM 73:2185?2194 (2015)
% susceptibility of fat is from https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.27640
% susceptibility of muscle is from https://www.researchgate.net/figure/Magnetic-susceptibility-of-different-matter_tbl1_288160127
% susceptibility of venous and arterial blood from https://www.researchgate.net/figure/Magnetic-susceptibility-of-different-matter_tbl1_288160127
% susceptibility of cartilage is from https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.26596
tissue_sus = struct('sinus',0.2e-6,'teeth',-12.3e-6,'bone',-11.1e-6,'thalamus',-9.02e-6, ...
    'caudate_nucleus',-8.99e-6,'putamen',-8.96e-6,'globus_pallidus',-8.89e-6,'GM',-9.03e-6,'WM',-9.0830e-6,'CSF',-9.05e-6,...
    'air',0.35e-6,'water',-9.05e-6,'fat',-8.39e-6,'muscle',-9.03e-6,'venous_blood',-7.8e-6,'arterial_blood',-9.3e-6,...
    'cartilage',-9.055e-6);

phantom_sus = struct('outside_phantom',tissue_sus.air,'skin',tissue_sus.water,'cerebral_fluid',tissue_sus.CSF,'spinal_cord',(tissue_sus.GM+tissue_sus.WM)/2,...
    'skull',tissue_sus.bone,'spine',tissue_sus.bone,'dens_of_axis',tissue_sus.bone,'jaw_bone',tissue_sus.bone,'parotid_gland',tissue_sus.water,...
    'skeletal_muscle',tissue_sus.muscle,'lacrimal_glands',tissue_sus.water,'spinal_canal',tissue_sus.CSF,'hard_palate',tissue_sus.bone,...
    'cerebellum',tissue_sus.GM,'tongue',tissue_sus.muscle,'pharynx',tissue_sus.sinus,'esophagus',tissue_sus.water,'horn_of_mandible',tissue_sus.bone,...
    'nasal_septum',tissue_sus.bone,'white_matter',tissue_sus.WM,'superior_sagittal_sinus',tissue_sus.venous_blood,...
    'medulla_oblongata',(tissue_sus.GM+tissue_sus.WM)/2,'fat1',tissue_sus.fat,'fat2',tissue_sus.fat,'blood_pool',(tissue_sus.venous_blood+tissue_sus.arterial_blood)/2,...
    'frontal_lobes',tissue_sus.GM,'bone_marrow',tissue_sus.CSF,'pons',tissue_sus.WM,'third_ventricle',tissue_sus.CSF,'trachea',tissue_sus.sinus,...
    'cartilage',tissue_sus.cartilage,'occipital_lobes',tissue_sus.GM,'hippocampus',tissue_sus.GM,'pituitary_gland',tissue_sus.water,'ear_bones',tissue_sus.bone,...
    'turbinates',(tissue_sus.bone+tissue_sus.water)/2,'caudate_nucleus',tissue_sus.caudate_nucleus,'zygoma',tissue_sus.bone,...
    'insula_cortex',tissue_sus.GM,'sinuses_mouth_cavity',tissue_sus.sinus,'putamen',tissue_sus.putamen,'optic_nerve',tissue_sus.WM,...
    'internal_capsule',tissue_sus.WM,'septum_pellucidium',(tissue_sus.WM+tissue_sus.GM)/2,'thalamus',tissue_sus.thalamus,'eyeball',tissue_sus.water,...
    'corpus_collosum',tissue_sus.WM,'special_region_frontal_lobes',tissue_sus.GM,'cerebral_falx',tissue_sus.GM,'temporal_lobes',tissue_sus.GM,...
    'fourth_ventricle',tissue_sus.CSF,'frontal_portion_eyes',tissue_sus.water,'parietal_lobes',tissue_sus.GM,'amygdala',tissue_sus.GM,...
    'eye',tissue_sus.water,'globus_pallidus',tissue_sus.globus_pallidus,'lens',tissue_sus.water,'cerebral_aquaduct',tissue_sus.CSF,...
    'lateral_ventricles',tissue_sus.CSF,'prefrontal_lobes',tissue_sus.GM,'teeth',tissue_sus.teeth,'sigmoid_sinus',tissue_sus.venous_blood);

% T2* values at 3T from 
% https://index.mirasmart.com/ISMRM2019/PDFfiles/4509.html
% https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.26809
% T2* CSF is set to 50% of the T2 of CSF (http://mri-q.com/why-is-t1--t2.html)
% T2* water is set to 50% of the T2 of water (http://mri-q.com/why-is-t1--t2.html)
% T2* of teeth is set to 50% of the T2 of teeth (https://www.karger.com/Article/Abstract/501901)
% T2* of fat is set to 50% of the T2 of fat (http://mri-q.com/why-is-t1--t2.html)
% T2* of muscle is set to 50% of the T2 of muscle (http://mri-q.com/why-is-t1--t2.html)
% T2* of venous and arterial blood is from:
% https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.21342, table 1 using
% (1-Ya) = 0.02 and (1-Yv) = 0.39 and a Hct level of 0.44
% T2* of the "frontal portion eyes" is set to the T2 of the cornea (~ 50
% ms) from https://onlinelibrary.wiley.com/doi/full/10.1002/jmri.21017
% -> replaced with T2* of cornea from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4160095/
% T2* of the "eyeball" is set to the T2 of the Chorioretina (~ 75 ms) from
% https://onlinelibrary.wiley.com/doi/full/10.1002/jmri.21017
% -> replaced with T2* of the sclera from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4160095/
% T2* of the "eye" is set to the 50% of T2 of the vitrous humour (~ 750 ms) from https://onlinelibrary.wiley.com/doi/full/10.1002/jmri.21017
% T2* of the "lens" is set to the 50% of T2 of the lens nucleus (~ 1130 ms) from https://onlinelibrary.wiley.com/doi/full/10.1002/jmri.21017
% T2* of the "globus pallidus" is from https://pubs.rsna.org/doi/full/10.1148/radiol.2522081399
% T2* of cartilage is from https://openmedicinejournal.com/VOLUME/5/PAGE/119/FULLTEXT/
% T2* of the skin is set to 50% of T2 of the skin at 1.5T from https://www.sciencedirect.com/science/article/pii/S0022202X9190213A
% T2* of the esophagus is from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5447643/
% T2* of the parotid gland is set to 50% of T2 of the parotid gland from https://pubmed.ncbi.nlm.nih.gov/23166360/
% T2* of the pituitary gland is set to 50% of T2 of the pituitary gland from https://www.ajronline.org/doi/full/10.2214/ajr.175.6.1751567
tissue_t2s = struct('sinus',0,'teeth',0.5*150e-3,'bone',3.3e3,'thalamus',58e-3, ...
    'caudate_nucleus',58e-3,'putamen',40e-3,'globus_pallidus',27e-3,'GM',66e-3,'WM',61e-3,'CSF',0.5*2,'air',0,'water',0.5*2, ...
    'fat',0.5*70e-3,'muscle',0.5*50e-3,'venous_blood',55e-3,'arterial_blood',20e-3,'frontal_portion_eyes',25e-3,'eyeball',10e-3,...
    'eye',0.5*750e-3,'lens',0.5*1130e-3,'cartilage',25e-3,'skin',0.5*30e-3,'esophagus',17e-3,'parotid_gland',0.5*120e-3,'pituitary_gland',0.5*90e-3);

% notes: 
% turbinates are made of bone and soft tissue, so the average of bone and water is used
% T2* of cerebral_falx is set to that of GM
% T2* of the lacrimal glands is set to that of the parotid gland
phantom_t2s = struct('outside_phantom',tissue_t2s.air,'skin',tissue_t2s.skin,'cerebral_fluid',tissue_t2s.CSF,'spinal_cord',(tissue_t2s.GM+tissue_t2s.WM)/2,...
    'skull',tissue_t2s.bone,'spine',tissue_t2s.bone,'dens_of_axis',tissue_t2s.bone,'jaw_bone',tissue_t2s.bone,'parotid_gland',tissue_t2s.parotid_gland,...
    'skeletal_muscle',tissue_t2s.muscle,'lacrimal_glands',tissue_t2s.parotid_gland,'spinal_canal',tissue_t2s.CSF,'hard_palate',tissue_t2s.bone,...
    'cerebellum',tissue_t2s.GM,'tongue',tissue_t2s.muscle,'pharynx',tissue_t2s.sinus,'esophagus',tissue_t2s.esophagus,'horn_of_mandible',tissue_t2s.bone,...
    'nasal_septum',tissue_t2s.bone,'white_matter',tissue_t2s.WM,'superior_sagittal_sinus',tissue_t2s.venous_blood,...
    'medulla_oblongata',(tissue_t2s.GM+tissue_t2s.WM)/2,'fat1',tissue_t2s.fat,'fat2',tissue_t2s.fat,'blood_pool',(tissue_t2s.venous_blood+tissue_t2s.arterial_blood)/2,...
    'frontal_lobes',tissue_t2s.GM,'bone_marrow',tissue_t2s.CSF,'pons',tissue_t2s.WM,'third_ventricle',tissue_t2s.CSF,'trachea',tissue_t2s.sinus,...
    'cartilage',tissue_t2s.cartilage,'occipital_lobes',tissue_t2s.GM,'hippocampus',tissue_t2s.GM,'pituitary_gland',tissue_t2s.pituitary_gland,'ear_bones',tissue_t2s.bone,...
    'turbinates',(tissue_t2s.bone+tissue_t2s.water)/2,'caudate_nucleus',tissue_t2s.caudate_nucleus,'zygoma',tissue_t2s.bone,...
    'insula_cortex',tissue_t2s.GM,'sinuses_mouth_cavity',tissue_t2s.sinus,'putamen',tissue_t2s.putamen,'optic_nerve',tissue_t2s.WM,...
    'internal_capsule',tissue_t2s.WM,'septum_pellucidium',(tissue_t2s.WM+tissue_t2s.GM)/2,'thalamus',tissue_t2s.thalamus,'eyeball',tissue_t2s.frontal_portion_eyes,...
    'corpus_collosum',tissue_t2s.WM,'special_region_frontal_lobes',tissue_t2s.GM,'cerebral_falx',tissue_t2s.GM,'temporal_lobes',tissue_t2s.GM,...
    'fourth_ventricle',tissue_t2s.CSF,'frontal_portion_eyes',tissue_t2s.frontal_portion_eyes,'parietal_lobes',tissue_t2s.GM,'amygdala',tissue_t2s.GM,...
    'eye',tissue_t2s.eye,'globus_pallidus',tissue_t2s.globus_pallidus,'lens',tissue_t2s.lens,'cerebral_aquaduct',tissue_t2s.CSF,...
    'lateral_ventricles',tissue_t2s.CSF,'prefrontal_lobes',tissue_t2s.GM,'teeth',tissue_t2s.teeth,'sigmoid_sinus',tissue_t2s.venous_blood);

% create suscceptibility phantom
head_sus = head_phantom;
head_sus(head_phantom==phantom_id.outside_phantom) = phantom_sus.outside_phantom;
head_sus(head_phantom==phantom_id.skin) = phantom_sus.skin;
head_sus(head_phantom==phantom_id.cerebral_fluid) = phantom_sus.cerebral_fluid;
head_sus(head_phantom==phantom_id.spinal_cord) = phantom_sus.spinal_cord;
head_sus(head_phantom==phantom_id.skull) = phantom_sus.skull;
head_sus(head_phantom==phantom_id.spine) = phantom_sus.spine;
head_sus(head_phantom==phantom_id.dens_of_axis) = phantom_sus.dens_of_axis;
head_sus(head_phantom==phantom_id.jaw_bone) = phantom_sus.jaw_bone;
head_sus(head_phantom==phantom_id.parotid_gland) = phantom_sus.parotid_gland;
head_sus(head_phantom==phantom_id.skeletal_muscle) = phantom_sus.skeletal_muscle;
head_sus(head_phantom==phantom_id.lacrimal_glands) = phantom_sus.lacrimal_glands;
head_sus(head_phantom==phantom_id.spinal_canal) = phantom_sus.spinal_canal;
head_sus(head_phantom==phantom_id.hard_palate) = phantom_sus.hard_palate;
head_sus(head_phantom==phantom_id.cerebellum) = phantom_sus.cerebellum;
head_sus(head_phantom==phantom_id.tongue) = phantom_sus.tongue;
head_sus(head_phantom==phantom_id.pharynx) = phantom_sus.pharynx;
head_sus(head_phantom==phantom_id.esophagus) = phantom_sus.esophagus;
head_sus(head_phantom==phantom_id.horn_of_mandible) = phantom_sus.horn_of_mandible;
head_sus(head_phantom==phantom_id.nasal_septum) = phantom_sus.nasal_septum;
head_sus(head_phantom==phantom_id.white_matter) = phantom_sus.white_matter;
head_sus(head_phantom==phantom_id.superior_sagittal_sinus) = phantom_sus.superior_sagittal_sinus;
head_sus(head_phantom==phantom_id.medulla_oblongata) = phantom_sus.medulla_oblongata;
head_sus(head_phantom==phantom_id.blood_pool) = phantom_sus.blood_pool;
head_sus(head_phantom==phantom_id.frontal_lobes) = phantom_sus.frontal_lobes;
head_sus(head_phantom==phantom_id.bone_marrow) = phantom_sus.bone_marrow;
head_sus(head_phantom==phantom_id.pons) = phantom_sus.pons;
head_sus(head_phantom==phantom_id.third_ventricle) = phantom_sus.third_ventricle;
head_sus(head_phantom==phantom_id.trachea) = phantom_sus.trachea;
head_sus(head_phantom==phantom_id.cartilage) = phantom_sus.cartilage;
head_sus(head_phantom==phantom_id.occipital_lobes) = phantom_sus.occipital_lobes;
head_sus(head_phantom==phantom_id.hippocampus) = phantom_sus.hippocampus;
head_sus(head_phantom==phantom_id.pituitary_gland) = phantom_sus.pituitary_gland;
head_sus(head_phantom==phantom_id.fat1) = phantom_sus.fat1;
head_sus(head_phantom==phantom_id.fat2) = phantom_sus.fat2;
head_sus(head_phantom==phantom_id.ear_bones) = phantom_sus.ear_bones;
head_sus(head_phantom==phantom_id.turbinates) = phantom_sus.turbinates;
head_sus(head_phantom==phantom_id.caudate_nucleus) = phantom_sus.caudate_nucleus;
head_sus(head_phantom==phantom_id.zygoma) = phantom_sus.zygoma;
head_sus(head_phantom==phantom_id.insula_cortex) = phantom_sus.insula_cortex;
head_sus(head_phantom==phantom_id.sinuses_mouth_cavity) = phantom_sus.sinuses_mouth_cavity;
head_sus(head_phantom==phantom_id.putamen) = phantom_sus.putamen;
head_sus(head_phantom==phantom_id.optic_nerve) = phantom_sus.optic_nerve;
head_sus(head_phantom==phantom_id.internal_capsule) = phantom_sus.internal_capsule;
head_sus(head_phantom==phantom_id.septum_pellucidium) = phantom_sus.septum_pellucidium;
head_sus(head_phantom==phantom_id.thalamus) = phantom_sus.thalamus;
head_sus(head_phantom==phantom_id.eyeball) = phantom_sus.eyeball;
head_sus(head_phantom==phantom_id.corpus_collosum) = phantom_sus.corpus_collosum;
head_sus(head_phantom==phantom_id.special_region_frontal_lobes) = phantom_sus.special_region_frontal_lobes;
head_sus(head_phantom==phantom_id.cerebral_falx) = phantom_sus.cerebral_falx;
head_sus(head_phantom==phantom_id.temporal_lobes) = phantom_sus.temporal_lobes;
head_sus(head_phantom==phantom_id.fourth_ventricle) = phantom_sus.fourth_ventricle;
head_sus(head_phantom==phantom_id.frontal_portion_eyes) = phantom_sus.frontal_portion_eyes;
head_sus(head_phantom==phantom_id.parietal_lobes) = phantom_sus.parietal_lobes;
head_sus(head_phantom==phantom_id.amygdala) = phantom_sus.amygdala;
head_sus(head_phantom==phantom_id.eye) = phantom_sus.eye;
head_sus(head_phantom==phantom_id.globus_pallidus) = phantom_sus.globus_pallidus;
head_sus(head_phantom==phantom_id.lens) = phantom_sus.lens;
head_sus(head_phantom==phantom_id.cerebral_aquaduct) = phantom_sus.cerebral_aquaduct;
head_sus(head_phantom==phantom_id.lateral_ventricles) = phantom_sus.lateral_ventricles;
head_sus(head_phantom==phantom_id.prefrontal_lobes) = phantom_sus.prefrontal_lobes;
head_sus(head_phantom==phantom_id.teeth) = phantom_sus.teeth;
head_sus(head_phantom==phantom_id.sigmoid_sinus) = phantom_sus.sigmoid_sinus;

%nii_vol = make_nii(flip(head_sus,3));
nii_vol = make_nii(head_sus);
save_nii(nii_vol, ['zubal_head_sus' '.nii']);

Bdz = fourier_based_field_est(head_sus, [1.1 1.1 1.4], [256 256 128]);

%nii_vol = make_nii((flip(real(Bdz),3)));
nii_vol = make_nii(real(Bdz));
save_nii(nii_vol,'Bdz_zubal_head_sus.nii');