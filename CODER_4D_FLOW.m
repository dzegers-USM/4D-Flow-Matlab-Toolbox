%codegen

function CODER_4D_FLOW(handles) %#codegen
    % handles = struct('id_unwrappping', 0, ...
    %     'VENC', 0, ...
    %     'voxel_MR', [0 0 0], ...
    %     'type', 'UNK');
    
    % THIS FUNCTION LOAD THE INFORMATION TO BE PROCESSED 
    handles.id_unwrappping = 0;
    %set(handles.uipanel1,'Visible','off')
    %cla(handles.uipanel1,'reset');

    % TODO: Load these folders
    % path(path,['IO_CODES',filesep]) % cambiar
    % path(path,'iso2mesh/')

    % TODO: Ask user for folder and read dynamically
   
    % reading files names matlab, par-rec, dcm, dat
    % *.mat
    % *.PAR
    % *.REC
    % *.dcm
    % *.dat

    files_names_mat = {'data/data.mat' 'data.mat'};

    files_names_par = [];
    files_names_rec = [];
    files_names_dcm = [];
    
    cont_mat = 1;
    cont_par = 1;
    cont_rec = 1;
    cont_dcm = 1;
    
    cont_dat_vd_1 = 1;
    cont_dat_vd_2 = 1;
    cont_dat_vd_3 = 1;
    
    files_names_dat_vd_1 = [];
    files_names_dat_vd_2 = [];
    files_names_dat_vd_3 = [];
    files_names_txt_hd = [];
    files_names_dat_CD = [];
    
    % reading other files
    % for k = 1 : numberOfFolders
    % 
    %     if isempty(files_names_mat)==1
    % 
    % 
    %         thisFolder = listOfFolderNames{k};
    % 
    %         filePattern_par = sprintf('%s/*.PAR', thisFolder);
    %         baseFileNames_par = dir(filePattern_par);
    %         numberOfFiles_par = length(baseFileNames_par);
    % 
    %         filePattern_rec = sprintf('%s/*.REC', thisFolder);
    %         baseFileNames_rec = dir(filePattern_rec);
    %         numberOfFiles_rec = length(baseFileNames_rec);
    % 
    %         filePattern_dcm = sprintf('%s/*.dcm', thisFolder);
    %         baseFileNames_dcm = dir(filePattern_dcm);
    %         numberOfFiles_dcm = length(baseFileNames_dcm);
    % 
    %         filePattern_dat = sprintf('%s/*.dat', thisFolder);
    %         baseFileNames_dat = dir(filePattern_dat);
    %         numberOfFiles_dat = length(baseFileNames_dat);
    % 
    %         if numberOfFiles_par>=1 
    %             for n = 1:numberOfFiles_par
    %                 files_names_par{cont_par} = [thisFolder,filesep,baseFileNames_par(n).name];
    %                 cont_par = cont_par + 1;
    %             end
    %         end
    % 
    %         if numberOfFiles_rec>=1 
    %             for n = 1:numberOfFiles_rec
    %                 files_names_rec{cont_rec} = [thisFolder,filesep,baseFileNames_rec(n).name];
    %                 cont_rec = cont_rec + 1;
    %             end
    %         end
    % 
    %         if numberOfFiles_dcm>=1 
    %             for n = 1:numberOfFiles_dcm
    %                 files_names_dcm{cont_dcm} = [thisFolder,filesep,baseFileNames_dcm(n).name];
    %                 cont_dcm = cont_dcm + 1;
    %             end
    %         end
    % 
    %         if numberOfFiles_dat>=1 
    %             files_names_txt_hd = [thisFolder,filesep,'pcvipr_header.txt'];
    %             files_names_dat_CD = [thisFolder,filesep,'CD.dat'];
    %             for n = 1:numberOfFiles_dat
    %                 if strncmp(flip(baseFileNames_dat(n).name),'tad.1_dv_',9)==1 && strncmp(baseFileNames_dat(n).name,'ph_',3)==1
    %                     files_names_dat_vd_1{cont_dat_vd_1} = [thisFolder,filesep,baseFileNames_dat(n).name];
    %                     cont_dat_vd_1 = cont_dat_vd_1 + 1;
    %                 elseif  strncmp(flip(baseFileNames_dat(n).name),'tad.2_dv_',9)==1 && strncmp(baseFileNames_dat(n).name,'ph_',3)==1
    %                     files_names_dat_vd_2{cont_dat_vd_2} = [thisFolder,filesep,baseFileNames_dat(n).name];
    %                     cont_dat_vd_2 = cont_dat_vd_2 + 1;
    %                 elseif  strncmp(flip(baseFileNames_dat(n).name),'tad.3_dv_',9)==1 && strncmp(baseFileNames_dat(n).name,'ph_',3)==1
    %                     files_names_dat_vd_3{cont_dat_vd_3} = [thisFolder,filesep,baseFileNames_dat(n).name];
    %                     cont_dat_vd_3 = cont_dat_vd_3 + 1;
    %                 end
    %             end
    %         end
    %     end
    % end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % READING MATLAB FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(files_names_mat)==0
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % remove hidden files
        % mat_temp = [];
        % cont=1;
        % for ff = 1:length(files_names_mat)
        %     if isempty(strfind(files_names_mat{ff},'._'))==1
        %         mat_temp{cont}=files_names_mat{ff};
        %         cont = cont + 1;
        %     end
        % end
        % files_names_mat = mat_temp;
        
        
        % loading structure data
        % h = msgbox({'Please wait ...','Loading structure data ...'});
        
        fullFileName = files_names_mat{1};
        fileData = load(fullFileName);
        data = fileData.data;
        
        handles.VENC = data.VENC;
        handles.voxel_MR = data.voxel_MR;
        handles.heart_rate = data.heart_rate;
        % handles.type = data.type;
        handles.MR_FFE_FH = data.MR_FFE_FH;
        handles.MR_FFE_AP = data.MR_FFE_AP;
        handles.MR_FFE_RL = data.MR_FFE_RL;
        handles.MR_PCA_FH = data.MR_PCA_FH;
        handles.MR_PCA_AP = data.MR_PCA_AP;
        handles.MR_PCA_RL = data.MR_PCA_RL;

        [~,~,~,d] = size(data.MR_FFE_FH);
        IPCMRA = sqrt((1/d)*sum( ((data.MR_FFE_FH).^2).*((data.MR_PCA_FH).^2 + (data.MR_PCA_AP).^2 + (data.MR_PCA_RL).^2),4));
        handles.IPCMRA = IPCMRA;

        [a,b,c,d] = size(data.MR_FFE_FH);
        handles.a = a;
        handles.b = b;
        handles.c = c;
        handles.d = d;
        handles
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load data offset error JSOTELO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.type = 'MAT';

        input = struct( ...
            'VENC', handles.VENC, ...
            'voxel_MR', handles.voxel_MR, ...
            'heart_rate', handles.heart_rate, ...
            'type', handles.type, ...
            'MR_FFE_FH', handles.MR_FFE_FH, ...
            'MR_FFE_AP', handles.MR_FFE_AP, ...
            'MR_FFE_RL', handles.MR_FFE_RL, ...
            'MR_PCA_FH', handles.MR_PCA_FH, ...
            'MR_PCA_AP', handles.MR_PCA_AP, ...
            'MR_PCA_RL', handles.MR_PCA_RL, ...
            'IPCMRA', handles.IPCMRA, ...
            'id', 0, ...
            'id_while', 0);
        % input.VENC = handles.VENC;
        % input.voxel_MR = handles.voxel_MR;
        % input.heart_rate = handles.heart_rate
        % input.type = handles.type;
        % input.MR_FFE_FH = handles.MR_FFE_FH;
        % input.MR_FFE_AP = handles.MR_FFE_AP;
        % input.MR_FFE_RL = handles.MR_FFE_RL;
        % input.MR_PCA_FH = handles.MR_PCA_FH;
        % input.MR_PCA_AP = handles.MR_PCA_AP;
        % input.MR_PCA_RL = handles.MR_PCA_RL;
        % input.id = 0;
        % input.id_while = 0;  

        % TODO: Err and noise masking
        % id_while = 0;
        % while(1)
        %     while(id_while == 0)
        %         OFFSET_ERR_AND_NOISE_MASKING(input)
        %         input.VENC = getappdata(0,'VENC');
        %         input.voxel_MR = getappdata(0,'voxel_MR');
        %         input.heart_rate = getappdata(0,'heart_rate');
        %         input.type = getappdata(0,'type');
        %         input.MR_FFE_FH = getappdata(0,'MR_FFE_FH');
        %         input.MR_FFE_AP = getappdata(0,'MR_FFE_AP');
        %         input.MR_FFE_RL = getappdata(0,'MR_FFE_RL');
        %         input.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
        %         input.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
        %         input.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
        %         input.IPCMRA = getappdata(0,'IPCMRA');
        %         input.id = getappdata(0,'id');
        %         id_while = getappdata(0,'id_while');
        %         pause(0.05)
        %     end
        %     handles.VENC = getappdata(0,'VENC');
        %     handles.voxel_MR = getappdata(0,'voxel_MR');
        %     handles.heart_rate = getappdata(0,'heart_rate');
        %     %handles.type = getappdata(0,'type');
        %     handles.MR_FFE_FH = getappdata(0,'MR_FFE_FH');
        %     handles.MR_FFE_AP = getappdata(0,'MR_FFE_AP');
        %     handles.MR_FFE_RL = getappdata(0,'MR_FFE_RL');
        %     handles.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
        %     handles.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
        %     handles.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
        %     handles.IPCMRA = getappdata(0,'IPCMRA');
        %     pause(0.05)
        %     break
        % end
        % load data offset error JSOTELO
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [a,b,c,d] = size(handles.MR_FFE_FH);
        handles.a = a;
        handles.b = b;
        handles.c = c;
        handles.d = d;
        MR_FFE_FH_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_FFE_AP_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_FFE_RL_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_PCA_FH_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_PCA_AP_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_PCA_RL_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        IPCMRA_n = zeros(handles.a+2,handles.b+2,handles.c+2);
        MR_FFE_FH_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE_FH;
        MR_FFE_AP_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE_AP;
        MR_FFE_RL_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE_RL;
        MR_PCA_FH_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_FH;
        MR_PCA_AP_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_AP;
        MR_PCA_RL_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_RL;
        IPCMRA_n(2:end-1,2:end-1,2:end-1)  = handles.IPCMRA;
        handles.MR_FFE_FH   = MR_FFE_FH_n;
        handles.MR_FFE_AP   = MR_FFE_AP_n;
        handles.MR_FFE_RL   = MR_FFE_RL_n;
        handles.MR_PCA_FH   = MR_PCA_FH_n;
        handles.MR_PCA_AP   = MR_PCA_AP_n;
        handles.MR_PCA_RL   = MR_PCA_RL_n;
        handles.IPCMRA   = IPCMRA_n;
        [a,b,c,d] = size(handles.MR_FFE_FH);
        handles.a = a;
        handles.b = b;
        handles.c = c;
        handles.d = d;
%         handles.IPCMRA = (1/d)*sum( (handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);
        handles.IPCMRA = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % READING PAR-REC FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO: READ PAR-REC

    % TODO: Read DCM (Philips)

    % TODO: Read DAT_VD
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    val = X*handles.voxel_MR(1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2)
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    % axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    % axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    % axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step = zeros(2);
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    % set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    % set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    % set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    % set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    % set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    % set(handles.slider3,'visible','on')
    % set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    % set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    % set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
    
    % set(handles.pushbutton2, 'Units', 'pixels');
    % handles.pushbutton2_size = get(handles.pushbutton2, 'Position');
    % set(handles.pushbutton2, 'Units', 'normalized');
    % idx = mod(min(handles.pushbutton2_size(3:4)),2)>1;
    % w = floor(min(handles.pushbutton2_size(3:4)));
    % w(idx) = w(idx)+1;
    % im = imread('Symbols/P3.png');
    % g = double(imresize(im,[w-4 w-4],'method','nearest')>0);
    % set(handles.pushbutton2,'CData',g,'visible','on')
    % set(handles.pushbutton3,'CData',g,'visible','on')
    % set(handles.pushbutton4,'CData',g,'visible','on')
    % set(handles.pushbutton5,'CData',g)
    % set(handles.uipanel1,'Visible','on')
    % set(handles.uipanel6,'Visible','on')
    % set(handles.pushbutton1,'visible','on','BackgroundColor',[0 0 0])
    % set(handles.pushbutton15,'visible','on','BackgroundColor',[0.2 0.2 0.2])
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    % TODO: bwlabeln
    % [L,NUM] = bwlabeln(handles.SEG,6);
    % handles.L = L;
    % handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
        
    handles.id_mesh = 0;
    handles.id_unwrappping = 0;
    handles.id_vel = 0;
    handles.id_wss = 0;
    handles.id_osi = 0;
    handles.id_vor = 0;
    handles.id_hd   = 0;
    handles.id_rhd  = 0;
    handles.id_vd   = 0;
    handles.id_el   = 0;
    handles.id_ke   = 0;
    handles.id_lap  = 0;
    handles.id_cen  = 0;
    handles.id_rad  = 0;
    handles.id_dia  = 0;
    handles.id_auv  = 0;
    handles.id_cuv  = 0;
    handles.id_wssa = 0;
    handles.id_wssc = 0;
    handles.id_aan  = 0;
    handles.id_fve  = 0;
    handles.id_bve  = 0;
    handles.id_ref  = 0;
    handles.id_cenf = 0;
    handles.id_ecc  = 0;
    handles.id_cur  = 0; % Julio Sotelo 28-05-2019
    handles.id_ell  = 0; % Julio Sotelo 28-05-2019
    handles.id_len  = 0; % Julio Sotelo 28-05-2019
    handles.id_cir  = 0; % Julio Sotelo 28-05-2019
    handles.id_fov  = 0; % Julio Sotelo 28-05-2019
    handles.id_fla  = 0; % Julio Sotelo 28-05-2019
    handles.id_are  = 0; % Julio Sotelo 28-05-2019
    handles.id_aci  = 0; % Julio Sotelo 28-05-2019
    
    % UNUSED
    % handles.Lrgb_vel      = [];
    % handles.Lrgb_fve      = [];
    % handles.Lrgb_bve      = [];
    % handles.Lrgb_aan      = [];
    % handles.Lrgb_vor      = [];
    % handles.Lrgb_hd       = [];
    % handles.Lrgb_rhd      = [];
    % handles.Lrgb_vd       = [];
    % handles.Lrgb_el       = [];
    % handles.Lrgb_ke       = [];
    % handles.Lrgb_fov       = [];% Julio Sotelo 28-05-2019
        
    % handles.time            = [];
    % handles.peak_flow       = [];
    % handles.peak_flow_ori   = [];
    % handles.flow            = [];
    % handles.net_flow        = [];
    % handles.max_velocity    = [];
    % handles.min_velocity    = [];
    % 
    % handles.veset = [];
    % handles.VOR = [];
    % handles.WSS = [];
    % handles.mags_vel = [];
    % handles.mags_wss = [];
    % handles.mags_osi = [];
    % handles.mags_vor = [];
    % handles.mags_hd = [];
    % handles.mags_rhd = [];
    % handles.mags_vd = [];
    % handles.mags_el = [];
    % handles.mags_ke = [];
    % 
    % handles.Laplace = [];                         
    % handles.centerline = [];       
    % handles.centerline_lapid = [];
    % handles.radius = [];                        
    % handles.diameter = [];                        
    % handles.axial_unit_vectors = [];           
    % handles.circumferential_unit_vectors = [];   
    % handles.WSS_A = [];                         
    % handles.WSS_C = [];                          
    % handles.mag_WSS_A = [];                      
    % handles.mag_WSS_C = [];                     
    % handles.angle_axial_direction = [];           
    % handles.forward_velocity = [];               
    % handles.backward_velocity = [];               
    % handles.mag_forward_velocity = [];           
    % handles.mag_backward_velocity = [];           
    % handles.regurgitant_flow = [];               
    % handles.centerline_flow = [];                
    % handles.eccentricity = [];
    % 
    % handles.curvature  = []; % Julio Sotelo 28-05-2019
    % handles.ellipticity  = []; % Julio Sotelo 28-05-2019
    % handles.length_vessel  = []; % Julio Sotelo 28-05-2019
    % handles.circulation  = []; % Julio Sotelo 28-05-2019
    % handles.forward_vortex  = []; % Julio Sotelo 28-05-2019
    % handles.flattening  = []; % Julio Sotelo 28-05-2019
    % handles.area  = []; % Julio Sotelo 28-05-2019
    % handles.axial_circulation  = []; % Julio Sotelo 28-05-2019
    
    handles.save_id_SEG_mat = 0;
    handles.save_id_SEG_vti = 0;
    handles.save_id_IPCMRA_mat = 1;
    handles.save_id_IPCMRA_vti = 1;
    handles.save_id_MR_PCA_mat = 1;
    handles.save_id_MR_PCA_vti = 1;
    handles.save_id_MR_FFE_mat = 1;
    handles.save_id_MR_FFE_vti = 1;
    handles.save_id_mesh_mat = 0;
    handles.save_id_vel_mat = 0;
    handles.save_id_wss_mat = 0;
    handles.save_id_osi_mat = 0;
    handles.save_id_vor_mat = 0;
    handles.save_id_hd_mat = 0;
    handles.save_id_rhd_mat = 0;
    handles.save_id_vd_mat = 0;
    handles.save_id_el_mat = 0;
    handles.save_id_ke_mat = 0;
    handles.save_id_mesh_vtu = 0;
    handles.save_id_vel_vtu = 0;
    handles.save_id_wss_vtu = 0;
    handles.save_id_osi_vtu = 0;
    handles.save_id_vor_vtu = 0;
    handles.save_id_hd_vtu = 0;
    handles.save_id_rhd_vtu = 0;
    handles.save_id_vd_vtu = 0;
    handles.save_id_el_vtu = 0;
    handles.save_id_ke_vtu = 0;
    
    handles.save_id_lap_mat     = 0;
    handles.save_id_cen_mat     = 0;
    handles.save_id_rad_mat     = 0;
    handles.save_id_dia_mat     = 0;
    handles.save_id_auv_mat     = 0;
    handles.save_id_cuv_mat     = 0;
    handles.save_id_wssa_mat    = 0;
    handles.save_id_wssc_mat    = 0;
    handles.save_id_aan_mat     = 0;
    handles.save_id_fve_mat     = 0;
    handles.save_id_bve_mat     = 0;
    handles.save_id_ref_mat     = 0;
    handles.save_id_cebf_mat    = 0;
    handles.save_id_ecc_mat     = 0;

    handles.save_id_cur_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_ell_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_len_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_cir_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_fov_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_fla_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_are_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_aci_mat     = 0;% Julio Sotelo 28-05-2019
    
    handles.save_id_tim_mat     = 0; % Julio Sotelo 28-05-2019 time
    handles.save_id_flo_mat     = 0; % Julio Sotelo 28-05-2019 flow
    handles.save_id_nfl_mat     = 0; % Julio Sotelo 28-05-2019 net_flow
    handles.save_id_mav_mat     = 0; % Julio Sotelo 28-05-2019 max_velocity
    handles.save_id_miv_mat     = 0; % Julio Sotelo 28-05-2019 min_velocity
    
    
    handles.save_id_lap_vtu     = 0;
    handles.save_id_cen_vtu     = 0;
    handles.save_id_rad_vtu     = 0;
    handles.save_id_dia_vtu     = 0;
    handles.save_id_auv_vtu     = 0;
    handles.save_id_cuv_vtu     = 0;
    handles.save_id_wssa_vtu    = 0;
    handles.save_id_wssc_vtu    = 0;
    handles.save_id_aan_vtu     = 0;
    handles.save_id_fve_vtu     = 0;
    handles.save_id_bve_vtu     = 0;
    handles.save_id_ref_vtu     = 0;
    handles.save_id_cebf_vtu    = 0;
    handles.save_id_ecc_vtu     = 0;
    
    handles.save_id_cur_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_ell_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_len_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_cir_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_fov_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_fla_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_are_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_aci_vtu     = 0;% Julio Sotelo 28-05-2019
    
    handles.save_id_vel_csv   = 0;
    handles.save_id_wss_csv   = 0;
    handles.save_id_osi_csv   = 0;
    handles.save_id_vor_csv   = 0;
    handles.save_id_hd_csv    = 0;
    handles.save_id_rhd_csv   = 0;
    handles.save_id_vd_csv    = 0;
    handles.save_id_el_csv    = 0;
    handles.save_id_ke_csv    = 0;
    handles.save_id_rad_csv   = 0;
    handles.save_id_dia_csv   = 0;
    handles.save_id_wssa_csv  = 0;
    handles.save_id_wssc_csv  = 0;
    handles.save_id_aan_csv   = 0;
    handles.save_id_fve_csv   = 0;
    handles.save_id_bve_csv   = 0;
    handles.save_id_ref_csv   = 0;
    handles.save_id_ecc_csv   = 0;
    
    handles.save_id_cur_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_ell_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_len_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_cir_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_fov_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_fla_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_are_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_aci_csv     = 0;% Julio Sotelo 28-05-2019
    
    handles.save_id_tim_csv     = 0; % Julio Sotelo 28-05-2019 time
    handles.save_id_flo_csv     = 0; % Julio Sotelo 28-05-2019 flow
    handles.save_id_nfl_csv     = 0; % Julio Sotelo 28-05-2019 net_flow
    handles.save_id_mav_csv     = 0; % Julio Sotelo 28-05-2019 max_velocity
    handles.save_id_miv_csv     = 0; % Julio Sotelo 28-05-2019 min_velocity
    
%     UNUSED
%     handles.faces = [];
%     handles.elem = [];
%     handles.nodes = [];
    
    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
%     UNUSED
%     handles.position_sag_cor = [];
%     handles.position_axi_cor = [];
%     handles.position_cor_cor = [];
    
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [r,c,v] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    handles.id_ang = 1;
    handles.id_mag = 0;
    
    handles.id_csv_save = 0;
%     UNUSED
%     handles.SECTIONS_S = [];
%     handles.SECTIONS_V = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Showing the logo image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % txt = 'Symbols/VIEWS.png';
    % imlogo = imread(txt);
    % [av,bv,~]=size(imlogo);
    % windows_screen_size = get(0,'ScreenSize');
    % imlogo = imresize(imlogo,[round(windows_screen_size(4)*av/(av+bv)) round(windows_screen_size(4)*bv/(av+bv))]);
    % Sz= size(imlogo);
    % flogo = figure('Position',[10000 10000 Sz(2) + 4 Sz(1) + 4],'name','VIEWS','numbertitle','off','menubar','none');
    % movegui(flogo,'center');
    % set(flogo,'Units', 'pixels');
    % image(imlogo(:,:,1:3))
    % set(gca,'Visible','off','Units','pixels','Position', [2 2 Sz(2) Sz(1)]);
    % waitfor(flogo)
        
    % Julio Sotelo 18092022
    handles.exec_fe_mesh = 0;

% handles.output = hObject
% guidata(hObject, handles);