clear;clc;close all;

% % load ground truth depth
load('nyu_depth_v2_labeled.mat', 'depths')
load('nyu_depth_v2_labeled.mat', 'images')
load('nyu_depth_v2_labeled.mat', 'rawDepths')
load('list_Trans_str.mat')
load('list_Trans.mat')
depths = double(1000*depths);
filename = 'val_list.txt'
idx = readmatrix(filename);
normalize = false;

filepath = './result_test';
midname = '/result_';

filepath2 = './paper_results';

% % NLSPN
p2 = '/NLSPN_epoch0020';  

% %  HED data path
filepath3 = './HED_examples';
% %  HED file name header
midname2 = '/rgb_'; 
% % Load image

I_NLSPN = (readall(imageDatastore([filepath2  p2 '/*.png'])));



I_our_tep = ((imageDatastore([filepath  '/*.png'])));
I_our = {};
I_hed = {};
I_rgb = {};
I_depth = {};
our_missed= [];

idx = readmatrix(filename);

for i = 1:length(idx)
    if  isfile([filepath midname num2str(idx(i)) '.png']) 
        I_our{i} = double(imread([filepath midname num2str(idx(i)) '.png']));
        I_hed{i} = imread([filepath3 midname2 num2str(idx(i)) '.png']);
        I_rgb{i} = images(:,:,:,idx(i));
        I_depth{i} = double(rawDepths(:,:,idx(i)));
    else
         I_our{i}=[];
         I_hed{i}=[];
        NLSPN_NYU = [our_missed;i];
        I_depth{i} =[];
    end
   
end
% % gt collect
I_gt = {};
for o = 1:length(idx)
    I_gt{o} = depths(:,:,idx(o));
end

% % downsize
I_our_resize = {};
I_NYU_VAL = {};

I_hed_resize = {};
I_gt_resize = {};
I_depth_resize = {};

for s = 1:length(I_our)
    if length(I_our{s}) > 0
        I_temp = imresize(I_our{s},0.5);
        I_temp = I_temp(7:234, 9:312);
        % % normalize

        I_temp_hed = imresize(I_hed{s},0.5);
        I_temp_hed = I_temp_hed(7:234, 9:312,1:3);
        I_hed_resize{s} =I_temp_hed;  

        I_gt_temp =  imresize(I_gt{s},0.5);
        I_gt_temp = I_gt_temp(7:234, 9:312);

        I_depth_temp = imresize(I_depth{s},0.5);
        I_depth_temp = I_depth_temp(7:234, 9:312);

        % % normalize

        I_our_resize{s} = I_temp/max(max(I_temp));
        I_NYU_VAL_tmp = double(I_NLSPN{s})/max(max(double(I_NLSPN{s})));
        I_NYU_VAL{s} = I_NYU_VAL_tmp;


        I_gt_resize{s} = I_gt_temp;
        I_depth_resize{s} = I_depth_temp;
      
    else
        I_our_resize{s} = [];
        I_NYU_VAL{s} = [];

        I_hed_resize{s} = [];
        I_depth_resize{s} = [];
        I_gt_resize{s} = [];
    end
end



% % % 


bw_region = {};
RMSE_NYU = {};
RMSE_our = {};

RMSE_data = {};
MAE_NYU = {};
MAE_our = {};

MAE_data = {};
iRMSE_NYU = {};
iRMSE_our = {};

iRMSE_data = {};
iMAE_NYU = {};
iMAE_our = {};

iMAE_data = {};
Rel_NYU = {};
Rel_our = {};

Rel_data = {};
thre_NYU = {};
thre_our = {};

thre_data = {};
m_count = {};
idx_mix_gather = {};
I_gradient = {};
parfor g = 1:length(I_our_resize)

    I_NYU_test = I_NYU_VAL{g};

    I_our_test = I_our_resize{g};
    I_gt_test = I_gt_resize{g};
    I_input =  I_hed_resize{g};
    I_depth_input = I_depth_resize{g};
    if length(I_our_resize{g}) > 0
    % % load region

    I_for_level = rgb2gray(I_input);
    I_for_level = adapthisteq(I_for_level);
    level=graythresh(I_for_level);
    BW = imbinarize(I_for_level,level);
    [label_hed_surface n_hs]=bwlabel(BW,4);
    label_re=reshape(label_hed_surface,[],1);
    BWALL={};
    n_hs=1:n_hs;
        for hjj=1:length(n_hs)
            check_hs=zeros(size(label_hed_surface));
            check_hs(find(label_re==n_hs(hjj)))=1;
            check_hs = bwmorph(check_hs,'erode',5);
            BWALL{hjj}=check_hs;
        end
        
    bw_region {g}= BWALL;
    m_count_test = {};
    RMSE_NYU_test ={};
    RMSE_our_test = {};

    RMSE_data_test = {};

    MAE_NYU_test = {};

    MAE_our_test = {};
    MAE_data_test = {};

    iRMSE_NYU_test ={};
    iRMSE_our_test = {};
    
    iRMSE_data_test = {};    

    iMAE_NYU_test = {};

    iMAE_our_test = {};
    iMAE_data_test = {};

    Rel_NYU_test = {};
    Rel_our_test = {};

    Rel_data_test = {};

    thre_NYU_test = {};
    thre_our_test = {};

    thre_data_test = {};    
    idx_mix_tmp = {};
    I_gradient_test = {};
        for ss = 1:length(n_hs)
            I_NYU_test_input = I_NYU_test.*BWALL{ss};

            I_our_test_input = I_our_test .*BWALL{ss};
            
            I_gt_test_input = I_gt_test.*BWALL{ss};
            I_depth_test_input = I_depth_input.*BWALL{ss};
            
            I_NYU_test_input = I_NYU_test_input/max(max(I_NYU_test_input));

            I_our_test_input = I_our_test_input/max(max(I_our_test_input));
            I_gt_test_input = I_gt_test_input/max(max(I_gt_test_input));
            I_gradient_test {ss} = gradient(gradient(I_gt_test_input.*BWALL{ss}));

% %  available data
            idx_NYU  = find(I_NYU_test_input~=0);
            idx_our = find(I_our_test_input~=0);

            if length(idx_NYU)<length(idx_our) 
                idx_mix = idx_NYU;

            else
                idx_mix = idx_our;
            end
            idx_mix_tmp{ss} = idx_mix;
           ava_count = length(( (idx_mix)));


    % %         RMSE
            m = find(reshape(BWALL{ss},[],1)==1);
            m_count_test{ss} = length(m);
            RMSE_NYU_test{ss} =   ( sqrt(sum(sum(((I_NYU_test_input( (idx_mix))-I_gt_test_input( (idx_mix))).^2)))/ava_count));
            RMSE_our_test{ss} =   ( sqrt(sum(sum(((I_our_test_input( (idx_mix))-I_gt_test_input( (idx_mix))).^2)))/ava_count));
           
            RMSE_data_test{ss} =  ( sqrt(sum(sum(((I_depth_test_input( (idx_mix))-I_gt_test_input( (idx_mix))).^2)))/ava_count));
           
           
    % %     MAE
            MAE_NYU_test{ss} =   ( sum(sum(abs(I_NYU_test_input((idx_mix))-I_gt_test_input( (idx_mix)))))/ava_count) ;
            MAE_our_test{ss} =   ( sum(sum(abs(I_our_test_input( (idx_mix))-I_gt_test_input( (idx_mix)))))/ava_count);
         
            MAE_data_test{ss} =  ( sum(sum(abs(I_depth_test_input( (idx_mix))-I_gt_test_input( (idx_mix)))))/ava_count);
% % irmse

            iRMSE_NYU_test{ss} =    sqrt(1/ava_count*sum(sum((1./I_NYU_test_input( (idx_mix))-1./I_gt_test_input( (idx_mix))).^2)));
            iRMSE_our_test{ss} =    sqrt(1/ava_count*sum(sum((1./I_our_test_input( (idx_mix))-1./I_gt_test_input( (idx_mix))).^2)));
       
            iRMSE_data_test{ss} =   sqrt(1/ava_count*sum(sum((1./I_depth_test_input( (idx_mix))-1./I_gt_test_input( (idx_mix))).^2)));
% % iMAE
            iMAE_NYU_test{ss} =  1/ava_count*(sum(abs(1./I_NYU_test_input( (idx_mix))-1./I_gt_test_input( (idx_mix)))));
            iMAE_our_test{ss} =  1/ava_count*(sum(abs(1./I_our_test_input( (idx_mix))-1./I_gt_test_input( (idx_mix)))));
       
            iMAE_data_test{ss} =  1/ava_count*(sum(abs(1./I_depth_test_input( (idx_mix))-1./I_gt_test_input( (idx_mix)))));

    % %     Rel
            Rel_NYU_test{ss} = sum(abs(I_NYU_test_input( (idx_mix))-I_gt_test_input( (idx_mix)))./I_gt_test_input( (idx_mix)))/ava_count;
            Rel_our_test{ss} =  sum(abs(I_our_test_input( (idx_mix))-I_gt_test_input( (idx_mix)))./I_gt_test_input( (idx_mix)))/ava_count;
     
            Rel_data_test{ss} =  sum(abs(I_depth_test_input( (idx_mix))-I_gt_test_input( (idx_mix)))./I_gt_test_input( (idx_mix)))/ava_count;

% %  threshold accuracy
                thre_NYU_testp = max([I_NYU_test_input( (idx_mix))./ I_gt_test_input( (idx_mix)) I_gt_test_input( (idx_mix))./I_NYU_test_input( (idx_mix))]')';
                thre_NYU_test{ss} = [sum(thre_NYU_testp < 1.25)/length(thre_NYU_testp) sum(thre_NYU_testp < 1.25^2)/length(thre_NYU_testp)...
                    sum(thre_NYU_testp < 1.25^3)/length(thre_NYU_testp)];
               
                thre_our_testp = max([I_our_test_input( (idx_mix))./ I_gt_test_input( (idx_mix)) I_gt_test_input( (idx_mix))./I_our_test_input( (idx_mix))]')';
                thre_our_test{ss} =  [sum(thre_our_testp < 1.25)/length(thre_our_testp) sum(thre_our_testp < 1.25^2)/length(thre_our_testp)...
                    sum(thre_our_testp < 1.25^3)/length(thre_our_testp)];
               

               
                thre_data_testp = max([I_depth_test_input( (idx_mix))./ I_gt_test_input( (idx_mix)) I_gt_test_input( (idx_mix))./I_depth_test_input( (idx_mix))]')';
                thre_data_test{ss} = [sum(thre_data_testp < 1.25)/length(thre_data_testp) sum(thre_data_testp < 1.25^2)/length(thre_data_testp)...
                    sum(thre_data_testp < 1.25^3)/length(thre_data_testp)];

        end
    
    
    RMSE_NYU{g} = RMSE_NYU_test;
    RMSE_our{g} = RMSE_our_test;

    RMSE_data{g} = RMSE_data_test;

    m_count{g} = m_count_test;

    MAE_NYU{g} = MAE_NYU_test;
    MAE_our{g} = MAE_our_test;

    MAE_data{g} = MAE_data_test;

    iRMSE_NYU{g} = iRMSE_NYU_test;
    iRMSE_our{g} = iRMSE_our_test;

    iRMSE_data{g} = iRMSE_data_test;
        
    iMAE_NYU{g} = iMAE_NYU_test;
    iMAE_our{g} = iMAE_our_test;

    iMAE_data{g} = iMAE_data_test;

    Rel_NYU{g} = Rel_NYU_test;
    Rel_our{g} = Rel_our_test;

    Rel_data{g} = Rel_data_test;
    
    thre_NYU{g} = thre_NYU_test;
    thre_our{g} = thre_our_test;

    thre_data{g} = thre_data_test;

    idx_mix_gather{g} = idx_mix_tmp;
    else
        RMSE_NYU{g} ={};
        RMSE_our{g} = {};

        RMSE_data{g} = {};

        m_count{g} = {};

        MAE_NYU{g} = {};
        MAE_our{g} = {};
 
        MAE_data{g} = {};

        iRMSE_NYU{g} = {};
        iRMSE_our{g} = {};
 
        iRMSE_data{g} = {}; 

        iMAE_NYU{g} = {};
        iMAE_our{g} = {};

        iMAE_data{g} = {};

        Rel_NYU{g} = {};
        Rel_our{g} = {};

        Rel_data{g} = {};
        
        thre_NYU{g} = {};
        thre_our{g} = {};

        thre_data{g} = {};
        idx_mix_gather{g} = {};

    end 

    
    end

    
