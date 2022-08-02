clear all;close all;fclose all;clc;

this_script = matlab.desktop.editor.getActiveFilename;
this_dir = split(this_script, '\');
this_dir = strjoin(this_dir(1:end-1),'\');
cd(this_dir);


data_dir = [this_dir '\Raw Neutron Data'];
Eu_data_old = [this_dir '\Eu_data.mat'];
addpath(this_dir, data_dir)

cd(data_dir);
ls_txt = ls('*.txt');
[num_files a] = size(ls_txt);

iradon_dataset = cell(1, num_files);
raw_dataset = cell(1, num_files);
Anton_Eu1d_dataset = cell(1, num_files);
Eu1d_dataset = cell(1, num_files);
avgC_dataset = cell(1, num_files);
z_slice_dataset = (cell(1, num_files));

max_C = 0;
min_C = Inf;

start_file = 1;
end_file = 12;


menu=[
0   % Read original data
0   % Process files from raw data
1   % Show 2D results
1   % Show results with raw data
1   % Produce 1D, averaged profiles
1   % Just print raw data, uncropped
];

ampoule_wall = [109 114]; % pixels where ampoule wall are

clims = [0 10];

num_angles = 90;
raw_data_side = 'merged'; % Read side 'A' (left), side 'B' (right), or 'merged' (average the two)

dir_name = ['IRad_export_a-' num2str(num_angles), '_', raw_data_side, '_cyl_b'];
save_file = [this_dir, '\', dir_name, '\', dir_name, '.mat'];
save_file_1D = [this_dir, '\',  dir_name, '\', dir_name, '_1D.mat'];

px_wd_1d = 60; % pixel width 1d: number of pixels from the center to average in order to get a pseudo-1D profile
px_pitch = 6e-3/106; % m
r_px = 106; % radius in pixels

% Important variables:
% z_slice_data: 268 (l) x 268 (w) x 379 (depth) 3D matrix 

%% Process each file
% Menu 1
if logical(menu(1)) 
    for k = start_file:end_file %num_files
        t = tic;
        current_file = [data_dir '\' ls_txt(k, :)];

        fprintf('File %i...', k);
        fprintf('cropping...');
        cropped_data = Eu_read(current_file, raw_data_side);

        avgC_dataset{k} = cropped_data;
        
        fprintf('sample thickness...')

        for col = 1:r_px
            sample_thickness = 2*sqrt(r_px^2 - (col-1)^2);
            cropped_data(:, col) = cropped_data(:, col)*sample_thickness;
        end
        
        raw_dataset{k} = cropped_data;
        if logical(menu(2)) 
            
            %% Menu 2
            fprintf('iradon...')
            [irad_data z_slice_data] = Eu_irad(cropped_data, 'L', num_angles);

            iradon_dataset{k}  = irad_data;
            z_slice_dataset{k} = z_slice_data;
            
        end
        
        s = toc(t);

        fprintf('DONE!  Elapsed time: %2.1fs\n\n', s);

    end

    % Save results
    if ~exist([this_dir dir_name], 'dir')
       mkdir([this_dir dir_name]);
    end
    save(save_file, 'iradon_dataset', 'raw_dataset', 'z_slice_dataset', 'avgC_dataset', '-v7.3');
end



saved_irad_data = load(save_file); 
V = VideoWriter('.\Vid_Process_iRadon.avi');
V.FrameRate = 2;
open(V);
%% Menu 3
if logical(menu(3)) 
    figure
    for k = start_file:end_file%num_files
        
%         f = figure('visible', 'off');
        title_txt = ['t = ', num2str(22*k-22), ' minutes'];
        
        % If not processing fresh data, just plot original data
        if ~logical(menu(2))
            sgtitle(title_txt);
            imagesc(saved_irad_data.raw_dataset{k})
            axis([0 size(saved_irad_data.raw_dataset{k}, 2) 0 size(saved_irad_data.raw_dataset{k}, 1)])
            axis equal

            set(gcf, 'Position', [ 100 100 300 600]);
            colormap(flipud(gray))
            
            xlim([0 150])
            ylim([0 375])
            
            set(gca,'Ytick',[]) 
            set(gca,'Xtick',[]) 
            set(gca,'Yticklabel',[]) 
            set(gca,'Xticklabel',[])
            
            frame = getframe(gcf);
            writeVideo(V,frame);
            
%             keyboard
        % Plot result next to original
        elseif logical(menu(4))  
%             sgtitle(title_txt);
            subplot(1, 2, 1)
            imagesc(saved_irad_data.raw_dataset{k})
            axis([0 size(saved_irad_data.iradon_dataset{k}, 2) 0 size(saved_irad_data.iradon_dataset{k}, 1)])
            axis equal
            colorbar
            subplot(1, 2, 2)
            imagesc(saved_irad_data.iradon_dataset{k})
            axis([0 size(saved_irad_data.iradon_dataset{k}, 2) 0 size(saved_irad_data.iradon_dataset{k}, 1)])
            axis equal
            colorbar
            caxis(clims)

            set(gcf, 'Position', [ 100 100 900 600]);
            colormap(flipud(gray))
            
        else % Only plot processed result
            imagesc(saved_irad_data.iradon_dataset{k})
    %         colormap('gray')
            title_txt = ['t = ', num2str(22*k), ' minutes'];
%             title(title_txt);
            axis equal;
            caxis(clims);
            axis([0 size(saved_irad_data.iradon_dataset{k}, 2) 0 size(saved_irad_data.iradon_dataset{k}, 1)]);
            set(gcf, 'Position', [ 100 100 300 600]);
            colorbar;
            title(title_txt);
            colormap(flipud(gray))
            drawnow;
        end
        
        
        datafilename = strsplit(ls_txt(k, :), '.');
        datafilename = datafilename{1};
        figfilename = [this_dir, '\', dir_name, '\',  datafilename '.png'];
        saveas(gcf,figfilename);
        fprintf('Figure saved as %s\n', figfilename);


%         close;
    end
    
    close(V);
end

if logical(menu(5)) 
    Eu1d_dataset = zeros(length(saved_irad_data.raw_dataset{12}), 12);
    Anton_Eu1d_dataset = zeros(length(saved_irad_data.raw_dataset{12}), 12);
    
    % Plot 1D profiles near center
    figure
    for k = start_file:end_file%num_files
        % Produce 1D
        
        Eu1d_dataset(:, k) = 100*fliplr(Eu_get1D_cyl(saved_irad_data.z_slice_dataset{k}, saved_irad_data.raw_dataset{k}, px_wd_1d));
        Anton_Eu1d_dataset(:, k) = 100*fliplr(mean(saved_irad_data.avgC_dataset{k}(:,1:px_wd_1d)'));

        subplot(3, 4, k);

        % inverse-radon processed 1D data set
        plot(Eu1d_dataset(:, k), 'b');hold on

        % dataset we got from Anton
        plot(Anton_Eu1d_dataset(:, k), 'r');
        title_txt = ['t = ', num2str(22*k-k), ' minutes'];
        title(title_txt);
        ylim([0 1000])

    end
    
    save(save_file_1D, 'Eu1d_dataset', 'Anton_Eu1d_dataset')
end

if logical(menu(6)) 
    V2 = VideoWriter('.\Vid_Process_rawuncut.avi');
    V2.FrameRate = 2;
    open(V2);
    for k = start_file:end_file
        
        title_txt = ['t = ', num2str(22*k-22), ' minutes'];

        
        current_file = [data_dir '\' ls_txt(k, :)];
        Eu_data_rawuncut = readmatrix(current_file);
        Eu_data_rawuncut(isnan(Eu_data_rawuncut)) = 0;
        figure
        imagesc(Eu_data_rawuncut)
        sgtitle(title_txt)
%         axis([0 2*size(saved_irad_data.raw_dataset{k}, 2) 0 size(saved_irad_data.raw_dataset{k}, 1)])
        axis equal
        set(gcf, 'Position', [ 100 100 300 500]);
        colormap(flipud(gray))

        xlim([0 280])
        ylim([0 550])

        set(gca,'Ytick',[]) 
        set(gca,'Xtick',[]) 
        set(gca,'Yticklabel',[]) 
        set(gca,'Xticklabel',[])

        frame = getframe(gcf);
        writeVideo(V2,frame);

    end
    
    close(V2);
end

