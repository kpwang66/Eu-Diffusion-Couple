function Eu_1D_prism_data = Eu_get1D_cyl(z_slice_data, rad_data, n_pix)

size_z_slice = size(z_slice_data);
size_rad = size(rad_data);
Eu_1D_prism_data = zeros(1, size_rad(1));

x = -size_z_slice(2)/2 +.5 :size_z_slice(2)/2 -.5;
y = -size_z_slice(1)/2 + .5 :size_z_slice(1)/2 -.5;
[xx yy] = meshgrid(x,y);
tot_circ_pix = sum(sum((xx.^2+yy.^2) < n_pix^2))/2;

% u((xx.^2+yy.^2)<n_pix^2)=0;   % radius 100, center at the origin

for row = 1:size_rad(1)
   
    % Take one layer of z-slice
    z_slice_row_dat = z_slice_data(:,:, row);

% Test plot    
%     if row == round(size_rad(1)/2)
%         figure
%         imagesc(z_slice_row_dat, [0 12])
%         colormap(gray)
%         colorbar
%         axis equal
%         drawnow;
% %         clims([0 12])
%     end
    
    % All the pixels that fall in the range of n_pix are 0.  (Hollow
    % cylinder)
    tot_circ_pix = sum(sum((xx.^2+yy.^2) < n_pix^2));
    z_slice_row_dat((xx.^2+yy.^2) < n_pix^2)=0;

% Test plot   
%     imagesc(z_slice_row_dat)
%     if row == round(size_rad(1)/2)
%         figure
%         imagesc(z_slice_row_dat, [0 12])
%         colormap(gray)
%         colorbar
%         axis equal
% %         clims([0 12])
%         drawnow;
% %         keyboard
%     end

    % Crop  
    z_slice_row_dat_cropped = z_slice_row_dat(:, round(size_z_slice(2)/2):round(size_z_slice(2)/2) + n_pix - 1);
    
    
%     z_slice_rad1 = radon(z_slice_row_dat_cropped, 90);
%     
%     figure
%     plot(z_slice_rad1)
%     title('z_slice_rad1')
    
    z_slice_rad2 = sum(z_slice_row_dat_cropped);

%     title('z_slice_rad2')
    
    rad_data_cropped = rad_data(row, 1:n_pix);
    
%     if row == round(size_rad(1)/2)
%         figure
%         plot(rad_data_cropped, 'LineWidth', 2)
%         hold all;
%         plot(z_slice_rad2, 'LineWidth', 2)
%         plot(rad_data_cropped - z_slice_rad2, 'LineWidth', 2)
% %         keyboard
%     end
    
    Eu_1D_prism_data(row) = sum(rad_data_cropped - z_slice_rad2)/tot_circ_pix;
%     Eu_1D_prism_data(row) = mean(rad_data_cropped - z_slice_rad2);
	
end