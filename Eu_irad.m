% This function takes a energy-resolved neutron image, takes each row, and
% performs the inverse-radon transformation on it to recover a true 2D
% z-slice data.  

function [iradon_data z_slice_data]= Eu_irad(radon_data, direction, num_deg)

data_size = size(radon_data);

row_times = zeros(1,data_size(1));
z_slice_data = zeros(2*data_size(2), 2*data_size(2), data_size(1));

for row = 1:data_size(1)
    radon_row = zeros( ceil(data_size(2)*sqrt(2)), 1 );
    radon_row(1:length(radon_data(row, :)')) = radon_data(row, :)';
    radon_row = [flipud(radon_row);
            radon_row];
%     tic;
    % Inverse radon image
    IRI = iradon(repmat(radon_row, 1, num_deg), linspace(0, 179, num_deg),'spline');

%     keyboard
    
    [Zr z_slice_avg]= aziavg_spline(IRI);
   
    iradon_data(row, :) = Zr;
    z_slice_data(:,:, row) = z_slice_avg;

end

% keyboard
end