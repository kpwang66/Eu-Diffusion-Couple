% Reads raw neutron image data, crops, and returns it in a matrix

function cropped_data = Eu_read(filepath, option)

    Eu_data = readmatrix(filepath);
    Eu_data(isnan(Eu_data)) = 0;

    figure
    imagesc(Eu_data)

    upper_right = [1, size(Eu_data, 2)];
    bottom_left = [379, 1];
    
    cropped_data = Eu_data(1:bottom_left(1), bottom_left(2):upper_right(2));

    figure
    imagesc(cropped_data);

    axis equal
    colormap(flipud(gray));

    colorbar

    cropped_size = size(cropped_data);

    side_A = cropped_data(:, round(cropped_size(2)/2) + 1:end);
    side_B = fliplr(cropped_data(:, 1:round(cropped_size(2)/2)));

    % Average both sides of data
    merged_data = (side_A + side_B)./2;

    figure
    imagesc(merged_data);
    axis equal
    colormap(flipud(gray));
    colorbar

    if option=='A'
        cropped_data = side_A;
    elseif option == 'B'
        cropped_data = side_B;
    elseif option == 'merged'
        cropped_data = merged_data;
    end

%     keyboard
end