% Generate data from cell count macro ; JEC 1-22-14

clearvars 

plots = false; % plots 2 graphs if true

[files,pathname] = uigetfile('*.csv', 'Select CSV files', ...
   'MultiSelect', 'on');
                filenames_combine = [files(:)'];
                
x_i = inputdlg('Enter identifying text', 'Filename identifier', [1 100]);
    data_id = x_i{:};

% setup variable/field names for structured arrays
vars_names = genvarname(files);
% vars_names_tmp1 = strrep(files, '-', '_');
% vars_names_tmp2 = strrep(vars_names_tmp1, ' ', '_');
% vars_names_csv = strrep(vars_names_tmp2, '.tif', ''); % if?
% vars_names = strcat({'id_'},strrep(vars_names_csv,'CELLS.csv',''));

% put all the CSV data and x-y coordinates into structured arrays ('csv',
% 'ref_coords', and 'coords') and save to one MAT file
for i=1:numel(filenames_combine)
%     fld = genvarname(vars_names{i});
    fld = vars_names{i};
    source_file = char(filenames_combine(1,i));
    csv.(fld) = importdata([pathname source_file]);               
end 

for i=1:numel(filenames_combine) 
    fld = vars_names{i};
    
    ref_coords.x.(fld) = csv.(fld).data(1,4)';
    ref_coords.y.(fld) = csv.(fld).data(1,5)';
    
    coords.x.(fld) = csv.(fld).data(2:end,4); % skip the first row, containsthe center ref (area = 900)
    coords.y.(fld) = csv.(fld).data(2:end,5); % skip the first row, contains the center ref (area = 900)
    x_y=[coords.x coords.y];
end
save([pathname data_id '_cellScatterData.mat'],'files','csv','vars_names','ref_coords','coords'); 



% get x-y coordinates for each image, calculate distance matrix values, and
% save to MAT files
for i=1:numel(filenames_combine)
    x = x_y(1,1).(vars_names{i});
    y = x_y(1,2).(vars_names{i});

    array_xy=[x y];
    fld = vars_names{i};
    
    if isempty(array_xy) == 1
        array_xy = NaN;
    end
    
    %calculate distance matrix (distmat) and Delaunay
    [dmat,opt] = distmat(array_xy);
        dmat_all_distances=dmat(2:end,1); % for saving to individual MAT files
        dmat_all_distances_csv.(fld)=dmat(2:end,1); % for saving to one CSV file
            dist_matrix_mean=mean(dmat_all_distances);
            dist_matrix_std=std(dmat_all_distances);
    if ~isnan(array_xy) == 1
        dt = DelaunayTri(array_xy);
        save([pathname vars_names{i}],'array_xy','dmat_all_distances','dist_matrix_mean','dist_matrix_std','dt');
    else
        save([pathname vars_names{i}],'array_xy','dmat_all_distances','dist_matrix_mean','dist_matrix_std');
    end
    % dn = delaunayn(x_y);
    

end   

struct2csv(dmat_all_distances_csv,[data_id '_dist_matrix.csv']);

% plot data
if plots == true
    
    figure(1); % Scatter plot of all distances
        plot(ones(length(dmat_all_distances)),dmat_all_distances,'co'); hold on;
    %             plot(1,dist_matrix_mean,'ks', 'MarkerSize', 8);
            errorbar(dist_matrix_mean,dist_matrix_std/sqrt(length(dmat_all_distances)),'gs','LineWidth',0.5);
        boxplot(dmat_all_distances);

    figure(2);
        triplot(dt,'m-'); hold on;
        plot(array_xy(:,1),array_xy(:,2),'bo');
    %         plot(x_y(4:6,1),x_y(4:6,2),'ko');

%     figure(3)
%         boxplot(dmat_all_distances);
  
    hold off;
%     pause % wait for user input (press any key to continue)
else
end



% *************************IGNORE BELOW - EXTRA STUF FROM PREVIOUS****************************************               
%     b=mean(dmat(:,1));
%     bb=std(dmat(:,1));
%     dmat_b=dmat(:,1);
%     figure(1)
%     boxplot(dmat_all_distances);
%     figure(2)
%     boxplot(dmat_b(:,1));

%     hist(dmat(:,1),20)
    
%     h = findobj(gca,'Type','patch');
%     set(h,'FaceColor','b')
 
    %%%%%%%%%%





% 
% figure(2);
%     hist(dt(:,1));
% figure(3);
%     hist(dt(:,2));
% figure(4);
%     hist(dt(:,3));
    

    
%     figure(1);
%     [theta,rho] = cart2pol(xcoord,ycoord);
%         polar(theta,rho,'b-');
%         hold all


% rgb=imread(rgb);
% imshow(rgb,'initialmagnification',100);
% hold on
% M = size(rgb,1);
% N = size(rgb,2);

% Plot cells
%scatter(polex,poley,'ro'); % sets origin
%figure(2);
% for n = 1:length(xcoord)
% %     plot(roix1(n),roiy1(n)+polex,'yo','MarkerSize',6);
% %     plot(xcoord(n),ycoord(n),'ys','MarkerFace','y','MarkerSize',10);
% % plot(xcoord(n)+polex,ycoord(n)+poley,'ws','MarkerFace','g','MarkerSize',8);
% plot(xcoord(:),ycoord(:),'bs','MarkerSize',8);
% end
% hold all

% % 
% % % % Generate tif and eps files in current directory
% % % saveas(gcf, [filesave], 'tif');
% % % print('-depsc',filesave)
% % % 
% % % disp(M);
% % % disp(N);
% % % hold off


