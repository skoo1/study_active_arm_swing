clear all;close all;clc;
%% read data
fileID = fopen('1_upper_limb_active.json','r');
mytxt = fscanf(fileID,'%s');
fclose(fileID);
mystruct1 = jsondecode(mytxt);
fldname = fieldnames(mystruct1);

fileID = fopen('2_upper_limb_fixed.json','r');
mytxt = fscanf(fileID,'%s');
fclose(fileID);
mystruct2 = jsondecode(mytxt);

fileID = fopen('3_upper_limb_passive.json','r');
mytxt = fscanf(fileID,'%s');
fclose(fileID);
mystruct3 = jsondecode(mytxt);
%% GRF detection
grf_temp = mystruct1.GRF_r_z;
grf_temp(:,2) = [diff(grf_temp);0];
idx1 = (grf_temp(:,1) == 0);
idx2 = (grf_temp(:,2) > 0);
hs_idx_temp = find((idx1+idx2) == 2);
j = 1;
for i = 1:length(hs_idx_temp)
    if hs_idx_temp(j) < 1
        hs_idx_temp(j) = [];
        continue;
    end
    if sum(grf_temp(hs_idx_temp(j)-1:(hs_idx_temp(j)), 1)) ~= 0
        hs_idx_temp(j) = [];
        j = j-1;
    end
    j = j+1;
end
hs_idx = hs_idx_temp;
n_stance = length(hs_idx) - 1;

figure(901);
plot(grf_temp(:,1))
y_lim = ylim;
for i = 1:n_stance
    line([hs_idx_temp(i) hs_idx_temp(i)], [y_lim(1) y_lim(2)], 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
end

hs_idx1 = hs_idx_temp;
n_stance1 = n_stance;
%% GRF - 1
name = sprintf('GRF_r_z');
name_x = sprintf('GRF_r_x');
name_y = sprintf('GRF_r_y');

% lowpass(mystruct1.(name), 15, 100)

i = 1;
for j = 1:n_stance1
    frame_temp = hs_idx(j):hs_idx(j+1);
    data_grf1_temp = spline(frame_temp, mystruct1.(name)(frame_temp), linspace(frame_temp(1), frame_temp(end), 101));
    data_grf1_x_temp = spline(frame_temp, mystruct1.(name_x)(frame_temp), linspace(frame_temp(1), frame_temp(end), 101));
    data_grf1_y_temp = spline(frame_temp, mystruct1.(name_y)(frame_temp), linspace(frame_temp(1), frame_temp(end), 101));
    
    x_diff = mystruct1.global_x(frame_temp(end)) - mystruct1.global_x(frame_temp(1));
    y_diff = mystruct1.global_y(frame_temp(end)) - mystruct1.global_y(frame_temp(1));
    myangle = atan2(y_diff, x_diff);
    data_grf1_x_temp_rot = data_grf1_x_temp*cos(myangle) + data_grf1_y_temp*sin(myangle);
    data_grf1_y_temp_rot = data_grf1_y_temp*cos(myangle) - data_grf1_x_temp*sin(myangle);
    
    %     if find(data_grf1_temp > 3)
    %         n_stance1 = n_stance1 - 1;
    %         continue;
    %     end
    frame1(i, 1:2) = [hs_idx(j), hs_idx(j+1)];
    data_grf1(:,i) = data_grf1_temp;
    data_grf1_x(:,i) = data_grf1_x_temp_rot;
    data_grf1_y(:,i) = data_grf1_y_temp_rot;
    i = i+1;
end

data_grf1(find(data_grf1 < 0)) = 0;

figure(101);hold on;
% for istance = 1:length(data_grf1(1,:))
%     plot(0:100, data_grf1(:,istance)*1000, 'Color', [0.75 0.75 0.75], 'LineWidth', 1)
%     plot(0:100, data_grf1_x(:,istance)*1000, 'Color', [1.0 0.5 0.5], 'LineWidth', 1)
%     plot(0:100, data_grf1_y(:,istance)*1000, 'Color', [0.5 0.5 1.0], 'LineWidth', 1)
% end
plot(0:100, mean(data_grf1(:,:)*1000, 2), 'Color', [0 0 0], 'LineWidth', 4)
plot(0:100, mean(data_grf1_x(:,:)*1000, 2), 'Color', [0.75 0 0], 'LineWidth', 4)
plot(0:100, mean(data_grf1_y(:,:)*1000, 2), 'Color', [0 0 0.75], 'LineWidth', 4)
axis([0 40 -500 2000])
set(gca, 'FontSize', 20)
%% GRF detection - 2
grf_temp = mystruct2.GRF_r_z;
grf_temp(:,2) = [diff(grf_temp);0];
idx1 = (grf_temp(:,1) == 0);
idx2 = (grf_temp(:,2) > 0);
hs_idx_temp = find((idx1+idx2) == 2);
j = 1;
for i = 1:length(hs_idx_temp)
    if hs_idx_temp(j) < 10
        hs_idx_temp(j) = [];
        continue;
    end
    if sum(grf_temp(hs_idx_temp(j)-10:(hs_idx_temp(j)), 1)) ~= 0
        hs_idx_temp(j) = [];
        j = j-1;
    end
    j = j+1;
end
hs_idx = hs_idx_temp;
n_stance = length(hs_idx) - 1;

figure(902);
plot(grf_temp(:,1))
y_lim = ylim;
for i = 1:n_stance
    line([hs_idx_temp(i) hs_idx_temp(i)], [y_lim(1) y_lim(2)], 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
end

hs_idx2 = hs_idx_temp;
n_stance2 = n_stance;
%% GRF - 2
name = sprintf('GRF_r_z');
name_x = sprintf('GRF_r_x');
name_y = sprintf('GRF_r_y');

% lowpass(mystruct1.(name), 15, 100)

i = 1;
for j = 1:n_stance2
    frame_temp = hs_idx(j):hs_idx(j+1);
    data_grf2_temp = spline(frame_temp, mystruct2.(name)(frame_temp), linspace(frame_temp(1), frame_temp(end), 101));
    data_grf2_x_temp = spline(frame_temp, mystruct2.(name_x)(frame_temp), linspace(frame_temp(1), frame_temp(end), 101));
    data_grf2_y_temp = spline(frame_temp, mystruct2.(name_y)(frame_temp), linspace(frame_temp(1), frame_temp(end), 101));
    %     if find(data_grf2_temp > 3)
    %         n_stance2 = n_stance2 - 1;
    %         continue;
    %     end
    
    x_diff = mystruct2.global_x(frame_temp(end)) - mystruct2.global_x(frame_temp(1));
    y_diff = mystruct2.global_y(frame_temp(end)) - mystruct2.global_y(frame_temp(1));
    myangle = atan2(y_diff, x_diff);
    data_grf2_x_temp_rot = data_grf2_x_temp*cos(myangle) + data_grf2_y_temp*sin(myangle);
    data_grf2_y_temp_rot = data_grf2_y_temp*cos(myangle) - data_grf2_x_temp*sin(myangle);
    
    frame2(i, 1:2) = [hs_idx(j), hs_idx(j+1)];
    data_grf2(:,i) = data_grf2_temp;
    data_grf2_x(:,i) = data_grf2_x_temp_rot;
    data_grf2_y(:,i) = data_grf2_y_temp_rot;
    i = i+1;
end

data_grf2(find(data_grf2 < 0)) = 0;

figure(102);hold on;
% for istance = 1:length(data_grf2(1,:))
%     plot(0:100, data_grf2(:,istance)*1000, 'Color', [0.75 0.75 0.75], 'LineWidth', 1)
%     plot(0:100, data_grf2_x(:,istance)*1000, 'Color', [1.0 0.5 0.5], 'LineWidth', 1)
%     plot(0:100, data_grf2_y(:,istance)*1000, 'Color', [0.5 0.5 1.0], 'LineWidth', 1)
% end
plot(0:100, mean(data_grf2(:,:)*1000, 2), 'Color', [0 0 0], 'LineWidth', 4)
plot(0:100, mean(data_grf2_x(:,:)*1000, 2), 'Color', [0.75 0 0], 'LineWidth', 4)
plot(0:100, mean(data_grf2_y(:,:)*1000, 2), 'Color', [0 0 0.75], 'LineWidth', 4)
axis([0 40 -500 2000])
set(gca, 'FontSize', 20)
%% GRF detection - 3
grf_temp = mystruct3.GRF_r_z;
grf_temp(:,2) = [diff(grf_temp);0];
idx1 = (grf_temp(:,1) == 0);
idx2 = (grf_temp(:,2) > 0);
hs_idx_temp = find((idx1+idx2) == 2);
j = 1;
for i = 1:length(hs_idx_temp)
    if hs_idx_temp(j) < 10
        hs_idx_temp(j) = [];
        continue;
    end
    if sum(grf_temp(hs_idx_temp(j)-9:(hs_idx_temp(j)), 1)) ~= 0
        hs_idx_temp(j) = [];
        j = j-1;
    end
    j = j+1;
end
hs_idx = hs_idx_temp;
n_stance = length(hs_idx) - 1;

figure(903);
plot(grf_temp(:,1))
y_lim = ylim;
for i = 1:n_stance
    line([hs_idx_temp(i) hs_idx_temp(i)], [y_lim(1) y_lim(2)], 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
end

hs_idx3 = hs_idx_temp;
n_stance3 = n_stance;
%% GRF - 3
name = sprintf('GRF_r_z');
name_x = sprintf('GRF_r_x');
name_y = sprintf('GRF_r_y');

% lowpass(mystruct1.(name), 15, 100)

i = 1;
for j = 1:n_stance3
    frame_temp = hs_idx(j):hs_idx(j+1);
    data_grf3_temp = spline(frame_temp, mystruct3.(name)(frame_temp), linspace(frame_temp(1), frame_temp(end), 101));
    data_grf3_x_temp = spline(frame_temp, mystruct3.(name_x)(frame_temp), linspace(frame_temp(1), frame_temp(end), 101));
    data_grf3_y_temp = spline(frame_temp, mystruct3.(name_y)(frame_temp), linspace(frame_temp(1), frame_temp(end), 101));
    %     if find(data_grf3_temp > 3)
    %         n_stance3 = n_stance3 - 1;
    %         continue;
    %     end
    
    x_diff = mystruct3.global_x(frame_temp(end)) - mystruct3.global_x(frame_temp(1));
    y_diff = mystruct3.global_y(frame_temp(end)) - mystruct3.global_y(frame_temp(1));
    myangle = atan2(y_diff, x_diff);
    data_grf3_x_temp_rot = data_grf3_x_temp*cos(myangle) + data_grf3_y_temp*sin(myangle);
    data_grf3_y_temp_rot = data_grf3_y_temp*cos(myangle) - data_grf3_x_temp*sin(myangle);
    
    frame3(i, 1:2) = [hs_idx(j), hs_idx(j+1)];
    data_grf3(:,i) = data_grf3_temp;
    data_grf3_x(:,i) = data_grf3_x_temp_rot;
    data_grf3_y(:,i) = data_grf3_y_temp_rot;
    i = i+1;
end

data_grf3(find(data_grf3 < 0)) = 0;

figure(103);hold on;
% for istance = 1:length(data_grf3(1,:))
%     plot(0:100, data_grf3(:,istance)*1000, 'Color', [0.75 0.75 0.75], 'LineWidth', 1)
%     plot(0:100, data_grf3_x(:,istance)*1000, 'Color', [1.0 0.5 0.5], 'LineWidth', 1)
%     plot(0:100, data_grf3_y(:,istance)*1000, 'Color', [0.5 0.5 1.0], 'LineWidth', 1)
% end
plot(0:100, mean(data_grf3(:,:)*1000, 2), 'Color', [0 0 0], 'LineWidth', 4)
plot(0:100, mean(data_grf3_x(:,:)*1000, 2), 'Color', [0.75 0 0], 'LineWidth', 4)
plot(0:100, mean(data_grf3_y(:,:)*1000, 2), 'Color', [0 0 0.75], 'LineWidth', 4)
axis([0 40 -500 2000])
set(gca, 'FontSize', 20)
close all;
%% upper limb kineamtics
for i = 4:7
    name = sprintf('gc%d',i-1);
    
    for j = 1:n_stance1
        frame_temp = frame1(j,1):frame1(j,2);
        data1(:,j,i) = spline(frame_temp, mystruct1.(name)(frame_temp), linspace(frame_temp(1), frame_temp(end), 101));
    end
    if i == 4
        data2(1:101,1:35,i) = 0;
    end
    if i == 5
        data2(1:101,1:35,i) = 35*pi/180;
    end
    if i == 6
        data2(1:101,1:35,i) = 0;
    end
    if i == 7
        data2(1:101,1:35,i) = 115*pi/180;
    end
    for j = 1:n_stance3
        frame_temp = frame3(j,1):frame3(j,2);
        data3(:,j,i) = spline(frame_temp, mystruct3.(name)(frame_temp), linspace(frame_temp(1), frame_temp(end), 101));
    end
    
    figure(200+i);hold on;
   
    plot(mean(data1(:,:,i)*180/pi, 2), 'Color', [0 0 0], 'LineWidth', 4)
    plot(mean(data2(:,:,i)*180/pi, 2), 'Color', [0.75 0 0], 'LineWidth', 4)
    plot(mean(data3(:,:,i)*180/pi, 2), 'Color', [0 0 0.75], 'LineWidth', 4)
    set(gca, 'FontSize', 20)
end
