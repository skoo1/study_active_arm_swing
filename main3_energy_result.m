clear all;close all;clc;
%% read data
fileID = fopen('1_energy_active.json','r');
mytxt = fscanf(fileID,'%s');
fclose(fileID);
mystruct1 = jsondecode(mytxt);
fldname = fieldnames(mystruct1);

fileID = fopen('2_energy_fixed.json','r');
mytxt = fscanf(fileID,'%s');
fclose(fileID);
mystruct2 = jsondecode(mytxt);

fileID = fopen('3_energy_passive.json','r');
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
%%
distnace1 = 0;
distnace2 = 0;
distnace3 = 0;
for itime = frame1(5,1):frame1(35,2)
    x_diff1 = (mystruct1.global_x(itime+1) - mystruct1.global_x(itime));
    y_diff1 = (mystruct1.global_y(itime+1) - mystruct1.global_y(itime));
    distnace1 = distnace1 + sqrt(x_diff1^2 + y_diff1^2);
end

for itime = frame2(5,1):frame2(35,2)
    x_diff2 = (mystruct2.global_x(itime+1) - mystruct2.global_x(itime));
    y_diff2 = (mystruct2.global_y(itime+1) - mystruct2.global_y(itime));
    distnace2 = distnace2 + sqrt(x_diff2^2 + y_diff2^2);
end

for itime = frame3(5,1):frame3(35,2)
    x_diff3 = (mystruct3.global_x(itime+1) - mystruct3.global_x(itime));
    y_diff3 = (mystruct3.global_y(itime+1) - mystruct3.global_y(itime));
    distnace3 = distnace3 + sqrt(x_diff3^2 + y_diff3^2);
end

vel1 = distnace1/(frame1(35,2) - frame1(5,1))/0.01;
vel2 = distnace2/(frame2(35,2) - frame2(5,1))/0.01;
vel3 = distnace3/(frame3(35,2) - frame3(5,1))/0.01;
close all;
%% activation
run('ReadMuslceName.m');
Name_lower = Name;
clear Name
run('ReadMuslceName_upper_limb.m');
Name_upper = Name;
activation1_lower = zeros(frame1(35,2)-frame1(5,1)+1, 1);
activation2_lower = zeros(frame2(35,2)-frame2(5,1)+1, 1);
activation3_lower = zeros(frame3(35,2)-frame3(5,1)+1, 1);

activation1_torso = zeros(frame1(35,2)-frame1(5,1)+1, 1);
activation2_torso = zeros(frame2(35,2)-frame2(5,1)+1, 1);
activation3_torso = zeros(frame3(35,2)-frame3(5,1)+1, 1);

activation1_upper = zeros(frame1(35,2)-frame1(5,1)+1, 1);
activation2_upper = zeros(frame2(35,2)-frame2(5,1)+1, 1);
activation3_upper = zeros(frame3(35,2)-frame3(5,1)+1, 1);

c1 = 0.25;
c2 = -1.2;
for i = 1:86
    name1 = [Name_lower{i}, '_velocity'];
    name2 = [Name_lower{i}, '_force'];
    
    frame1_positive = find(mystruct1.(name1)(frame1(5,1):frame1(35,2)) < 0);
    frame2_positive = find(mystruct1.(name1)(frame2(5,1):frame2(35,2)) < 0);
    frame3_positive = find(mystruct1.(name1)(frame3(5,1):frame3(35,2)) < 0);
    
    frame1_negative = find(mystruct1.(name1)(frame1(5,1):frame1(35,2)) > 0);
    frame2_negative = find(mystruct1.(name1)(frame2(5,1):frame2(35,2)) > 0);
    frame3_negative = find(mystruct1.(name1)(frame3(5,1):frame3(35,2)) > 0);
    
    activation1_lower(frame1_positive) = activation1_lower(frame1_positive) + abs(mystruct1.(name1)(frame1_positive).*mystruct1.(name2)(frame1_positive))/c1*0.01/72*0.769;
    activation2_lower(frame2_positive) = activation2_lower(frame2_positive) + abs(mystruct2.(name1)(frame2_positive).*mystruct2.(name2)(frame2_positive))/c1*0.01/72*0.769;
    activation3_lower(frame3_positive) = activation3_lower(frame3_positive) + abs(mystruct3.(name1)(frame3_positive).*mystruct3.(name2)(frame3_positive))/c1*0.01/72;
    
    activation1_lower(frame1_negative) = activation1_lower(frame1_negative) + abs(mystruct1.(name1)(frame1_negative).*mystruct1.(name2)(frame1_negative))/c2*0.01/72*0.769;
    activation2_lower(frame2_negative) = activation2_lower(frame2_negative) + abs(mystruct2.(name1)(frame2_negative).*mystruct2.(name2)(frame2_negative))/c2*0.01/72*0.769;
    activation3_lower(frame3_negative) = activation3_lower(frame3_negative) + abs(mystruct3.(name1)(frame3_negative).*mystruct3.(name2)(frame3_negative))/c2*0.01/72;
end

iter = 1;
for j = 5:35
    frame_temp1 = (frame1(j,1) - frame1(5,1) + 1):(frame1(j,2) - frame1(5,1) + 1);
    frame_temp2 = (frame2(j,1) - frame2(5,1) + 1):(frame2(j,2) - frame2(5,1) + 1);
    frame_temp3 = (frame3(j,1) - frame3(5,1) + 1):(frame3(j,2) - frame3(5,1) + 1);
    lower1(:,iter) = spline(frame_temp1, activation1_lower(frame_temp1), linspace(frame_temp1(1), frame_temp1(end), 101));
    lower2(:,iter) = spline(frame_temp2, activation2_lower(frame_temp2), linspace(frame_temp2(1), frame_temp2(end), 101));
    lower3(:,iter) = spline(frame_temp3, activation3_lower(frame_temp3), linspace(frame_temp3(1), frame_temp3(end), 101));
    iter = iter+1;
end

for i = 87:92
    name1 = [Name_lower{i}, '_velocity'];
    name2 = [Name_lower{i}, '_force'];
    
    frame1_positive = find(mystruct1.(name1)(frame1(5,1):frame1(35,2)) < 0);
    frame2_positive = find(mystruct1.(name1)(frame2(5,1):frame2(35,2)) < 0);
    frame3_positive = find(mystruct1.(name1)(frame3(5,1):frame3(35,2)) < 0);
    
    frame1_negative = find(mystruct1.(name1)(frame1(5,1):frame1(35,2)) > 0);
    frame2_negative = find(mystruct1.(name1)(frame2(5,1):frame2(35,2)) > 0);
    frame3_negative = find(mystruct1.(name1)(frame3(5,1):frame3(35,2)) > 0);
    
    activation1_torso(frame1_positive) = activation1_torso(frame1_positive) + abs(mystruct1.(name1)(frame1_positive).*mystruct1.(name2)(frame1_positive))/c1*0.01/72*0.769;
    activation2_torso(frame2_positive) = activation2_torso(frame2_positive) + abs(mystruct2.(name1)(frame2_positive).*mystruct2.(name2)(frame2_positive))/c1*0.01/72*0.769;
    activation3_torso(frame3_positive) = activation3_torso(frame3_positive) + abs(mystruct3.(name1)(frame3_positive).*mystruct3.(name2)(frame3_positive))/c1*0.01/72*0.769;
    
    activation1_torso(frame1_negative) = activation1_torso(frame1_negative) + abs(mystruct1.(name1)(frame1_negative).*mystruct1.(name2)(frame1_negative))/c2*0.01/72*0.769;
    activation2_torso(frame2_negative) = activation2_torso(frame2_negative) + abs(mystruct2.(name1)(frame2_negative).*mystruct2.(name2)(frame2_negative))/c2*0.01/72*0.769;
    activation3_torso(frame3_negative) = activation3_torso(frame3_negative) + abs(mystruct3.(name1)(frame3_negative).*mystruct3.(name2)(frame3_negative))/c2*0.01/72*0.769;
end

iter = 1;
for j = 5:35
    frame_temp1 = (frame1(j,1) - frame1(5,1) + 1):(frame1(j,2) - frame1(5,1) + 1);
    frame_temp2 = (frame2(j,1) - frame2(5,1) + 1):(frame2(j,2) - frame2(5,1) + 1);
    frame_temp3 = (frame3(j,1) - frame3(5,1) + 1):(frame3(j,2) - frame3(5,1) + 1);
    torso1(:,iter) = spline(frame_temp1, activation1_torso(frame_temp1), linspace(frame_temp1(1), frame_temp1(end), 101));
    torso2(:,iter) = spline(frame_temp2, activation2_torso(frame_temp2), linspace(frame_temp2(1), frame_temp2(end), 101));
    torso3(:,iter) = spline(frame_temp3, activation3_torso(frame_temp3), linspace(frame_temp3(1), frame_temp3(end), 101));
    iter = iter+1;
end

for i = 93:150
    name1 = [Name_upper{i-92}, '_velocity'];
    name2 = [Name_upper{i-92}, '_force'];
    
    frame1_positive = find(mystruct1.(name1)(frame1(5,1):frame1(35,2)) < 0);
    frame2_positive = find(mystruct1.(name1)(frame2(5,1):frame2(35,2)) < 0);
    frame3_positive = find(mystruct1.(name1)(frame3(5,1):frame3(35,2)) < 0);
    
    frame1_negative = find(mystruct1.(name1)(frame1(5,1):frame1(35,2)) > 0);
    frame2_negative = find(mystruct1.(name1)(frame2(5,1):frame2(35,2)) > 0);
    frame3_negative = find(mystruct1.(name1)(frame3(5,1):frame3(35,2)) > 0);
    
    activation1_upper(frame1_positive) = activation1_upper(frame1_positive) + abs(mystruct1.(name1)(frame1_positive).*mystruct1.(name2)(frame1_positive))/c1*0.01/72*0.769;
    activation2_upper(frame2_positive) = activation2_upper(frame2_positive) + 0*0.769;
    activation3_upper(frame3_positive) = activation3_upper(frame3_positive) + 0*0.769;
    
    activation1_upper(frame1_negative) = activation1_upper(frame1_negative) + abs(mystruct1.(name1)(frame1_negative).*mystruct1.(name2)(frame1_negative))/c2*0.01/72*0.769;
    activation2_upper(frame2_negative) = activation2_upper(frame2_negative) + 0*0.769;
    activation3_upper(frame3_negative) = activation3_upper(frame3_negative) + 0*0.769;
end

iter = 1;
for j = 5:35
    frame_temp1 = (frame1(j,1) - frame1(5,1) + 1):(frame1(j,2) - frame1(5,1) + 1);
    frame_temp2 = (frame2(j,1) - frame2(5,1) + 1):(frame2(j,2) - frame2(5,1) + 1);
    frame_temp3 = (frame3(j,1) - frame3(5,1) + 1):(frame3(j,2) - frame3(5,1) + 1);
    upper1(:,iter) = spline(frame_temp1, activation1_upper(frame_temp1), linspace(frame_temp1(1), frame_temp1(end), 101));
    upper2(:,iter) = spline(frame_temp2, activation2_upper(frame_temp2), linspace(frame_temp2(1), frame_temp2(end), 101));
    upper3(:,iter) = spline(frame_temp3, activation3_upper(frame_temp3), linspace(frame_temp3(1), frame_temp3(end), 101));
    iter = iter+1;
end

basal1 = ones(length(activation1_lower), 1)*1.2*0.01;
basal2 = ones(length(activation2_lower), 1)*1.2*0.01;
basal3 = ones(length(activation3_lower), 1)*1.2*0.01;

lower1 = (sum(activation1_lower(:)))/distnace1;
lower2 = (sum(activation2_lower(:)))/distnace2;
lower3 = (sum(activation3_lower(:)))/distnace3;

torso1 = (sum(activation1_torso(:)))/distnace1;
torso2 = (sum(activation2_torso(:)))/distnace2;
torso3 = (sum(activation3_torso(:)))/distnace3;

upper1 = (sum(activation1_upper(:)))/distnace1;
upper2 = (sum(activation2_upper(:)))/distnace2;
upper3 = (sum(activation3_upper(:)))/distnace3;

basal1 = sum(basal1)/distnace1;
basal2 = sum(basal2)/distnace2;
basal3 = sum(basal3)/distnace3;

result = [upper1, upper2, upper3;
    torso1, torso2, torso3;
    lower1, lower2, lower3;
    basal1, basal2, basal3]