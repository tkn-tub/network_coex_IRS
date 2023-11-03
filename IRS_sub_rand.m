clear all
close all
clc

%%%% Define number of operators and subsurfaces
U=5;
n_subs=5;
%%%% Define num_drops for different placements and trials for the rand
%%%% combinations generator
n_drops=20;
trials=20;

fc = 28e9;
c = physconst('lightspeed');
lambda = c/fc;

% Each network operates at a slighly different frequency
f=zeros(1,n_subs);
f= fc+1e8*(1:n_subs);

% Setup surface
Nr = 50;
Nc = 50;
dr = 0.5*lambda;
dc = 0.5*lambda;

% construct subsurfaces
for i=1:n_subs
ris{i} = helperRISSurface('Size',[Nr Nc],'ElementSpacing',[dr dc],...
    'ReflectorElement',phased.IsotropicAntennaElement,'OperatingFrequency',f(i));
end


%%%% For each drop, new locations
for l=1:n_drops
    v = zeros(3,1);
    pos_ap = zeros(U,3);
    for u = 1:U
        while true
            % Generate random angles for spherical coordinates (phi and theta)
            phi = rand(1) * 2 * pi;     % Azimuthal angle (0 to 2*pi)
            theta = rand(1) * pi;       % Polar angle (0 to pi)

            % Generate a random radius between 25 and 100
            r = 25 + rand(1) * 100;

            % Convert spherical coordinates to Cartesian coordinates
            x = r * sin(theta) * cos(phi);
            y = r * sin(theta) * sin(phi);
            z = r * cos(theta);
    
            % Check if the point is at least 25 units away from the origin
            if norm([x, y, z]) >= 25
                pos_ap(u, :) = [x, y, z];
                break;
            end
        end
    end


    pos_ue = zeros(U, 3);

    for u = 1:U
        while true
            % Generate random angles for spherical coordinates (phi and theta)
            phi = rand(1) * 2 * pi;     % Azimuthal angle (0 to 2*pi)
            theta = rand(1) * pi;       % Polar angle (0 to pi)

            % Generate a random radius between 25 and 100
            r = 25 + rand(1) * 100;

            % Convert spherical coordinates to Cartesian coordinates
            x = r * sin(theta) * cos(phi);
            y = r * sin(theta) * sin(phi);
            z = r * cos(theta);

            % Check if the point is at least 25 units away from the corresponding base station
            if norm([x, y, z] - pos_ap(u, :)) >= 25
                pos_ue(u, :) = [x, y, z];
                break;
            end
        end
    end


    dbr =  Nc*dc*(0:n_subs-1);
    % The coordinates of the subsurfaces
    pos_ris = [dbr;zeros(2,U)]; %%% Assuming there is a large IRS, subsurfaced
    
    pos_ap=pos_ap';
    pos_ue=pos_ue';



    % compute the range and angle of the RIS from the base station and the UE
    % Initialize arrays to store the results
    r_ap_ris = zeros(U, U);
    ang_ap_ris = cell(U, U);
    r_ue_ris = zeros(U, U);
    ang_ue_ris = cell(U, U);

    for i = 1:U
        for j = 1:n_subs
            % Calculate range and angle for pos_ap and pos_ris pairs
            [r_ap_ris(i, j), ang_ap_ris{i, j}] = rangeangle(pos_ap(:, i), pos_ris(:, j));

            % Calculate range and angle for pos_ue and pos_ris pairs
            [r_ue_ris(i, j), ang_ue_ris{i, j}] = rangeangle(pos_ue(:, i), pos_ris(:, j));
        end
    end

    % signal
    fs = 10e6;
    x = 2*randi(2,[100 1])-3;
    tx = phased.Transmitter('PeakPower',50e-3,'Gain',0);
    xt = tx(x);
    N0dB = -60-30;

    % channels --- LOS between AP and UE is ignored
    chanAPToRIS = phased.FreeSpace('SampleRate',fs,'PropagationSpeed',c,'MaximumDistanceSource','Property','MaximumDistance',500);
    chanRISToUE = phased.FreeSpace('SampleRate',fs,'PropagationSpeed',c,'MaximumDistanceSource','Property','MaximumDistance',500);

    rand_assignments_length=trials;

    RIS_comb_assignments_rand=zeros(rand_assignments_length,U);
    for g=1:rand_assignments_length
        RIS_comb_assignments_rand(g,:)=[randperm(U)];
    end

    SNRriso = zeros(rand_assignments_length, U);
    SNRriso_e = zeros(rand_assignments_length, U);
    SNR_f = zeros(rand_assignments_length, U);
    rcoeff_ris=cell(U, U);

    for m=1:rand_assignments_length
        csu=RIS_comb_assignments_rand(m,:);
        for i=1:U
            yriso=zeros(100,1);
            yriso_e=zeros(100,1);
            for j = 1:n_subs

                operator=i;

                % Extract the positions for the current pair of AP and UE
                current_pos_ap = pos_ap(:, operator);
                current_pos_ue = pos_ue(:, operator);
                if csu(j) ~= operator

                    stv = getSteeringVector(ris{j});

                    % Compute optimal phase control for the current pair

                    g = db2mag(-fspl(r_ap_ris(csu(j), j), lambda)) * exp(1i * 2 * pi * r_ap_ris(csu(j), j) / c) * stv(f(j), ang_ap_ris{csu(j), j});
                    hr = db2mag(-fspl(r_ue_ris(csu(j), j), lambda)) * exp(1i * 2 * pi * r_ue_ris(csu(j), j) / c) * stv(f(j), ang_ue_ris{csu(j),j});
                    rcoeff_ris{csu(j),j} = exp(1i * (-angle(hr) - angle(g)));

                    x_ris_in = chanAPToRIS(xt, current_pos_ap, pos_ris(:, j), v, v);

                    % Calculate x_ris_out for the current pair using the i-th RIS object
                    x_ris_out = ris{j}(x_ris_in, ang_ap_ris{operator, j}, ang_ue_ris{operator,j}, rcoeff_ris{csu(j),j});

                    % Calculate yriso for the current pair
                    yriso_e = yriso_e+chanRISToUE(x_ris_out, pos_ris(:, j), current_pos_ue, v, v);

                else
                    stv = getSteeringVector(ris{j});
                    g = db2mag(-fspl(r_ap_ris(operator, j), lambda)) * exp(1i * 2 * pi * r_ap_ris(operator, j) / c) * stv(f(j), ang_ap_ris{operator, j});
                    hr = db2mag(-fspl(r_ue_ris(operator, j), lambda)) * exp(1i * 2 * pi * r_ue_ris(operator, j) / c) * stv(f(j), ang_ue_ris{operator,j});
                    rcoeff_ris{operator,j} = exp(1i * (-angle(hr) - angle(g)));
                    x_ris_in = chanAPToRIS(xt, current_pos_ap, pos_ris(:, j), v, v);

                    % Calculate x_ris_out for the current pair using the i-th RIS object
                    x_ris_out = ris{j}(x_ris_in, ang_ap_ris{operator, j}, ang_ue_ris{operator,j}, rcoeff_ris{operator,j});

                    % Calculate yriso for the current pair
                    yriso = yriso+chanRISToUE(x_ris_out, pos_ris(:, j), current_pos_ue, v, v);

                end

            end
            SNRriso(m,i) = pow2db(bandpower(yriso)) - N0dB;
            SNRriso_e(m,i) = pow2db(bandpower(yriso_e)) - N0dB;
            SNR_f(m,i)=SNRriso_e(m,i)+SNRriso(m,i);
        end


    end


    rate=log2(1+db2pow(SNR_f));


    for k=1:length(RIS_comb_assignments_rand(:,1))
        sum_rate(l,k)=sum((rate(k,:)));
    end
end



hFig = figure;
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 101 700 400])
h3=cdfplot(min(sum_rate)');
set(h3,'LineWidth',2,'LineStyle','-.');
hold on;
h4=cdfplot(mean(sum_rate'));
set(h4,'LineWidth',2,'LineStyle','--');
hold on;
h5=cdfplot(max(sum_rate'));
set(h5,'LineWidth',2,'LineStyle','-','Color','black');
xlabel('\itSum Rate','FontSize', 15);
ylabel('\itFraction of different placements','FontSize', 15);
legend('Minimum possible outcome', 'Average possible outcome', 'Maximum possible outcome', 'FontAngle','italic', 'FontSize',13);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Serif','FontSize',15)
hold off;
