clear all;
close all;
clc

%%%% Define number of operators, subsurfaces and  n_drops for different placements 
U=5;
n_drops=100;
n_subs=5;

fc = 28e9;
c = physconst('lightspeed');
lambda = c/fc;

f= fc+1e8*(1:n_subs);

% Setup surface
Nr = 50;
Nc = 50;
dr = 0.5*lambda;
dc = 0.5*lambda;

% construct surfaces
for i=1:n_subs
    ris{i} = helperRISSurface('Size',[Nr Nc],'ElementSpacing',[dr dc],...
        'ReflectorElement',phased.IsotropicAntennaElement,'OperatingFrequency',f(i));
end

%%%% Each drop corresponds to new setup
for l=1:num_drops

    v = zeros(3,1);
    pos_ap = zeros(U,3);
    for u = 1:U
        while true
            % Generate random angles for spherical coordinates (phi and theta)
            phi = rand(1) * 2 * pi;     % Azimuthal angle (0 to 2*pi)
            theta = rand(1) * pi;       % Polar angle (0 to pi)

            % Generate a random radius between 25 and 100
            r = 25 + rand(1) * 75;

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
            r = 25 + rand(1) * 75;

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
    pos_ris = [dbr;zeros(2,U)];

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

    % channel
    chanAPToRIS = phased.FreeSpace('SampleRate',fs,'PropagationSpeed',c,'MaximumDistanceSource','Property','MaximumDistance',500);
    chanRISToUE = phased.FreeSpace('SampleRate',fs,'PropagationSpeed',c,'MaximumDistanceSource','Property','MaximumDistance',500);

    % LOS path propagation
    num_pairs = size(pos_ap, 2);  % Get the number of pairs (assuming pos_ap and pos_ue have the same number of columns)

    % Initialize arrays to store the results
    yref = cell(1, num_pairs);
    SNRref = zeros(1, num_pairs);


    SNRriso = zeros(1, U);
    SNRriso_extra = zeros(U, U);
    rcoeff_ris=cell(U, U);
    used_indices = zeros(1, U);
    %csu used to hold the postion of each operator placed one at a time
    csu=zeros(1,U); 
    %%



    for i=1:U
        max_snr = -Inf; 
        max_j = 0;

        for j = 1:n_subs

            operator=i;

            % Extract the positions for the current pair of AP and UE
            current_pos_ap = pos_ap(:, operator);
            current_pos_ue = pos_ue(:, operator);
            yriso_e=zeros(100,1);
            yriso=zeros(100,1);
            if csu(j) == 0

                stv = getSteeringVector(ris{j});

                % Compute optimal phase control for the current pair

                g = db2mag(-fspl(r_ap_ris(operator, j), lambda)) * exp(1i * 2 * pi * r_ap_ris(operator, j) / c) * stv(f(j), ang_ap_ris{operator, j});
                hr = db2mag(-fspl(r_ue_ris(operator, j), lambda)) * exp(1i * 2 * pi * r_ue_ris(operator, j) / c) * stv(f(j), ang_ue_ris{operator,j});
                rcoeff_ris{operator,j} = exp(1i * (-angle(hr) - angle(g)));

                x_ris_in = chanAPToRIS(xt, current_pos_ap, pos_ris(:, j), v, v);

                % Calculate x_ris_out for the current pair using the i-th RIS object
                x_ris_out = ris{j}(x_ris_in, ang_ap_ris{operator, j}, ang_ue_ris{operator,j}, rcoeff_ris{operator,j});

                % Calculate yriso for the current pair
                yriso = chanRISToUE(x_ris_out, pos_ris(:, j), current_pos_ue, v, v);

                for k=1:n_subs
                    if csu(k)~=0
                        stv = getSteeringVector(ris{k});

                        % Compute optimal phase control for the current pair

                        x_ris_in = chanAPToRIS(xt, current_pos_ap, pos_ris(:, k), v, v);

                        % Calculate x_ris_out for the current pair using the i-th RIS object
                        x_ris_out = ris{k}(x_ris_in, ang_ap_ris{operator, k}, ang_ue_ris{operator,k}, rcoeff_ris{csu(k),k});

                        % Calculate yriso for the current pair
                        yriso_e = yriso_e + chanRISToUE(x_ris_out, pos_ris(:, k), current_pos_ue, v, v);
                    end
                end
            end

            if sum(yriso_e)~=0
                SNRriso(operator, j) = pow2db(bandpower(yriso)) - N0dB+pow2db(bandpower(yriso_e))-N0dB;
            else
                SNRriso(operator, j) = pow2db(bandpower(yriso))-N0dB;
            end
            if SNRriso(operator, j) > max_snr && ~used_indices(j)
                max_snr = SNRriso(operator, j);
                max_j = j;  % Update max_j with the index that gives the maximum SNRriso for 'i'
                % In other words, choose the subsurface which leads to
                % least inter operator interference
            end


        end

        % Store the 'j' value corresponding to the maximum SNRriso for 'i' in the csu array
        csu(max_j) = operator;
        used_indices(max_j) = 1;

    end

    sum_rate(l)=compute_SR(csu,U,n_subs,pos_ap,pos_ue,ris,r_ap_ris,lambda,c,f,ang_ap_ris,chanAPToRIS,v,pos_ris,r_ue_ris,ang_ue_ris,rcoeff_ris,chanRISToUE,N0dB,xt);


end

      

