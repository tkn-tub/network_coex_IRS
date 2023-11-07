 
function [sum_rate]=compute_SR(csu,U,n_subs,pos_ap,pos_ue,ris,r_ap_ris,lambda,c,f,ang_ap_ris,chanAPToRIS,v,pos_ris,r_ue_ris,ang_ue_ris,rcoeff_ris,chanRISToUE,N0dB,xt)

for i=1:U
        yriso1=zeros(100,1);
        yriso1_e=zeros(100,1);
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
                yriso1_e = yriso1_e+chanRISToUE(x_ris_out, pos_ris(:, j), current_pos_ue, v, v);


            else
                stv = getSteeringVector(ris{j});
                g = db2mag(-fspl(r_ap_ris(operator, j), lambda)) * exp(1i * 2 * pi * r_ap_ris(operator, j) / c) * stv(f(j), ang_ap_ris{operator, j});
                hr = db2mag(-fspl(r_ue_ris(operator, j), lambda)) * exp(1i * 2 * pi * r_ue_ris(operator, j) / c) * stv(f(j), ang_ue_ris{operator,j});
                rcoeff_ris{operator,j} = exp(1i * (-angle(hr) - angle(g)));
                x_ris_in = chanAPToRIS(xt, current_pos_ap, pos_ris(:, j), v, v);

                % Calculate x_ris_out for the current pair using the i-th RIS object
                x_ris_out = ris{j}(x_ris_in, ang_ap_ris{operator, j}, ang_ue_ris{operator,j}, rcoeff_ris{csu(j),j});

                % Calculate yriso for the current pair
                yriso1 = yriso1+chanRISToUE(x_ris_out, pos_ris(:, j), current_pos_ue, v, v);

            end
        end

        SNRriso1(i) = pow2db(bandpower(yriso1)) - N0dB;
        SNRriso_e(i) = pow2db(bandpower(yriso1_e)) - N0dB;
        SNR_f(i)=SNRriso_e(i)+SNRriso1(i);
    end

    rate=log2(1+db2pow(SNR_f));
    utility=sum(log(rate));
    sum_rate=sum((rate));
    fairness_index=((sum(rate))^2)/(U*(sum((rate).^2)));