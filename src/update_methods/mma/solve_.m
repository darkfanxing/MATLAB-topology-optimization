% ==================================================================
%
%    Author: Krister Svanberg <krille@math.kth.se>
%            Department of Mathematics, KTH,
%            SE-10044 Stockholm, Sweden.
%
%    Public Date: Dec 2006.
%
% ==================================================================
%
% This function solves the MMA subproblem:
%         
%     minimize     sum_{ p_{0j} / (U_j - x_j) + q_{0j} / (x_j - L_j) }
%                  + sum_{ c_i * y_i + 0.5*d_i*(y_i)^2 }
%                  + a_0*z
%
%     subject to   sum_{ p_{ij} / (U_j - x_j) + q_{ij} / (x_j - L_j) } - a_i*z - y_i <= b_i,
%                  alpha_j <=  x_j <=  beta_j,
%                  y_i >= 0,
%                  z >= 0

function [ ...
    x, ...
    y, ...
    z, ...
    general_mma_constraint_lagrange_multipliers, ...
    alpha_lagrange_multipliers, ...
    beta_lagrange_multipliers, ...
    constraint_lagrange_multipliers, ...
    z_lagrange_multiplier, ...
    general_mma_constraint_slack ...
] = solve_(m, n, epsi_min, lower_bound, upper_bound, alpha, beta, p0, q0, p_i_vector, q_i_vector, a0, a_i, b, c, d)
    design_variable_unit_array = ones(n, 1);
    constraint_unit_array = ones(m, 1);
    
    epsi = 1;

    x = 0.5 * (alpha + beta);
    y = constraint_unit_array;
    z = 1;

    general_mma_constraint_lagrange_multipliers = constraint_unit_array;
    alpha_lagrange_multipliers = design_variable_unit_array ./ (x - alpha);
    alpha_lagrange_multipliers = max(alpha_lagrange_multipliers, design_variable_unit_array);
    beta_lagrange_multipliers = design_variable_unit_array ./ (beta - x);
    beta_lagrange_multipliers = max(beta_lagrange_multipliers, design_variable_unit_array);
    constraint_lagrange_multipliers  = max(constraint_unit_array, 0.5 * c);
    z_lagrange_multiplier = 1;
    general_mma_constraint_slack = constraint_unit_array;
    
    iteration_number = 0;
    while epsi > epsi_min
          epsvecn = epsi * design_variable_unit_array;
          epsvecm = epsi * constraint_unit_array;
    
          u_minus_x = upper_bound - x;
          x_minus_l = x - lower_bound;
          
          dl_dx = ...
              (p0 + p_i_vector' * general_mma_constraint_lagrange_multipliers) ./ (u_minus_x .* u_minus_x) ...
              - (q0 + q_i_vector' * general_mma_constraint_lagrange_multipliers) ./ (x_minus_l .* x_minus_l);

          df_dx = dl_dx - alpha_lagrange_multipliers + beta_lagrange_multipliers;
          df_dy = c + d .* y - constraint_lagrange_multipliers - general_mma_constraint_lagrange_multipliers;
          df_dz = a0 - z_lagrange_multiplier - a_i' * general_mma_constraint_lagrange_multipliers;
    
          df_dlambda = ...
              p_i_vector * (design_variable_unit_array ./ u_minus_x) ...
              + q_i_vector * (design_variable_unit_array ./ x_minus_l) ...
              - a_i * z ...
              - y ...
              + general_mma_constraint_slack ...
              - b;

          rexsi = alpha_lagrange_multipliers .* (x - alpha) - epsvecn;

          reeta = beta_lagrange_multipliers .* (beta - x) - epsvecn;
          remu = constraint_lagrange_multipliers .* y - epsvecm;
          rezet = z_lagrange_multiplier*z - epsi;
          res = general_mma_constraint_lagrange_multipliers.*general_mma_constraint_slack - epsvecm;
          residu1 = [df_dx' df_dy' df_dz]';
          residu2 = [df_dlambda' rexsi' reeta' remu' rezet res']';
          residu = [residu1' residu2']';
          residunorm = sqrt(residu'*residu);
          residumax = max(abs(residu));

          ittt = 0;
          while residumax > 0.9*epsi && ittt < 200
            ittt = ittt + 1;
            iteration_number = iteration_number + 1;
            u_minus_x = upper_bound-x;
            x_minus_l = x-lower_bound;
            u_minus_x_square = u_minus_x.*u_minus_x;
            x_minus_l_square = x_minus_l.*x_minus_l;
            ux3 = u_minus_x.*u_minus_x_square;
            xl3 = x_minus_l.*x_minus_l_square;
            uxinv1 = design_variable_unit_array./u_minus_x;
            xlinv1 = design_variable_unit_array./x_minus_l;
            uxinv2 = design_variable_unit_array./u_minus_x_square;
            xlinv2 = design_variable_unit_array./x_minus_l_square;
            plam = p0 + p_i_vector'*general_mma_constraint_lagrange_multipliers ;
            qlam = q0 + q_i_vector'*general_mma_constraint_lagrange_multipliers ;
            gvec = p_i_vector*uxinv1 + q_i_vector*xlinv1;
            GG = p_i_vector*spdiags(uxinv2,0,n,n) - q_i_vector*spdiags(xlinv2,0,n,n);
            dpsidx = plam./u_minus_x_square - qlam./x_minus_l_square ;
            delx = dpsidx - epsvecn./(x-alpha) + epsvecn./(beta-x);
            dely = c + d.*y - general_mma_constraint_lagrange_multipliers - epsvecm./y;
            delz = a0 - a_i'*general_mma_constraint_lagrange_multipliers - epsi/z;
            dellam = gvec - a_i*z - y - b + epsvecm./general_mma_constraint_lagrange_multipliers;
            diagx = plam./ux3 + qlam./xl3;
            diagx = 2*diagx + alpha_lagrange_multipliers./(x-alpha) + beta_lagrange_multipliers./(beta-x);
            diagxinv = design_variable_unit_array./diagx;
            diagy = d + constraint_lagrange_multipliers./y;
            diagyinv = constraint_unit_array./diagy;
            diaglam = general_mma_constraint_slack./general_mma_constraint_lagrange_multipliers;
            diaglamyi = diaglam+diagyinv;
            if m < n
              blam = dellam + dely./diagy - GG*(delx./diagx);
              bb = [blam' delz]';
              Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
              AA = [Alam     a_i
                    a_i'    -z_lagrange_multiplier/z ];
              solut = AA\bb;
              dlam = solut(1:m);
              dz = solut(m+1);
              dx = -delx./diagx - (GG'*dlam)./diagx;
            else
              diaglamyiinv = constraint_unit_array./diaglamyi;
              dellamyi = dellam + dely./diagy;
              Axx = spdiags(diagx,0,n,n) + GG'*spdiags(diaglamyiinv,0,m,m)*GG;
              azz = z_lagrange_multiplier/z + a_i'*(a_i./diaglamyi);
              axz = -GG'*(a_i./diaglamyi);
              bx = delx + GG'*(dellamyi./diaglamyi);
              bz  = delz - a_i'*(dellamyi./diaglamyi);
              AA = [Axx   axz
                    axz'  azz ];
              bb = [-bx' -bz]';
              solut = AA\bb;
              dx  = solut(1:n);
              dz = solut(n+1);
              dlam = (GG*dx)./diaglamyi - dz*(a_i./diaglamyi) + dellamyi./diaglamyi;
            end
        %
            dy = -dely./diagy + dlam./diagy;
            dxsi = -alpha_lagrange_multipliers + epsvecn./(x-alpha) - (alpha_lagrange_multipliers.*dx)./(x-alpha);
            deta = -beta_lagrange_multipliers + epsvecn./(beta-x) + (beta_lagrange_multipliers.*dx)./(beta-x);
            dmu  = -constraint_lagrange_multipliers + epsvecm./y - (constraint_lagrange_multipliers.*dy)./y;
            dzet = -z_lagrange_multiplier + epsi/z - z_lagrange_multiplier*dz/z;
            ds   = -general_mma_constraint_slack + epsvecm./general_mma_constraint_lagrange_multipliers - (general_mma_constraint_slack.*dlam)./general_mma_constraint_lagrange_multipliers;
            xx  = [ y'  z  general_mma_constraint_lagrange_multipliers'  alpha_lagrange_multipliers'  beta_lagrange_multipliers'  constraint_lagrange_multipliers'  z_lagrange_multiplier  general_mma_constraint_slack']';
            dxx = [dy' dz dlam' dxsi' deta' dmu' dzet ds']';
        %    
            stepxx = -1.01*dxx./xx;
            stmxx  = max(stepxx);
            stepalpha = -1.01*dx./(x-alpha);
            stmalpha = max(stepalpha);
            stepbeta = 1.01*dx./(beta-x);
            stmbeta = max(stepbeta);
            stmalbe  = max(stmalpha,stmbeta);
            stmalbexx = max(stmalbe,stmxx);
            stminv = max(stmalbexx,1);
            steg = 1/stminv;
        %
            xold   =   x;
            yold   =   y;
            zold   =   z;
            lamold =  general_mma_constraint_lagrange_multipliers;
            xsiold =  alpha_lagrange_multipliers;
            etaold =  beta_lagrange_multipliers;
            muold  =  constraint_lagrange_multipliers;
            zetold =  z_lagrange_multiplier;
            sold   =   general_mma_constraint_slack;
        %
            itto = 0;
            resinew = 2*residunorm;
            while resinew > residunorm & itto < 50
                itto = itto+1;
                x   =   xold + steg*dx;
                y   =   yold + steg*dy;
                z   =   zold + steg*dz;
                general_mma_constraint_lagrange_multipliers = lamold + steg*dlam;
                alpha_lagrange_multipliers = xsiold + steg*dxsi;
                beta_lagrange_multipliers = etaold + steg*deta;
                constraint_lagrange_multipliers  = muold  + steg*dmu;
                z_lagrange_multiplier = zetold + steg*dzet;
                general_mma_constraint_slack   =   sold + steg*ds;
        
                u_minus_x = upper_bound-x;
                x_minus_l = x-lower_bound;
                u_minus_x_square = u_minus_x.*u_minus_x;
                x_minus_l_square = x_minus_l.*x_minus_l;
        
                uxinv1 = design_variable_unit_array./u_minus_x;
                xlinv1 = design_variable_unit_array./x_minus_l;
        
                plam = p0 + p_i_vector'*general_mma_constraint_lagrange_multipliers ;
                qlam = q0 + q_i_vector'*general_mma_constraint_lagrange_multipliers ;
                gvec = p_i_vector*uxinv1 + q_i_vector*xlinv1;
                dpsidx = plam./u_minus_x_square - qlam./x_minus_l_square ;
                rex = dpsidx - alpha_lagrange_multipliers + beta_lagrange_multipliers;
                rey = c + d.*y - constraint_lagrange_multipliers - general_mma_constraint_lagrange_multipliers;
                rez = a0 - z_lagrange_multiplier - a_i'*general_mma_constraint_lagrange_multipliers;
                relam = gvec - a_i*z - y + general_mma_constraint_slack - b;
                rexsi = alpha_lagrange_multipliers.*(x-alpha) - epsvecn;
                reeta = beta_lagrange_multipliers.*(beta-x) - epsvecn;
                remu = constraint_lagrange_multipliers.*y - epsvecm;
                rezet = z_lagrange_multiplier*z - epsi;
                res = general_mma_constraint_lagrange_multipliers.*general_mma_constraint_slack - epsvecm;
                residu1 = [rex' rey' rez]';
                residu2 = [relam' rexsi' reeta' remu' rezet res']';
                residu = [residu1' residu2']';
                resinew = sqrt(residu'*residu);
                steg = steg/2;
            end

            residunorm = resinew;
            residumax = max(abs(residu));
            steg = 2 * steg;
          end

          epsi = 0.1*epsi;
    end
end