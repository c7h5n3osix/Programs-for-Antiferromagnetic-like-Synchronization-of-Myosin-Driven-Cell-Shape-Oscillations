function phase_diagram_D_vs_E_Final_Corrected
%% 0. 环境初始化
    tic;
    clear; 
    % close all; 
    clc;
    
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool; 
    end
    
    %% 1. 物理参数 
    p.k_off   = 1;      
    p.eta     = 0.05;   
    p.alpha   = 0.01;     
    p.sigma_L = 30;      
    p.K       = 50;     
    
    fixed_sigma = 14; 
    Cs = 1 / (1 + p.k_off); 
    
    param_str = sprintf('\\sigma_0=%g, \\eta=%g, \\alpha=%g, k_{off}=%g, K=%g', ...
                        fixed_sigma, p.eta, p.alpha, p.k_off, p.K);
    
    %% 2. 向量化准备
    resolution = 500; 
    D_vec = linspace(0.00, 3, resolution)'; 
    n_D = length(D_vec);
    
    k_val_vec = 0.00:0.05:15; 
    k_vec = k_val_vec * pi;   
    K2_row = k_vec.^2;     
    n_k = length(K2_row);
    
    E_Hopf_Global    = zeros(n_D, 1); 
    E_Collapse_Global = zeros(n_D, 1); 
    
    %% 3. 核心计算
    
    p_k_off = p.k_off; p_eta = p.eta; 
    p_alpha = p.alpha; p_sigma_L = p.sigma_L; p_K = p.K;
    
    parfor i = 1:n_D
        current_D = D_vec(i);

        Ak_base_noE = 1 + p_k_off + current_D * K2_row; 
        Bk          = 1 + p_eta * K2_row;               
        
        m2_term = (-Cs * K2_row) * fixed_sigma; 
        m1_term = (-Cs * p_K * K2_row) * fixed_sigma;
        
        term_common = Ak_base_noE + p_K - p_alpha * p_sigma_L;
        
        % 系数构造: Coeff = Slope * E + Intercept
        S_A3 = zeros(1, n_k); I_A3 = Bk;
        S_A2 = K2_row;        I_A2 = m2_term + Bk .* term_common;
        S_A1 = K2_row .* term_common; I_A1 = m1_term + Ak_base_noE .* Bk * p_K;
        S_A0 = p_K * K2_row .* Ak_base_noE; I_A0 = zeros(1, n_k);
        
        % ----------------------------------------------------
        % 逐模式求解 (Mode-by-Mode Analysis)
        % ----------------------------------------------------
        max_hopf_val = 0;
        max_col_val  = 0;
        
        for k = 1:n_k
            % --- 1. 求解 Hopf 阈值 (E_h) ---
            sa1 = S_A1(k); ia1 = I_A1(k);
            sa2 = S_A2(k); ia2 = I_A2(k);
            sa0 = S_A0(k); ia3 = I_A3(k); 
            
            QA = sa1 * sa2;
            QB = sa1 * ia2 + ia1 * sa2 - sa0 * ia3;
            QC = ia1 * ia2;
            
            Delta_H = QB^2 - 4 * QA * QC;
            
            e_h = 0;
            if Delta_H >= 0
                sqrt_D = sqrt(Delta_H);
                r1 = (-QB + sqrt_D)/(2*QA);
                r2 = (-QB - sqrt_D)/(2*QA);
                candidates = [r1, r2];
                pos_r = candidates(candidates > 0);
                if ~isempty(pos_r)
                    e_h = max(pos_r); 
                end
            end
            
            % --- 2. 求解 Real Root 阈值 (E_d) ---
            p_a = ia3;
            p_b = [sa2, ia2];
            p_c = [sa1, ia1];
            p_d = [sa0, 0];
            
            % 构造判别式多项式
            bd = [sa2*sa0, ia2*sa0, 0];
            ac = [sa1*p_a, ia1*p_a];
            b2 = [sa2^2, 2*sa2*ia2, ia2^2];
            c2 = [sa1^2, 2*sa1*ia1, ia1^2];
            d2 = [sa0^2, 0, 0];
            
            T1 = 18 * conv(ac, bd); 
            b3 = conv(p_b, b2);
            T2 = -4 * conv(b3, p_d); 
            T3 = conv(b2, c2); 
            c3 = conv(p_c, c2);
            T4 = -4 * p_a * c3; 
            T5 = -27 * p_a^2 * d2; 
            
            Coeffs = zeros(1, 5);
            Coeffs = add_poly(Coeffs, T1);
            Coeffs = add_poly(Coeffs, T2);
            Coeffs = add_poly(Coeffs, T3);
            Coeffs = add_poly(Coeffs, T4);
            Coeffs = add_poly(Coeffs, T5);
            
            r_delta = roots(Coeffs);
            valid_r = r_delta(imag(r_delta)==0 & real(r_delta)>0);
            
            e_d = 0;
            if ~isempty(valid_r)
                e_d = min(valid_r); 
            end
            

            current_col_k = 0;
            if e_h > 0
                if e_d > 0
                    current_col_k = min(e_h, e_d);
                else
                    current_col_k = 0; 
                end
            end
            
            max_hopf_val = max(max_hopf_val, e_h);
            max_col_val  = max(max_col_val, current_col_k);
        end
        
        E_Hopf_Global(i)     = max_hopf_val;
        E_Collapse_Global(i) = max_col_val;
    end
    
    %% 4. 绘图
    plot_limit_E = 10;
    
    figure('Color', 'white', 'Position', [100, 100, 750, 600]);
    hold on; box on;
    
    Color_Stable   = [0.90, 0.96, 1.00]; 
    Color_Osc      = [0.60, 0.80, 0.95]; 
    Color_Collapse = [0.20, 0.40, 0.60]; 
    
    Y_fill = [D_vec; flipud(D_vec)];
    
    % 1. Stable (E > Hopf)
    X_Stable = [E_Hopf_Global; ones(size(D_vec))*plot_limit_E];
    fill(X_Stable, Y_fill, Color_Stable, 'EdgeColor', 'none', 'DisplayName', 'Stable');
    
    % 2. Oscillatory (Collapse < E < Hopf)
    X_Osc = [E_Collapse_Global; flipud(E_Hopf_Global)];
    fill(X_Osc, Y_fill, Color_Osc, 'EdgeColor', 'none', 'DisplayName', 'Oscillatory');
    
    % 3. Collapse (0 < E < Collapse)
    X_Col = [zeros(size(D_vec)); flipud(E_Collapse_Global)];
    fill(X_Col, Y_fill, Color_Collapse, 'EdgeColor', 'none', 'DisplayName', 'Collapse');
    
    % 边界线
    plot(E_Hopf_Global, D_vec, 'b--', 'LineWidth', 2, 'DisplayName', 'Hopf Boundary');
    plot(E_Collapse_Global, D_vec, 'r-', 'LineWidth', 2, 'DisplayName', 'Collapse Boundary');
    
    xlim([0, plot_limit_E]); 
    ylim([0, max(D_vec)]);
    xlabel('Elastic Modulus E', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Diffusion Coefficient D', 'FontSize', 14, 'FontWeight', 'bold');
    title({'Phase Diagram (D vs E)', ['\fontsize{10}\color[rgb]{0.3,0.3,0.3}' param_str]}, ...
          'Interpreter', 'tex', 'FontWeight', 'bold');
    
    % 标注
    text(8, 1.5, 'Stable', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0 0.2 0.4]);
    text(3.5, 1.0, 'Oscillatory', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'white');
    text(2, 1.5, 'Collapse', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'white');
    
    toc;
end

function P_out = add_poly(P_target, P_add)
    len_t = length(P_target);
    len_a = length(P_add);
    shift = len_t - len_a;
    P_out = P_target;
    P_out(shift+1:end) = P_out(shift+1:end) + P_add;
end