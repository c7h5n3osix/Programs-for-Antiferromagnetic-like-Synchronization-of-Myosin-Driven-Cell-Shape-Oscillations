function phase_diagram_RH_vs_Delta_Final
%% 0. 环境初始化
tic
    clear; 
    % close all; 
    clc;
    
    % 检查并行池
    if isempty(gcp('nocreate'))
        parpool; 
    end
    
    %% 1. 物理参数
    p.k_off   = 1;    
    p.E       = 5;    
    p.eta     = 0.05; 
    p.alpha   = 0.01;   
    p.sigma_L = 30;    
    p.K       = 50;   
    
    Cs = 1 / (1 + p.k_off); 
    plot_limit = 50;
    param_str = sprintf('E=%g, \\eta=%g, \\alpha=%g, k_{off}=%g, K=%g \\sigma_L=%g', ...
                        p.E, p.eta, p.alpha, p.k_off, p.K, p.sigma_L);
    
    %% 2. 扫描参数
    resolution = 500; 
    D_vec = linspace(0.00, 5, resolution)'; 
    k_val_vec = 0.00:0.05:15; 
    k_vec = k_val_vec * pi;   
    K2_vec = k_vec.^2; 
    
    Sig_RH_min    = zeros(length(D_vec), 1); 
    Sig_Delta_max = zeros(length(D_vec), 1); 
    
    %% 3. 核心计算循环 
    
    p_k_off = p.k_off; p_E = p.E; p_eta = p.eta; 
    p_alpha = p.alpha; p_sigma_L = p.sigma_L; p_K = p.K;
    
    parfor i = 1:length(D_vec)
        current_D = D_vec(i);
        
        % ==========================================================
        % Part A: Hopf (RH)
        % ==========================================================
        Ak_base_all = 1 + p_k_off + current_D * K2_vec;
        Bk_all      = 1 + p_eta * K2_vec;
        
        m2_vec = -Cs * K2_vec; 
        c2_vec = Bk_all .* (Ak_base_all + p_K - p_alpha * p_sigma_L) + p_E * K2_vec;
        
        m1_vec = -Cs * p_K * K2_vec;
        c1_vec = p_E * K2_vec .* (Ak_base_all + p_K - p_alpha * p_sigma_L) + Ak_base_all .* Bk_all * p_K;
        
        c0_vec = p_E * p_K * K2_vec .* Ak_base_all;
        c3_vec = Bk_all;

        quad_a = m1_vec .* m2_vec;
        quad_b = m1_vec .* c2_vec + m2_vec .* c1_vec;
        quad_c = c1_vec .* c2_vec - c0_vec .* c3_vec;
        
        delta_rh = quad_b.^2 - 4 * quad_a .* quad_c;
        valid_mask = delta_rh >= 0;
        
        r1_all = (-quad_b(valid_mask) + sqrt(delta_rh(valid_mask))) ./ (2 * quad_a(valid_mask));
        r2_all = (-quad_b(valid_mask) - sqrt(delta_rh(valid_mask))) ./ (2 * quad_a(valid_mask));
        
        % 水平拼接，防止维度错误
        positive_roots = [r1_all(r1_all > 0), r2_all(r2_all > 0)];
        
        limit_val = 1000;
        if isempty(positive_roots)
            Sig_RH_min(i) = limit_val;
        else
            Sig_RH_min(i) = min(positive_roots);
        end
        
        % ==========================================================
        % Part B: Delta (Collapse) 
        % ==========================================================
        min_delta_k = limit_val;
        
        m2_sub = m2_vec; c2_sub = c2_vec;
        m1_sub = m1_vec; c1_sub = c1_vec;
        c0_sub = c0_vec; c3_sub = c3_vec;

        for j = 1:length(K2_vec)
            m2 = m2_sub(j); c2 = c2_sub(j);
            m1 = m1_sub(j); c1 = c1_sub(j);
            c0 = c0_sub(j); c3 = c3_sub(j);
            
            pA2 = [m2, c2];
            pA1 = [m1, c1];
            
            pA2_sq = conv(pA2, pA2);    
            pA2_cb = conv(pA2_sq, pA2); 
            pA1_sq = conv(pA1, pA1);    
            pA1_cb = conv(pA1_sq, pA1); 
            
            T1 = 18 * c3 * c0 * conv(pA2, pA1); 
            T2 = 4 * c0 * pA2_cb; 
            T3 = conv(pA2_sq, pA1_sq); 
            T4 = 4 * c3 * pA1_cb; 
            T5 = 27 * c3^2 * c0^2; 
            
            coeffs = zeros(1, 5);
            coeffs(3:5) = coeffs(3:5) + T1; 
            coeffs(2:5) = coeffs(2:5) - T2; 
            coeffs(1:5) = coeffs(1:5) + T3; 
            coeffs(2:5) = coeffs(2:5) - T4; 
            coeffs(5)   = coeffs(5)   - T5; 
            
            r_delta = roots(coeffs);
            valid_r = r_delta(imag(r_delta)==0 & real(r_delta)>0);
            
            if ~isempty(valid_r)
                % 对于当前 k，取最大正根 (max)，然后在所有 k 中取最小 (min)
                min_delta_k = min(min_delta_k, max(valid_r));
            end
        end
        Sig_Delta_max(i) = min_delta_k;
    end
    
    
    Sig_RH_min(Sig_RH_min > plot_limit) = plot_limit;
    Sig_Delta_max(Sig_Delta_max > plot_limit) = plot_limit;

    %% 4. 绘图与填充
    figure('Color', 'white', 'Position', [100, 100, 700, 550]);
    hold on;
    
    Color_Stable   = [0.88, 0.94, 0.99]; 
    Color_Osc      = [0.65, 0.78, 0.92]; 
    Color_Collapse = [0.28, 0.45, 0.65]; 
    y_fill = [D_vec; flipud(D_vec)];
    
    % 1. Stable
    fill([zeros(size(D_vec)); flipud(Sig_RH_min)], y_fill, Color_Stable, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % 2. Oscillatory 
    mask_osc = Sig_RH_min < Sig_Delta_max;
    if any(mask_osc)
        D_osc = D_vec(mask_osc);
        RH_osc = Sig_RH_min(mask_osc);
        Delta_osc = Sig_Delta_max(mask_osc);
        
        fill([RH_osc; flipud(Delta_osc)], [D_osc; flipud(D_osc)], Color_Osc, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    
    % 3. Collapse
    fill([Sig_Delta_max; ones(size(D_vec))*plot_limit], y_fill, Color_Collapse, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');

    % 边界线
    plot(Sig_RH_min, D_vec, 'b--', 'LineWidth', 2, 'DisplayName', 'Hopf (RH=0)');
    plot(Sig_Delta_max, D_vec, 'r-', 'LineWidth', 2, 'DisplayName', '\Delta=0');
    
    xlim([0, plot_limit]); ylim([0, max(D_vec)]);
    xlabel('Active Stress \sigma_0', 'FontSize', 14);
    ylabel('Diffusion Coefficient D', 'FontSize', 14);
    title({'Phase Diagram ', ['\fontsize{10}' param_str]}, 'Interpreter', 'tex');
    set(gca, 'Layer', 'top', 'LineWidth', 1.5, 'FontSize', 12);
    legend('Location', 'best');
    
    % 标注
    text(4, 1.5*mean(D_vec), 'Stable', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
    [max_w, idx_w] = max(Sig_Delta_max - Sig_RH_min);
    if max_w > 5 && Sig_RH_min(idx_w) < plot_limit
        text((Sig_RH_min(idx_w)+Sig_Delta_max(idx_w))/2-15, mean(D_vec), 'Oscillatory', ...
            'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');
    end
    text(plot_limit-15, 0.5*mean(D_vec), 'Collapse', 'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', 'white', 'HorizontalAlignment', 'center');
    toc
end