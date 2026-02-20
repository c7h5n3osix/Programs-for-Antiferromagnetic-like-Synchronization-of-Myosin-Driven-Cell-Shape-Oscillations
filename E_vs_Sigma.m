function phase_diagram_E_vs_Sigma_Robust_Final
%% 0. 环境初始化
    tic;
    clear;
    % close all; 
    clc;
    
    % 检查并启动并行池
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool; 
    end
    
    %% 1. 物理参数定义
    p.k_off   = 1;      
    p.eta     = 0.05;   
    p.alpha   = 0.01;     
    p.sigma_L = 30;      
    p.K       = 50;     
    
    const_D   = 1;     
    Cs = 1 / (1 + p.k_off); 
    plot_limit = 70;
    param_str = sprintf('D=%g, \\eta=%g, \\alpha=%g, k_{off}=%g, K=%g, \\sigma_L=%g', ...
                        const_D, p.eta, p.alpha, p.k_off, p.K, p.sigma_L);
    
    %% 2. 扫描参数设置
    resolution = 500; 
    
    E_vec = linspace(0.01, 20, resolution)'; 
    
    k_val_vec = 0.00:0.05:15; 
    k_vec = k_val_vec * pi;   
    K2_vec = (k_vec.^2)';     
    
    Sig_RH_min    = zeros(length(E_vec), 1); 
    Sig_Delta_max = zeros(length(E_vec), 1); 
    
    %% 3. 核心计算循环 (遍历 E)
    
    p_k_off = p.k_off; 
    p_eta = p.eta; 
    p_alpha = p.alpha; 
    p_sigma_L = p.sigma_L; 
    p_K = p.K;
    
    parfor i = 1:length(E_vec)
        current_E = E_vec(i); 
        
        % ==========================================================
        % Part A: Hopf (RH)
        % ==========================================================
        
        Ak_base_all = 1 + p_k_off + const_D * K2_vec;
        Bk_all      = 1 + p_eta * K2_vec;
        
        m2_vec = -Cs * K2_vec; 
        m1_vec = -Cs * p_K * K2_vec; 
        
        term_common = (Ak_base_all + p_K - p_alpha * p_sigma_L);
        
        c2_vec = Bk_all .* term_common + current_E * K2_vec;
        c1_vec = current_E * K2_vec .* term_common + Ak_base_all .* Bk_all * p_K;
        c0_vec = current_E * p_K * K2_vec .* Ak_base_all;
        c3_vec = Bk_all;
        
        % 计算 RH 判据
        quad_a = m1_vec .* m2_vec;
        quad_b = m1_vec .* c2_vec + m2_vec .* c1_vec;
        quad_c = c1_vec .* c2_vec - c0_vec .* c3_vec;
        
        delta_rh = quad_b.^2 - 4 * quad_a .* quad_c;
        valid_mask = delta_rh >= 0;
        
        % 向量化求根
        if any(valid_mask)
            r1_all = (-quad_b(valid_mask) + sqrt(delta_rh(valid_mask))) ./ (2 * quad_a(valid_mask));
            r2_all = (-quad_b(valid_mask) - sqrt(delta_rh(valid_mask))) ./ (2 * quad_a(valid_mask));
            positive_roots = [r1_all(r1_all > 0); r2_all(r2_all > 0)];
        else
            positive_roots = [];
        end
        
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
        
        % 准备子循环变量
        m2_sub = m2_vec; c2_sub = c2_vec;
        m1_sub = m1_vec; c1_sub = c1_vec;
        c0_sub = c0_vec; c3_sub = c3_vec;
        
        num_k = length(K2_vec);
        
        for j = 1:num_k
            % 提取标量
            m2 = m2_sub(j); c2 = c2_sub(j);
            m1 = m1_sub(j); c1 = c1_sub(j);
            c0 = c0_sub(j); c3 = c3_sub(j);
            
            pA2 = [m2, c2];
            pA1 = [m1, c1];
            
            % 卷积运算
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
            
            % 赋值操作
            len_T1 = length(T1); shift = 5 - len_T1; coeffs(shift+1:end) = coeffs(shift+1:end) + T1;
            len_T2 = length(T2); shift = 5 - len_T2; coeffs(shift+1:end) = coeffs(shift+1:end) - T2;
            len_T3 = length(T3); shift = 5 - len_T3; coeffs(shift+1:end) = coeffs(shift+1:end) + T3;
            len_T4 = length(T4); shift = 5 - len_T4; coeffs(shift+1:end) = coeffs(shift+1:end) - T4;
            coeffs(5) = coeffs(5) - T5; 
            
            r_delta = roots(coeffs);
            valid_r = r_delta(imag(r_delta)==0 & real(r_delta)>0);
            
            if ~isempty(valid_r)
                
                current_k_min_sigma = max(valid_r); 

                if current_k_min_sigma < min_delta_k
                    min_delta_k = current_k_min_sigma;
                end
            end
        end
        Sig_Delta_max(i) = min_delta_k;
    end
    
    %% 4. 绘图 (逻辑修正版)
   
    Sig_RH_min(Sig_RH_min > plot_limit) = plot_limit;
    Sig_Delta_max(Sig_Delta_max > plot_limit) = plot_limit;
    
   Sig_Collapse_Boundary = max(Sig_RH_min, Sig_Delta_max);
    
    figure('Color', 'white', 'Position', [100, 100, 750, 600]);
    hold on; box on;
    
    Color_Stable   = [0.90, 0.96, 1.00]; 
    Color_Osc      = [0.60, 0.80, 0.95]; 
    Color_Collapse = [0.20, 0.40, 0.60]; 
    
    y_fill = [E_vec; flipud(E_vec)]; 
    
    fill([zeros(size(E_vec)); flipud(Sig_RH_min)], y_fill, Color_Stable, ...
        'EdgeColor', 'none', 'DisplayName', 'Stable');

    fill([Sig_RH_min; flipud(Sig_Collapse_Boundary)], y_fill, Color_Osc, ...
        'EdgeColor', 'none', 'DisplayName', 'Oscillatory');

    fill([Sig_Collapse_Boundary; ones(size(E_vec))*plot_limit], y_fill, Color_Collapse, ...
        'EdgeColor', 'none', 'DisplayName', 'Collapse');
    
    % 4. 绘制边界线
    plot(Sig_RH_min, E_vec, 'b--', 'LineWidth', 2, 'DisplayName', 'Hopf Bifurcation');
    plot(Sig_Delta_max, E_vec, 'r-', 'LineWidth', 2, 'DisplayName', '\Delta=0 (Real Roots)');
    
    % 图像设置
    xlim([0, plot_limit]); 
    ylim([0, max(E_vec)]);
    xlabel('Active Stress \sigma_0', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Elastic Modulus E', 'FontSize', 14, 'FontWeight', 'bold');
    title({'Phase Diagram (E vs \sigma_0)', ['\fontsize{10}\color[rgb]{0.3,0.3,0.3}' param_str]}, ...
          'Interpreter', 'tex', 'FontWeight', 'bold');
    set(gca, 'Layer', 'top', 'LineWidth', 1.5, 'FontSize', 12);
    
    text(2, mean(E_vec), 'Stable', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0 0.2 0.4]);
    
    width_osc = Sig_Collapse_Boundary - Sig_RH_min;
    [max_w, idx_w] = max(width_osc);
    if max_w > 2 && Sig_RH_min(idx_w) < plot_limit
        text(Sig_RH_min(idx_w) + max_w/2, E_vec(idx_w), 'Oscillatory', ...
            'FontSize', 14, 'FontWeight', 'bold', 'Color', 'white', ...
            'HorizontalAlignment', 'center');
    end
    
    text(plot_limit*0.9, mean(E_vec), 'Collapse', 'FontSize', 14, 'FontWeight', 'bold', ...
        'Color', 'white', 'HorizontalAlignment', 'center', 'Rotation', 90);
    
    toc;
end