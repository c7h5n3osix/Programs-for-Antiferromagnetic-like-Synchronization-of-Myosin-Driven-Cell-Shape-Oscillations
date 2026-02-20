function phase_diagram_Eta_vs_Sigma_Final_Verified
%% 0. 环境初始化
    % 清理工作区变量，防止旧变量干扰并行计算
    clear; 
    clc;
    
    % 设置全局绘图默认属性 (黑色坐标轴/字体)
    set(groot, 'defaultAxesXColor', [0 0 0]);
    set(groot, 'defaultAxesYColor', [0 0 0]);
    set(groot, 'defaultAxesZColor', [0 0 0]);
    set(groot, 'defaultTextColor', [0 0 0]);
    set(groot, 'defaultLegendTextColor', [0 0 0]);
    set(groot, 'defaultLegendEdgeColor', [0 0 0]);
    set(groot, 'defaultAxesFontSize', 12);
    
    tic;
    % 检查并启动并行池
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool; 
    end
    
    %% 1. 物理参数定义
    p.k_off   = 1;      
    p.alpha   = 0.01;     
    p.sigma_L = 30;      
    p.K       = 50;     
    
    % --- 固定参数 ---
    const_D   = 1;      % 固定扩散系数 D
    const_E   = 5;      % 固定弹性模量 E
    
    % 导出标量供并行循环使用 (减少广播开销)
    p_Cs = 1 / (1 + p.k_off); 
    p_k_off = p.k_off; 
    p_alpha = p.alpha; 
    p_sigma_L = p.sigma_L; 
    p_K = p.K;
    
    plot_limit = 70;    % X轴 (Sigma) 上限
    
    % 标题参数字符串
    param_str = sprintf('D=%g, E=%g, \\alpha=%g, k_{off}=%g, K=%g', ...
                        const_D, const_E, p.alpha, p.k_off, p.K);
    
    %% 2. 扫描参数设置
    resolution = 500; 
    
    % Y轴: Viscosity (Eta)
    eta_max = 1.0;
    eta_vec = linspace(0.001, eta_max, resolution)'; 
    
    % 波向量 K 设置 (强制为行向量)
    k_val_vec = 0.05:0.05:15; 
    k_vec = k_val_vec * pi;   
    K2_vec = (k_vec.^2); % 1xN 行向量
    
    % 结果预分配
    Sig_RH_min    = zeros(length(eta_vec), 1); 
    Sig_Delta_max = zeros(length(eta_vec), 1); 
    
    %% 3. 核心计算循环 (遍历 Eta)
    
    limit_val = 10000; % 定义大数值代表无解

    parfor i = 1:length(eta_vec)
        current_eta = eta_vec(i); 
        
        % ==========================================================
        % Part A: Hopf (RH)
        % ==========================================================
        
        % 【核心修正】将 Ak/Bk 计算移入循环内部，确保维度为 1xN
        Ak_local = 1 + p_k_off + const_D * K2_vec; 
        Bk_local = 1 + current_eta * K2_vec;
        
        % 公共项
        term_common = (Ak_local + p_K - p_alpha * p_sigma_L);
        
        % 系数计算 (全为 1xN 向量运算)
        m2_vec = -p_Cs * K2_vec; 
        c2_vec = Bk_local .* term_common + const_E * K2_vec;
        
        m1_vec = -p_Cs * p_K * K2_vec; 
        c1_vec = const_E * K2_vec .* term_common + Ak_local .* Bk_local * p_K;
        
        c0_vec = const_E * p_K * K2_vec .* Ak_local;
        c3_vec = Bk_local;
        
        % RH 判据: b^2 - 4ac
        quad_a = m1_vec .* m2_vec;
        quad_b = m1_vec .* c2_vec + m2_vec .* c1_vec;
        quad_c = c1_vec .* c2_vec - c0_vec .* c3_vec;
        
        delta_rh_val = quad_b.^2 - 4 * quad_a .* quad_c;
        valid_mask = delta_rh_val >= 0;
        
        % 求解 Hopf 边界
        if any(valid_mask)
            r1 = (-quad_b(valid_mask) + sqrt(delta_rh_val(valid_mask))) ./ (2 * quad_a(valid_mask));
            r2 = (-quad_b(valid_mask) - sqrt(delta_rh_val(valid_mask))) ./ (2 * quad_a(valid_mask));
            positive_roots = [r1(r1 > 0), r2(r2 > 0)]; 
            
            if isempty(positive_roots)
                Sig_RH_min(i) = limit_val;
            else
                Sig_RH_min(i) = min(positive_roots);
            end
        else
            Sig_RH_min(i) = limit_val;
        end
        
        % ==========================================================
        % Part B: Delta (Collapse)
        % ==========================================================
        min_delta_k = limit_val;
        
        % 局部变量切片，避免重复索引
        m2_sub = m2_vec; c2_sub = c2_vec;
        m1_sub = m1_vec; c1_sub = c1_vec;
        c0_sub = c0_vec; c3_sub = c3_vec;
        
        num_k = length(K2_vec);
        
        for j = 1:num_k
            m2 = m2_sub(j); c2 = c2_sub(j);
            m1 = m1_sub(j); c1 = c1_sub(j);
            c0 = c0_sub(j); c3 = c3_sub(j);
            
            pA2 = [m2, c2];
            pA1 = [m1, c1];
            
            % 卷积运算 (Polynomial multiplication)
            pA2_sq = conv(pA2, pA2);    
            pA2_cb = conv(pA2_sq, pA2); 
            pA1_sq = conv(pA1, pA1);    
            pA1_cb = conv(pA1_sq, pA1); 
            
            % 判别式各项 (注意 conv 长度不同，需手动对齐)
            T1 = 18 * c3 * c0 * conv(pA2, pA1); % Length 3
            T2 = 4 * c0 * pA2_cb;               % Length 4
            T3 = conv(pA2_sq, pA1_sq);          % Length 5
            T4 = 4 * c3 * pA1_cb;               % Length 4
            T5 = 27 * c3^2 * c0^2;              % Scalar
            
            % 维度对齐到 1x5
            coeffs = zeros(1, 5);
            coeffs(1:5) = coeffs(1:5) + T3;
            coeffs(2:5) = coeffs(2:5) - T2 - T4;
            coeffs(3:5) = coeffs(3:5) + T1;
            coeffs(5)   = coeffs(5)   - T5; 
            
            r_delta = roots(coeffs);
            valid_r = r_delta(imag(r_delta)==0 & real(r_delta)>0);
            
            if ~isempty(valid_r)
                current_k_critical = max(valid_r); 
                if current_k_critical < min_delta_k
                    min_delta_k = current_k_critical;
                end
            end
        end
        Sig_Delta_max(i) = min_delta_k;
    end
    
    %% 4. 绘图 (Robust Version)
    
    % 1. 数据截断 (防止 infinity 破坏绘图)
    Sig_RH_min(Sig_RH_min > plot_limit) = plot_limit;
    Sig_Delta_max(Sig_Delta_max > plot_limit) = plot_limit;
    
    % 2. 定义绘图边界 (避免 fill 出现交叉多边形)
    % Stable 区右边界：取 RH 和 Delta 的较小值
    Boundary_Stable_Right = min(Sig_RH_min, Sig_Delta_max);
    
    % Oscillatory 区左边界：就是 Stable 的右边界
    % Oscillatory 区右边界：Delta
    % 如果 RH > Delta (无振荡)，则 左边界=Delta, 右边界=Delta，宽度为0 (自动隐藏)
    Boundary_Osc_Left = Boundary_Stable_Right;
    Boundary_Osc_Right = Sig_Delta_max;
    
    figure('Color', 'white', 'Position', [100, 100, 750, 600]);
    hold on; box off; % 使用 box off 风格
    
    % 3. 配色
    Color_Stable   = [0.92, 0.96, 1.00]; 
    Color_Osc      = [0.60, 0.80, 0.95]; 
    Color_Collapse = [0.25, 0.40, 0.60]; 
    
    y_fill = [eta_vec; flipud(eta_vec)]; 
    
    % 4. 执行填充
    % Stable
    fill([zeros(size(eta_vec)); flipud(Boundary_Stable_Right)], y_fill, Color_Stable, ...
        'EdgeColor', 'none', 'DisplayName', 'Stable');
    
    % Oscillatory (使用修正后的边界，防止伪影)
    fill([Boundary_Osc_Left; flipud(Boundary_Osc_Right)], y_fill, Color_Osc, ...
        'EdgeColor', 'none', 'DisplayName', 'Oscillatory');
    
    % Collapse
    fill([Sig_Delta_max; ones(size(eta_vec))*plot_limit], y_fill, Color_Collapse, ...
        'EdgeColor', 'none', 'DisplayName', 'Collapse');
    
    % 5. 绘制线条
    plot(Sig_RH_min, eta_vec, 'b--', 'LineWidth', 2, 'DisplayName', 'Hopf Bifurcation');
    plot(Sig_Delta_max, eta_vec, 'r-', 'LineWidth', 2, 'DisplayName', '\Delta=0');
    
    % 6. 装饰
    xlim([0, plot_limit]); 
    ylim([0, max(eta_vec)]);
    xlabel('Active Stress \sigma_0', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Viscosity \eta', 'FontSize', 14, 'FontWeight', 'bold');
    
    title({'Phase Diagram (\eta vs \sigma_0)', ['\fontsize{10}\color[rgb]{0.3,0.3,0.3}' param_str]}, ...
          'Interpreter', 'tex', 'FontWeight', 'bold');
    
    set(gca, 'Layer', 'top', 'LineWidth', 1.5, 'FontSize', 12, 'TickDir', 'in');
    
    legend('Location', 'northeast');
    
    % 7. 文本标注 (确保在图层最上方)
    uistack(findobj(gca, 'Type', 'text'), 'top');
    
    toc;
end