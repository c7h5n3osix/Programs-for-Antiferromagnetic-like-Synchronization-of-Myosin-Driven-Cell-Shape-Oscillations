function phase_diagram_D_vs_Eta_Corrected
%% 0. 环境初始化
    clear; 
    clc;
    close all;
    
    % 设置绘图风格
    set(groot, 'defaultAxesXColor', [0 0 0]);
    set(groot, 'defaultAxesYColor', [0 0 0]);
    set(groot, 'defaultAxesFontSize', 12);
    
    tic;
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool; 
    end
    
    %% 1. 物理参数定义
    p.k_off   = 1;      
    p.alpha   = 0.01;     
    p.sigma_L = 30;      
    p.K       = 50;     
    
    % --- 固定参数 (切面设置) ---
    fixed_E     = 5;      % 固定弹性模量 E
    fixed_sigma = 14;     % 固定活性应力 Sigma = 15
    
    % 导出标量
    p_Cs = 1 / (1 + p.k_off); 
    p_k_off = p.k_off; 
    p_alpha = p.alpha; 
    p_sigma_L = p.sigma_L; 
    p_K = p.K;
    
    % 绘图限制
    plot_limit_D = 2;    
    
    param_str = sprintf('E=%g, \\sigma_0=%g, \\alpha=%g, k_{off}=%g, K=%g', ...
                        fixed_E, fixed_sigma, p.alpha, p.k_off, p.K);
    
    %% 2. 扫描参数设置
    resolution = 300; 
    
    % X轴: Viscosity (Eta)
    eta_max = 1.0; % 根据您的图，关注 0-1 范围即可
    eta_vec = linspace(0.001, eta_max, resolution)'; 
    
    % 波向量 K 设置
    k_val_vec = 0.05:0.05:15; 
    k_vec = k_val_vec * pi;   
    K2_vec = (k_vec.^2); 
    
    % 结果预分配
    Boundary_D_Hopf  = zeros(length(eta_vec), 1); 
    Boundary_D_Delta = zeros(length(eta_vec), 1); 
    
    %% 3. 核心计算循环
    
    limit_val = 10000; 
    
    % 预计算与 D 无关的 Ak 部分
    % Ak = (1 + k_off) + D * k^2  => Ak_const + D * Ak_coef
    Ak_const = 1 + p_k_off;
    Ak_coef  = K2_vec;
    
    parfor i = 1:length(eta_vec)
        current_eta = eta_vec(i);
        Bk_vec = 1 + current_eta * K2_vec; % Bk 固定 (只与 eta 有关)
        
        % ==========================================================
        % Part A: Hopf Boundary (修复版 - 二次方程求解 D)
        % ==========================================================
        % 原始判据: A_q * sigma^2 + B_q * sigma + C_q = 0
        % 现在 sigma 固定，我们将其重写为关于 D 的方程:
        % Poly_A * D^2 + Poly_B * D + Poly_C = 0
        
        % 1. 将中间变量拆解为 linear form: Val = L0 + L1 * D
        
        % Term = Ak + K - alpha*sigma_L 
        Term_L0 = Ak_const + p_K - p_alpha * p_sigma_L;
        Term_L1 = Ak_coef;
        
        % c2 = Bk * Term + E * K2
        c2_L0 = Bk_vec .* Term_L0 + fixed_E * K2_vec;
        c2_L1 = Bk_vec .* Term_L1;
        
        % c1 = E * K2 * Term + Ak * Bk * K
        %    = E*K2*(T0 + T1*D) + (A0 + A1*D)*Bk*K
        c1_L0 = fixed_E * K2_vec .* Term_L0 + Ak_const .* Bk_vec * p_K;
        c1_L1 = fixed_E * K2_vec .* Term_L1 + Ak_coef  .* Bk_vec * p_K;
        
        % c0 = E * K * K2 * Ak
        c0_L0 = fixed_E * p_K * K2_vec .* Ak_const;
        c0_L1 = fixed_E * p_K * K2_vec .* Ak_coef;
        
        % c3 = Bk (Constant w.r.t D)
        c3 = Bk_vec;
        
        % m1, m2 (Constant w.r.t D)
        m2 = -p_Cs * K2_vec;
        m1 = -p_Cs * p_K * K2_vec;
        
        % 2. 代入 Hopf 判据方程 (Delta_RH = 0)
        % 判据为: Quad_B^2 - 4 * Quad_A * Quad_C = 0
        % 其中 Quad_A, Quad_B, Quad_C 是关于 sigma 的系数
        % 注意：这里 sigma 是固定的 (fixed_sigma)
        % 所以方程是: (m1*m2)*sigma^2 + (m1*c2 + m2*c1)*sigma + (c1*c2 - c0*c3) = 0
        
        Sig = fixed_sigma;
        
        % 展开项:
        % Term 1: (m1*m2) * Sig^2  -> 常数 (Part of Poly_C)
        T1 = (m1 .* m2) * Sig^2;
        
        % Term 2: (m1*c2 + m2*c1) * Sig
        %       = Sig * [ m1(c2_0 + c2_1 D) + m2(c1_0 + c1_1 D) ]
        %       = Sig*[m1 c2_0 + m2 c1_0] + D * Sig*[m1 c2_1 + m2 c1_1]
        T2_C = Sig .* (m1 .* c2_L0 + m2 .* c1_L0); % Constant part
        T2_D = Sig .* (m1 .* c2_L1 + m2 .* c1_L1); % D^1 part
        
        % Term 3: c1*c2 - c0*c3
        %       = (c1_0 + c1_1 D)(c2_0 + c2_1 D) - (c0_0 + c0_1 D)c3
        %       = [c1_0 c2_0 - c0_0 c3] + D [c1_0 c2_1 + c1_1 c2_0 - c0_1 c3] + D^2 [c1_1 c2_1]
        T3_C = c1_L0 .* c2_L0 - c0_L0 .* c3;       % Constant part
        T3_D = c1_L0 .* c2_L1 + c1_L1 .* c2_L0 - c0_L1 .* c3; % D^1 part
        T3_D2= c1_L1 .* c2_L1;                     % D^2 part
        
        % 3. 组合成 D 的二次方程系数: Poly_A D^2 + Poly_B D + Poly_C = 0
        Poly_A = T3_D2;
        Poly_B = T2_D + T3_D;
        Poly_C = T1 + T2_C + T3_C;
        
        % 4. 求解 D
        Disc_D = Poly_B.^2 - 4 * Poly_A .* Poly_C;
        valid_k = Disc_D >= 0;
        
        max_required_D = 0; % 默认 0 (如果不需 D 也能稳定)
        
        if any(valid_k)
            % 求根公式
            D_1 = (-Poly_B(valid_k) + sqrt(Disc_D(valid_k))) ./ (2 * Poly_A(valid_k));
            D_2 = (-Poly_B(valid_k) - sqrt(Disc_D(valid_k))) ./ (2 * Poly_A(valid_k));
            
            % 取正实根
            roots_all = [D_1, D_2];
            positive_roots = roots_all(roots_all > 0);
            
            if ~isempty(positive_roots)
                % 我们需要找到一个 D，使得对于 *所有* k，系统都稳定。
                % 对于每个 k，我们解出了稳定/不稳定的边界。
                % 通常需要 D 大于最大值才能压制所有 k 的不稳定性。
                max_required_D = max(positive_roots);
            end
        end
        Boundary_D_Hopf(i) = max_required_D;
        
        
        % ==========================================================
        % Part B: Delta Boundary (二分法求解 D)
        % ==========================================================
        % 保持原有逻辑，因为这部分之前是正确的（红线形状看起来对）
        low = 0; 
        high = plot_limit_D * 5; 
        
        % 快速检查
        val_low = get_critical_sigma(low, K2_vec, Ak_const, Ak_coef, Bk_vec, p_Cs, p_K, p_alpha, p_sigma_L, fixed_E);
        
        found_root = limit_val;
        
        if val_low > fixed_sigma
            found_root = 0; % D=0 时临界值已经很高，无需 D 即可避免 Collapse
        else
            % 二分
            for iter = 1:25
                mid = (low + high) / 2;
                val_mid = get_critical_sigma(mid, K2_vec, Ak_const, Ak_coef, Bk_vec, p_Cs, p_K, p_alpha, p_sigma_L, fixed_E);
                
                if val_mid > fixed_sigma
                    high = mid; % D 够大了
                else
                    low = mid;  % D 不够
                end
            end
            found_root = high;
        end
        Boundary_D_Delta(i) = found_root;
        
    end
    
    %% 4. 绘图与区域填充
    
    % 截断数据
    Boundary_D_Hopf(Boundary_D_Hopf > plot_limit_D) = plot_limit_D;
    Boundary_D_Delta(Boundary_D_Delta > plot_limit_D) = plot_limit_D;
    
    % 定义区域
    % Stable: D > max(Hopf, Delta)
    Boundary_Stable = max(Boundary_D_Hopf, Boundary_D_Delta);
    
    figure('Color', 'white', 'Position', [100, 100, 750, 600]);
    hold on; box off;
    
    Color_Stable   = [0.92, 0.96, 1.00]; 
    Color_Osc      = [0.60, 0.80, 0.95]; 
    Color_Collapse = [0.25, 0.40, 0.60]; 
    
    % 1. Stable (上方)
    fill([eta_vec; flipud(eta_vec)], [Boundary_Stable; ones(size(eta_vec))*plot_limit_D], Color_Stable, ...
        'EdgeColor', 'none', 'DisplayName', 'Stable');
        
    % 2. Oscillatory (中间: Delta < D < Hopf)
    % 仅当 Hopf > Delta 时存在振荡区
    Band_Bottom = Boundary_D_Delta;
    Band_Top    = max(Boundary_D_Delta, Boundary_D_Hopf);
    
    fill([eta_vec; flipud(eta_vec)], [Band_Bottom; flipud(Band_Top)], Color_Osc, ...
        'EdgeColor', 'none', 'DisplayName', 'Oscillatory');
        
    % 3. Collapse (下方: D < Delta)
    fill([eta_vec; flipud(eta_vec)], [zeros(size(eta_vec)); flipud(Boundary_D_Delta)], Color_Collapse, ...
        'EdgeColor', 'none', 'DisplayName', 'Collapse');
        
    % 线条
    plot(eta_vec, Boundary_D_Hopf, 'b--', 'LineWidth', 2, 'DisplayName', 'Hopf Boundary');
    plot(eta_vec, Boundary_D_Delta, 'r-', 'LineWidth', 2, 'DisplayName', '\Delta=0 Boundary');
    
    xlim([0, eta_max]); 
    ylim([0, plot_limit_D]);
    xlabel('Viscosity \eta', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Diffusion Coefficient D', 'FontSize', 14, 'FontWeight', 'bold');
    
    title({'Phase Diagram (D vs \eta)', ['\fontsize{10}\color[rgb]{0.3,0.3,0.3}' param_str]}, ...
          'Interpreter', 'tex', 'FontWeight', 'bold');
      
    legend('Location', 'northeast');
    set(gca, 'Layer', 'top', 'LineWidth', 1.5, 'FontSize', 12, 'TickDir', 'in');
    
    % 标注
    text(eta_max*0.7, plot_limit_D*0.8, 'Stable', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0 0.2 0.5]);
    text(eta_max*0.1, plot_limit_D*0.05, 'Collapse', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'white');
    
    toc;
end

%% 辅助函数 (保持不变)
function val = get_critical_sigma(D_val, K2, Ak_const, Ak_coef, Bk, Cs, pK, alpha, sigL, fixedE)
    % 重构 Ak
    Ak = Ak_const + D_val * Ak_coef;
    
    m2 = -Cs * K2; 
    
    term_common = Ak + pK - alpha * sigL;
    c2 = Bk .* term_common + fixedE * K2;
    
    m1 = -Cs * pK * K2;
    c1 = fixedE * K2 .* term_common + Ak .* Bk * pK;
    
    c0 = fixedE * pK * K2 .* Ak;
    c3 = Bk;
    
    min_sigma = 10000;
    
    for j = 1:length(K2)
        mm2 = m2(j); cc2 = c2(j);
        mm1 = m1(j); cc1 = c1(j);
        cc0 = c0(j); cc3 = c3(j);
        
        A2_sq_1 = mm2^2; A2_sq_2 = 2*mm2*cc2; A2_sq_3 = cc2^2;
        A1_sq_1 = mm1^2; A1_sq_2 = 2*mm1*cc1; A1_sq_3 = cc1^2;
        
        T3_1 = A2_sq_1 * A1_sq_1;
        T3_2 = A2_sq_1 * A1_sq_2 + A2_sq_2 * A1_sq_1;
        T3_3 = A2_sq_1 * A1_sq_3 + A2_sq_2 * A1_sq_2 + A2_sq_3 * A1_sq_1;
        T3_4 = A2_sq_2 * A1_sq_3 + A2_sq_3 * A1_sq_2;
        T3_5 = A2_sq_3 * A1_sq_3;
        
        termT2 = 4*cc0;
        T2_1 = termT2 * (mm2^3);
        T2_2 = termT2 * (3*mm2^2*cc2);
        T2_3 = termT2 * (3*mm2*cc2^2);
        T2_4 = termT2 * (cc2^3);
        
        termT4 = 4*cc3;
        T4_1 = termT4 * (mm1^3);
        T4_2 = termT4 * (3*mm1^2*cc1);
        T4_3 = termT4 * (3*mm1*cc1^2);
        T4_4 = termT4 * (cc1^3);
        
        termT1 = 18 * cc3 * cc0;
        T1_1 = termT1 * (mm2*mm1);
        T1_2 = termT1 * (mm2*cc1 + cc2*mm1);
        T1_3 = termT1 * (cc2*cc1);
        
        T5 = 27 * cc3^2 * cc0^2;
        
        C1 = T3_1;
        C2 = T3_2 - T2_1 - T4_1;
        C3 = T3_3 - T2_2 - T4_2 + T1_1;
        C4 = T3_4 - T2_3 - T4_3 + T1_2;
        C5 = T3_5 - T2_4 - T4_4 + T1_3 - T5;
        
        rr = roots([C1 C2 C3 C4 C5]);
        valid_r = rr(imag(rr)==0 & real(rr)>0);
        
        if ~isempty(valid_r)
            k_crit = max(valid_r);
            if k_crit < min_sigma
                min_sigma = k_crit;
            end
        end
    end
    val = min_sigma;
end